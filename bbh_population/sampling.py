"""GWTC-4.0 BP2P source-catalog sampling.

This module draws whole, correlated hyperparameter samples from the official
GWTC-4.0 ``popsummary`` file for the model
``BrokenPowerLawTwoPeaks + GaussianComponentSpins + PowerLawRedshift``.
It produces source catalogs only: detector selection is deliberately absent.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Tuple, Union

import h5py
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.cosmology import Planck18
from scipy.stats import truncnorm


OFFICIAL_MODEL_FILE = (
    "BBHMassSpinRedshift_BrokenPowerLawTwoPeaks_"
    "GaussianComponentSpins_PowerLawRedshift.h5"
)


@dataclass(frozen=True)
class BP2PParameters:
    """One joint GWTC-4.0 posterior sample in paper notation."""

    alpha_1: float
    alpha_2: float
    m_break: float
    delta_m1: float
    delta_m2: float
    lambda_0: float
    lambda_1: float
    kappa: float
    m1_low: float
    m2_low: float
    m_high: float
    mu_1: float
    mu_2: float
    sigma_1: float
    sigma_2: float
    beta_q: float
    rate_0: float
    mu_chi: float
    sigma_chi: float
    mu_tilt: float
    sigma_tilt: float
    zeta_tilt: float
    posterior_index: int


# Official HDF5 hyperparameter names -> clear catalog names.
HDF_PARAMETER_MAP = {
    "alpha_1": "alpha_1",
    "alpha_2": "alpha_2",
    "break_mass": "m_break",
    "delta_m_1": "delta_m1",
    "delta_m_2": "delta_m2",
    "lam_0": "lambda_0",
    "lam_1": "lambda_1",
    "lamb": "kappa",
    "mlow_1": "m1_low",
    "mlow_2": "m2_low",
    "mmax": "m_high",
    "mpp_1": "mu_1",
    "mpp_2": "mu_2",
    "sigpp_1": "sigma_1",
    "sigpp_2": "sigma_2",
    "beta": "beta_q",
    "rate": "rate_0",
    "mu_chi": "mu_chi",
    "sigma_chi": "sigma_chi",
    "mu_spin": "mu_tilt",
    "sigma_spin": "sigma_tilt",
    "xi_spin": "zeta_tilt",
}


def _decode(values: np.ndarray) -> Tuple[str, ...]:
    return tuple(value.decode() if isinstance(value, bytes) else str(value) for value in values)


def load_official_posterior(path: Union[str, Path]) -> Tuple[BP2PParameters, ...]:
    """Read the official joint posterior and retain paper-valid samples only."""
    path = Path(path)
    if path.name != OFFICIAL_MODEL_FILE:
        raise ValueError(f"expected the official default-model file {OFFICIAL_MODEL_FILE}")

    with h5py.File(path, "r") as handle:
        model_names = " ".join(_decode(handle.attrs["model_names"]))
        if "TwoPeakBrokenPowerLawSmoothedMassDistribution" not in model_names:
            raise ValueError("HDF5 file is not the GWTC-4.0 BP2P mass model")
        if "iid_spin_magnitude_gaussian" not in model_names:
            raise ValueError("HDF5 file is not the Gaussian Component Spins model")

        names = _decode(handle.attrs["hyperparameters"])
        samples = handle["posterior/hyperparameter_samples"][:]

    missing = set(HDF_PARAMETER_MAP) - set(names)
    if missing:
        raise ValueError(f"official posterior is missing fields: {sorted(missing)}")

    columns: Dict[str, np.ndarray] = {
        output_name: samples[:, names.index(input_name)]
        for input_name, output_name in HDF_PARAMETER_MAP.items()
    }
    variance = samples[:, names.index("variance")]
    finite = np.isfinite(variance) & (variance < 1.0)
    for values in columns.values():
        finite &= np.isfinite(values)

    return tuple(
        BP2PParameters(
            **{name: float(values[index]) for name, values in columns.items()},
            posterior_index=int(index),
        )
        for index in np.flatnonzero(finite)
    )


def planck_taper(mass: np.ndarray, low: float, width: float) -> np.ndarray:
    """The GWTC Planck taper S(m | low, width), Eq. (B21)."""
    mass = np.asarray(mass, dtype=float)
    if width <= 0:
        return (mass >= low).astype(float)
    result = np.zeros_like(mass)
    high = low + width
    above = mass >= high
    middle = (mass >= low) & ~above
    result[above] = 1.0
    x = mass[middle] - low
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        exponent = width / x + width / (x - width)
        result[middle] = 1.0 / (1.0 + np.exp(exponent))
    return result


def _inverse_cdf(grid: np.ndarray, density: np.ndarray, size: int, rng: np.random.Generator) -> np.ndarray:
    density = np.clip(np.asarray(density, dtype=float), 0.0, None)
    cdf = np.concatenate(([0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(grid))))
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("distribution has zero or non-finite normalization")
    cdf /= cdf[-1]
    return np.interp(rng.random(size), cdf, grid)


def sample_primary_mass(params: BP2PParameters, size: int, rng: np.random.Generator) -> np.ndarray:
    """Sample the BP2P primary-mass model in Eq. (B20)."""
    if not (params.m1_low < params.m_break < params.m_high):
        raise ValueError("invalid posterior sample: m1_low < m_break < m_high is required")
    if min(params.lambda_0, params.lambda_1, 1.0 - params.lambda_0 - params.lambda_1) < 0:
        raise ValueError("invalid posterior sample: mixture weights must be non-negative")

    grid = np.linspace(params.m1_low, params.m_high, 20_000)
    broken = np.where(
        grid < params.m_break,
        (grid / params.m_break) ** (-params.alpha_1),
        (grid / params.m_break) ** (-params.alpha_2),
    )
    peak_1 = np.exp(-0.5 * ((grid - params.mu_1) / params.sigma_1) ** 2)
    peak_2 = np.exp(-0.5 * ((grid - params.mu_2) / params.sigma_2) ** 2)

    def normalized(component: np.ndarray) -> np.ndarray:
        return component / np.trapz(component, grid)

    density = (
        params.lambda_0 * normalized(broken)
        + params.lambda_1 * normalized(peak_1)
        + (1.0 - params.lambda_0 - params.lambda_1) * normalized(peak_2)
    ) * planck_taper(grid, params.m1_low, params.delta_m1)
    return _inverse_cdf(grid, density, size, rng)


def _sample_power_law(lower: np.ndarray, exponent: float, rng: np.random.Generator) -> np.ndarray:
    """Draw q**exponent on element-wise intervals [lower, 1]."""
    if np.any(lower <= 0) or np.any(lower > 1):
        raise ValueError("invalid mass-ratio bounds")
    draw = rng.random(lower.size)
    power = exponent + 1.0
    if np.isclose(power, 0.0):
        return lower ** (1.0 - draw)
    return (lower**power + draw * (1.0 - lower**power)) ** (1.0 / power)


def sample_mass_ratio(m1: np.ndarray, params: BP2PParameters, rng: np.random.Generator) -> np.ndarray:
    """Sample p(q|m1) from Eq. (B22) through tabulated inverse CDFs.

    The paper normalizes this conditional distribution on a log-uniform
    primary-mass grid.  Direct CDF inversion avoids the arbitrarily long
    rejection-sampling tail caused by broad secondary-mass tapers.
    """
    m1 = np.asarray(m1, dtype=float)
    if m1.size == 0:
        return np.empty(0, dtype=float)
    grid_low = max(params.m1_low, params.m2_low * (1.0 + 1e-9), float(m1.min()))
    grid_high = min(params.m_high, float(m1.max()))
    if grid_high <= grid_low:
        grid_high = grid_low * (1.0 + 1e-8)
    m1_grid = np.geomspace(grid_low, grid_high, 1024)
    q_lower = params.m2_low / m1_grid
    unit_grid = np.linspace(0.0, 1.0, 2048)
    q_grid = q_lower[:, None] + (1.0 - q_lower[:, None]) * unit_grid[None, :]
    secondary_mass = m1_grid[:, None] * q_grid
    if params.delta_m2 <= 0:
        log_taper = np.where(secondary_mass >= params.m2_low, 0.0, -np.inf)
    else:
        x = secondary_mass - params.m2_low
        log_taper = np.full_like(secondary_mass, -np.inf)
        above = x >= params.delta_m2
        middle = (x >= 0.0) & ~above
        log_taper[above] = 0.0
        with np.errstate(divide="ignore", invalid="ignore"):
            exponent = params.delta_m2 / x[middle] + params.delta_m2 / (x[middle] - params.delta_m2)
        log_taper[middle] = -np.logaddexp(0.0, exponent)
    log_density = params.beta_q * np.log(q_grid) + log_taper
    row_max = np.max(log_density, axis=1, keepdims=True)
    density = np.exp(log_density - row_max)
    cdf = np.zeros_like(q_grid)
    cdf[:, 1:] = np.cumsum(
        0.5 * (density[:, 1:] + density[:, :-1]) * np.diff(q_grid, axis=1), axis=1
    )
    normalization = cdf[:, -1]
    if np.any(~np.isfinite(normalization)) or np.any(normalization <= 0):
        raise ValueError("non-finite mass-ratio normalization")
    cdf /= normalization[:, None]

    insertion = np.searchsorted(m1_grid, m1)
    insertion = np.clip(insertion, 1, len(m1_grid) - 1)
    left = insertion - 1
    nearest = np.where(
        np.abs(np.log(m1) - np.log(m1_grid[left]))
        <= np.abs(np.log(m1_grid[insertion]) - np.log(m1)),
        left,
        insertion,
    )
    q = np.empty(m1.size)
    for index in np.unique(nearest):
        mask = nearest == index
        q[mask] = np.interp(rng.random(mask.sum()), cdf[index], q_grid[index])
    return q


def _redshift_density(params: BP2PParameters, z_max: float) -> Tuple[np.ndarray, np.ndarray]:
    if z_max <= 0:
        raise ValueError("z_max must be positive")
    z = np.linspace(0.0, z_max, 20_000)
    dvc_dz = (4.0 * np.pi * u.sr * Planck18.differential_comoving_volume(z)).to_value(u.Gpc**3)
    return z, dvc_dz * (1.0 + z) ** (params.kappa - 1.0)


def expected_source_count(params: BP2PParameters, observation_time_yr: float, z_max: float) -> float:
    """Detector-frame source count from Eq. (B25), drawn as a Poisson mean."""
    z, density = _redshift_density(params, z_max)
    return params.rate_0 * observation_time_yr * np.trapz(density, z)


def sample_redshift(params: BP2PParameters, size: int, z_max: float, rng: np.random.Generator) -> np.ndarray:
    z, density = _redshift_density(params, z_max)
    return _inverse_cdf(z, density, size, rng)


def _sample_truncated_normal(
    mean: float, sigma: float, low: float, high: float, size: int, rng: np.random.Generator
) -> np.ndarray:
    if sigma <= 0:
        raise ValueError("truncated-normal width must be positive")
    return truncnorm.rvs(
        (low - mean) / sigma,
        (high - mean) / sigma,
        loc=mean,
        scale=sigma,
        size=size,
        random_state=rng,
    )


def sample_component_spins(
    params: BP2PParameters, size: int, rng: np.random.Generator
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sample Eq. (B26) and the correlated tilt mixture in Eq. (B27)."""
    a_1 = _sample_truncated_normal(params.mu_chi, params.sigma_chi, 0.0, 1.0, size, rng)
    a_2 = _sample_truncated_normal(params.mu_chi, params.sigma_chi, 0.0, 1.0, size, rng)
    gaussian_component = rng.random(size) < params.zeta_tilt
    cos_tilt_1 = rng.uniform(-1.0, 1.0, size)
    cos_tilt_2 = rng.uniform(-1.0, 1.0, size)
    n_gaussian = gaussian_component.sum()
    cos_tilt_1[gaussian_component] = _sample_truncated_normal(
        params.mu_tilt, params.sigma_tilt, -1.0, 1.0, n_gaussian, rng
    )
    cos_tilt_2[gaussian_component] = _sample_truncated_normal(
        params.mu_tilt, params.sigma_tilt, -1.0, 1.0, n_gaussian, rng
    )
    phi_1 = rng.uniform(0.0, 2.0 * np.pi, size)
    phi_2 = rng.uniform(0.0, 2.0 * np.pi, size)
    return a_1, a_2, cos_tilt_1, cos_tilt_2, phi_1, phi_2


def generate_catalog(
    params: BP2PParameters, observation_time_yr: float, z_max: float, rng: np.random.Generator
) -> pd.DataFrame:
    """Generate one source catalog from one complete official posterior sample."""
    expected = expected_source_count(params, observation_time_yr, z_max)
    count = rng.poisson(expected)
    m1 = sample_primary_mass(params, count, rng)
    q = sample_mass_ratio(m1, params, rng)
    a_1, a_2, cos_tilt_1, cos_tilt_2, phi_1, phi_2 = sample_component_spins(params, count, rng)
    redshift = sample_redshift(params, count, z_max, rng)
    catalog = pd.DataFrame(
        {
            "z": redshift,
            "m1_Msun": m1,
            "m2_Msun": m1 * q,
            "q": q,
            "a_1": a_1,
            "a_2": a_2,
            "cos_tilt_1": cos_tilt_1,
            "cos_tilt_2": cos_tilt_2,
            "phi_1_rad": phi_1,
            "phi_2_rad": phi_2,
            "t_c_detector_yr": rng.uniform(0.0, observation_time_yr, count),
            "ra_rad": rng.uniform(0.0, 2.0 * np.pi, count),
            "sin_dec": rng.uniform(-1.0, 1.0, count),
            "cos_iota": rng.uniform(-1.0, 1.0, count),
            "psi_rad": rng.uniform(0.0, np.pi, count),
        }
    )
    catalog.attrs.update(asdict(params))
    catalog.attrs["expected_source_count"] = expected
    return catalog

