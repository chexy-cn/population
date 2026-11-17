"""Reference rate and high-redshift sampler in manuscript Appendix A."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.cosmology import Planck18

from .sampling import (
    BP2PParameters, _inverse_cdf, sample_component_spins,
    sample_mass_ratio, sample_primary_mass,
)


@dataclass(frozen=True)
class ExtensionParameters:
    """Reference evolution, distinct from posterior kappa_j."""

    z_join: float = 1.5
    z_max: float = 30.0
    gamma: float = 3.2
    z_turn: float = 1.9
    kappa_high: float = 5.6
    grid_size: int = 30_000

    def __post_init__(self):
        values = (self.z_join, self.z_max, self.gamma, self.z_turn, self.kappa_high)
        if not np.isfinite(values).all():
            raise ValueError("extension parameters must be finite")
        if not (0 < self.z_join < self.z_max and self.z_turn > -1):
            raise ValueError("require 0 < z_join < z_max and z_turn > -1")
        if self.gamma <= 0 or self.kappa_high <= 0 or self.grid_size < 100:
            raise ValueError("positive slopes and grid_size >= 100 are required")


def md_shape(z, config: ExtensionParameters):
    """g(z), with high-z asymptotic slope -kappa_high."""
    z = np.asarray(z, dtype=float)
    return (1 + z) ** config.gamma / (
        1 + ((1 + z) / (1 + config.z_turn)) ** (config.gamma + config.kappa_high)
    )


def high_rate_shape(z, config: ExtensionParameters):
    """High branch matched to the reference low branch at z_join."""
    return (1 + config.z_join) ** config.gamma * md_shape(z, config) / md_shape(config.z_join, config)


def reference_rate_shape(z, config: ExtensionParameters):
    """Reference R/R0 curve; actual retained low rows use posterior kappa_j."""
    z = np.asarray(z, dtype=float)
    return np.where(z <= config.z_join, (1 + z) ** config.gamma, high_rate_shape(z, config))


def event_density(z, config: ExtensionParameters):
    """Full-sky observer-frame shape in Gpc^3, excluding R0 and duration."""
    dvc_dz = (4 * np.pi * u.sr * Planck18.differential_comoving_volume(z)).to_value(u.Gpc**3)
    return reference_rate_shape(z, config) * dvc_dz / (1 + z)


def integration_grid(z_min, z_max, size):
    """Historical log(1+z) grid retained for numerical reproducibility."""
    return np.expm1(np.linspace(np.log1p(z_min), np.log1p(z_max), size))


def rate_integrals(config: ExtensionParameters):
    """Return (I_low, I_high), with R0 factored out, in Gpc^3."""
    low = integration_grid(0, config.z_join, config.grid_size)
    high = integration_grid(config.z_join, config.z_max, config.grid_size)
    return (float(np.trapz(event_density(low, config), low)),
            float(np.trapz(event_density(high, config), high)))


def sample_high_z(size, config, rng):
    grid = integration_grid(config.z_join, config.z_max, config.grid_size)
    result = _inverse_cdf(grid, event_density(grid, config), size, rng)
    return np.maximum(result, np.nextafter(config.z_join, np.inf))


def high_z_catalog(params: BP2PParameters, size, config, years, rng) -> pd.DataFrame:
    """Draw new sources in the original RNG order and 15-column schema."""
    m1 = sample_primary_mass(params, size, rng)
    q = sample_mass_ratio(m1, params, rng)
    a1, a2, cos1, cos2, phi1, phi2 = sample_component_spins(params, size, rng)
    return pd.DataFrame({
        "z": sample_high_z(size, config, rng),
        "m1_Msun": m1, "m2_Msun": m1 * q, "q": q,
        "a_1": a1, "a_2": a2, "cos_tilt_1": cos1, "cos_tilt_2": cos2,
        "phi_1_rad": phi1, "phi_2_rad": phi2,
        "t_c_detector_yr": rng.uniform(0, years, size),
        "ra_rad": rng.uniform(0, 2 * np.pi, size),
        "sin_dec": rng.uniform(-1, 1, size),
        "cos_iota": rng.uniform(-1, 1, size),
        "psi_rad": rng.uniform(0, np.pi, size),
    })
