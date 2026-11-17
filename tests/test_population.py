"""Scientific and reproducibility checks using small synthetic inputs."""

import json
import tempfile
import unittest
from dataclasses import asdict
from pathlib import Path

import h5py
import numpy as np

from bbh_population.extension import ExtensionParameters, high_rate_shape, rate_integrals, sample_high_z
from bbh_population.pipeline import ensemble_draws, generate_ensemble
from bbh_population.sampling import BP2PParameters, HDF_PARAMETER_MAP, OFFICIAL_MODEL_FILE, load_official_posterior
from bbh_population.validation import compare_append, validate_run


def synthetic_posterior(directory):
    """Write a tiny joint posterior for pipeline tests, never as scientific data."""
    params = BP2PParameters(
        alpha_1=1.7, alpha_2=4.5, m_break=34, delta_m1=2, delta_m2=2,
        lambda_0=.4, lambda_1=.4, kappa=3.1, m1_low=5, m2_low=3,
        m_high=100, mu_1=10, mu_2=33, sigma_1=1, sigma_2=4,
        beta_q=1.2, rate_0=.05, mu_chi=.1, sigma_chi=.3,
        mu_tilt=.4, sigma_tilt=1, zeta_tilt=.6, posterior_index=0,
    )
    names = list(HDF_PARAMETER_MAP) + ["variance"]
    values = [getattr(params, HDF_PARAMETER_MAP[name]) for name in names[:-1]] + [.1]
    samples = np.tile(values, (4, 1))
    samples[1, names.index("lamb")] = 3.5
    samples[2, names.index("variance")] = 2  # rejected by historical variance cut
    samples[3, names.index("alpha_1")] = np.nan
    path = Path(directory) / OFFICIAL_MODEL_FILE
    with h5py.File(path, "w") as handle:
        handle.attrs["model_names"] = np.asarray([
            "TwoPeakBrokenPowerLawSmoothedMassDistribution", "iid_spin_magnitude_gaussian",
        ], dtype="S")
        handle.attrs["hyperparameters"] = np.asarray(names, dtype="S")
        handle.create_dataset("posterior/hyperparameter_samples", data=samples)
    return path


class PopulationTests(unittest.TestCase):
    def test_joint_draw_filter_and_distinct_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            rows = load_official_posterior(synthetic_posterior(directory))
            self.assertEqual([p.posterior_index for p in rows], [0, 1])
            draws = ensemble_draws(rows, 2, 123)
            self.assertEqual(len({p.posterior_index for p, _ in draws}), 2)
            self.assertEqual(draws, ensemble_draws(rows, 2, 123))

    def test_reference_continuity_and_integral(self):
        config = ExtensionParameters()
        self.assertAlmostEqual(float(high_rate_shape(config.z_join, config)), 2.5 ** 3.2, places=13)
        low, high = rate_integrals(config)
        self.assertAlmostEqual(high / low, 2.251671096, places=7)
        low50, high50 = rate_integrals(ExtensionParameters(z_max=50))
        self.assertLess(abs(high50 / low50 - high / low), 5e-7)
        draws = sample_high_z(10_000, config, np.random.default_rng(23))
        self.assertTrue(((draws > 1.5) & (draws <= 30)).all())

    def test_pipeline_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = synthetic_posterior(root)
            settings = json.loads(Path("configs/paper_append.json").read_text())
            settings["catalogs"] = 2
            output = root / "paper"
            generate_ensemble(settings, path, output, workers=1)
            report = validate_run(output)
            self.assertTrue(report["passed"])
            self.assertEqual(report["catalog_count"], 2)
            parallel = root / "paper-parallel"
            generate_ensemble(settings, path, parallel, workers=2)
            for component in ("lowz", "highz", "combined"):
                for catalog_id in (1, 2):
                    name = f"catalog_{catalog_id:03d}.csv"
                    self.assertEqual((output / component / name).read_bytes(),
                                     (parallel / component / name).read_bytes())
            with self.assertRaises(FileExistsError):
                generate_ensemble(settings, path, output)
            # Detect mutations in the appended suffix, not only a correct low prefix.
            combined = output / "combined/catalog_001.csv"
            with combined.open("ab") as handle:
                handle.write(b"unexpected\n")
            self.assertFalse(compare_append(output / "lowz/catalog_001.csv",
                                            output / "highz/catalog_001.csv", combined))

    def test_gwtc4_and_empty_realization(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = synthetic_posterior(root)
            settings = json.loads(Path("configs/gwtc4.json").read_text())
            settings.update(catalogs=1, years=1e-15)
            output = root / "empty"
            generate_ensemble(settings, path, output)
            report = validate_run(output)
            self.assertTrue(report["passed"])
            self.assertEqual(report["N_total"], 0)


if __name__ == "__main__":
    unittest.main()
