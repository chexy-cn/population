"""Optional regressions against the historical manuscript files (never modified)."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from bbh_population.extension import ExtensionParameters, high_z_catalog, rate_integrals
from bbh_population.pipeline import ensemble_draws
from bbh_population.sampling import OFFICIAL_MODEL_FILE, generate_catalog, load_official_posterior

POSTERIOR = Path("o4a-astro/data_release") / OFFICIAL_MODEL_FILE
ORIGINAL = Path("bbh_catalogs_gwtc4_bp2p")
HIGH = Path("catalog_gwtc4_highz_extension_z1p5_30")


@unittest.skipUnless(POSTERIOR.is_file() and (ORIGINAL / "catalog_001.csv").is_file()
                     and (HIGH / "catalog_001.csv").is_file(), "historical files are local-only")
class ReferenceTests(unittest.TestCase):
    def test_all_500_posterior_draws_and_high_catalog_001_bytes(self):
        posterior = load_official_posterior(POSTERIOR)
        draws = ensemble_draws(posterior, 500, 20260918)
        historical = pd.read_csv(ORIGINAL / "catalog_parameters_summary.csv")
        self.assertEqual([p.posterior_index for p, _ in draws], historical.posterior_index.tolist())
        params, seed = draws[0]
        frame = generate_catalog(params, 1, 1.5, np.random.default_rng(seed))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "regression.csv"
            frame.to_csv(path, index=False)
            original = pd.read_csv(ORIGINAL / "catalog_001.csv", float_precision="round_trip")
            # The legacy low catalog predates the current q inverse-CDF implementation.
            # Its full RNG history is not frozen: regenerated low rows are NOT byte-identical.
            np.testing.assert_array_equal(frame.m1_Msun, original.m1_Msun)
            config = ExtensionParameters()
            low, high = rate_integrals(config)
            rng = np.random.default_rng(20260920)
            seed = int(rng.integers(0, np.iinfo(np.uint64).max, size=500, dtype=np.uint64)[0])
            rng = np.random.default_rng(seed)
            count = int(rng.poisson(len(original) * high / low))
            extra = high_z_catalog(params, count, config, 1, rng)
            extra.to_csv(path, index=False)
            self.assertEqual(path.read_bytes(), (HIGH / "catalog_001.csv").read_bytes())
