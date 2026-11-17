# GWTC-4 BBH population catalogs

Reproducible **intrinsic source-population models** for stellar-mass binary black
holes, with two explicit configurations:

| Configuration | Prescription | Redshift range |
| --- | --- | --- |
| `configs/gwtc4.json` | Official GWTC-4.0 BP2P joint posterior, including each draw's local rate and redshift exponent | 0–1.5 |
| `configs/paper_append.json` | Same low-z model plus the Appendix A high-z extension used in *Constraining the dynamical Chern-Simons gravity with future gravitational-wave detectors* | 0–30 |

This repository covers population generation. Waveform, SNR, detector selection,
and Fisher-matrix calculations belong to the downstream analysis and are not
implemented here. This is an implementation of the named GWTC-4 BP2P model,
not every alternative population model investigated in GWTC-4.

## Installation

Python 3.9 or later:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

For the historical numerical stack, install `requirements-reproduction.txt`
before installing this project. The existing Conda `GW` environment can be used.
No personal paths, Conda environment name, or SLURM partition is required by the code.

## Official input

Obtain the GWTC-4.0 population release from
[Zenodo](https://doi.org/10.5281/zenodo.16911563), linked by the
[official GWOSC release](https://gwosc.org/GWTC-4.0/). Extract this file:

```text
BBHMassSpinRedshift_BrokenPowerLawTwoPeaks_GaussianComponentSpins_PowerLawRedshift.h5
```

Keep its name unchanged. Store it anywhere and pass `--posterior`. The 342 MB
HDF5 input is not redistributed in this source release. Its SHA-256, filtering
rules, and scientific references are in [docs/INPUTS.md](docs/INPUTS.md).

## Generate the GWTC-4 low-z model

```bash
python -m bbh_population generate \
  --config configs/gwtc4.json \
  --posterior /path/to/BBHMassSpinRedshift_BrokenPowerLawTwoPeaks_GaussianComponentSpins_PowerLawRedshift.h5 \
  --output-dir outputs/gwtc4 --workers 8
```

The preset selects 500 distinct joint hyperparameter rows without replacement.
Each catalog has its own posterior parameters and source RNG. Its event count
is drawn from the posterior local rate integrated over one observer year.

## Generate the manuscript model

From scratch, generate low-z sources and then append high-z sources:

```bash
python -m bbh_population generate \
  --config configs/paper_append.json \
  --posterior /path/to/BBHMassSpinRedshift_BrokenPowerLawTwoPeaks_GaussianComponentSpins_PowerLawRedshift.h5 \
  --output-dir outputs/paper --workers 8
```

To preserve the actual historical low-z realization used in the manuscript,
add `--original-dir /path/to/bbh_catalogs_gwtc4_bp2p`. That directory must include
the numbered CSVs and `catalog_parameters_summary.csv`. Every original source
byte is copied first; new high-z rows follow. Original files are read-only inputs.

**Reproducibility boundary:** the historical low-z files predate the currently
available conditional mass-ratio implementation. A fresh run reproduces the
specified population model, but does not reconstruct every historical random
source row. Reusing the historical originals provides exact row preservation;
the high-z generator has a byte-level regression against the manuscript's
`catalog_001`. See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

Change ensemble size with `--catalogs 100`; change independent realizations with
`--seed` and `--extension-seed`. Ensemble size affects posterior selection and
child RNG assignment. To reproduce only catalog 1 of the **500-member** paper
ensemble, use `--catalog-ids 1` without changing `--catalogs`.

Output directories must be new. Failed runs retain their partial output and
`run.json` status; select a new directory for a retry. There is no implicit resume.

## Outputs and validation

```text
outputs/paper/
  run.json                              # settings, versions, input hash and status
  lowz/                                 # omitted when existing originals are reused
    catalog_001.csv ...
    catalog_parameters_summary.csv
  highz/
    catalog_001.csv ...
  combined/
    catalog_001.csv ...
    extension_parameters_summary.csv    # posterior IDs, child seeds, counts and hashes
```

Validate **all rows and all 15 numeric columns**, plus exact low-z prefix and
high-z suffix bytes:

```bash
python -m bbh_population validate --run-dir outputs/paper --workers 8
```

Validation writes `validation/summary.json` and exits nonzero on failed checks.
For reference-rate and actual catalog histogram plots:

```bash
python -m pip install -e '.[plots]'
python -m bbh_population plot --run-dir outputs/paper
```

The source masses, spin parameters, sky position, orientation and time fields are
documented in [docs/MODELS.md](docs/MODELS.md). Luminosity distance, detector-frame
masses and aligned spin projections are derivable without redefining stored fields.

## Tests and cluster use

```bash
python -m unittest discover -s tests -v
```

Tests use synthetic posterior rows. An additional regression runs automatically
when the historical posterior and catalog directories are present locally.
See [examples/generate.slurm](examples/generate.slurm) for a portable single-node
SLURM job that respects the allocated CPUs and caps numerical-library threads.

## Data publication

Generated catalogs, third-party HDF5 files, manuscripts, local archives and
detector-filtered results are excluded by `.gitignore`. Small historical parameter
and validation tables under `reference/` are retained as provenance; they are not
large event catalogs or new GWTC-4 inference results.

The manuscript's 500 combined CSVs occupy **14,230,310,053 bytes (13.25 GiB)**.
See [docs/DATA_INVENTORY.md](docs/DATA_INVENTORY.md) for per-file sizes, compression
estimates. The authors have selected Zenodo for the manuscript catalog; upload
and publication are pending. See [docs/ZENODO.md](docs/ZENODO.md) for account,
packaging, metadata and upload steps. No dataset DOI has yet been assigned.

The software and its accompanying documentation are released under the
[MIT License](LICENSE). This does not relicense third-party inputs, manuscripts,
or event datasets. The dataset license will be specified separately in the
Zenodo deposit. Scientific references are in [docs/REFERENCES.md](docs/REFERENCES.md).
