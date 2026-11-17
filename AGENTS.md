# Repository guidelines

This repository generates intrinsic GWTC-4 BP2P BBH populations and the
manuscript's high-redshift append model. Reusable Python code is in
`bbh_population/`; JSON scientific presets are in `configs/`; reproducibility,
input provenance and publication notes are in `docs/`.

Use Python >=3.9. Install with `python -m pip install -e .`. The historical
numerical stack is recorded in `requirements-reproduction.txt`. Run tests with
`python -m unittest discover -s tests -v`; optional local regression tests need
the historical posterior and source catalogs. Source changes must be verified
against the supported scientific prescription, not just self-consistency.

Preserve complete correlated posterior draws. Never silently substitute marginal
medians or a different posterior. Keep source-frame masses distinct from derived
detector-frame masses. Distinguish the posterior-specific low-z kappa_j from the
fixed reference gamma used only to anchor the manuscript extension.

Existing catalog and third-party input directories are local-only, ignored by
Git, and must not be overwritten. Output directories are exclusively created.
Keep archived historical scripts in `_local_archive/`; do not reintroduce their
machine-specific paths/defaults into the public CLI. Generated catalogs and
SLURM logs do not belong in commits. Add provenance for any new static data.

Use four-space indentation, snake_case and concise public-helper docstrings.
State physical units. Do not change scientific assumptions during code cleanup.
Document the difference between model regeneration and exact replay of saved
historical realizations (see `docs/REPRODUCIBILITY.md`). The population code
does not implement downstream waveform/SNR/Fisher analyses.
