# Reproducing models and historical realizations

The manuscript supplied on 2026-11-17 identifies the append ensemble in its
Table I and Appendix A. There are two different reproducibility targets:

1. **Model reproduction:** `configs/gwtc4.json` and `configs/paper_append.json`
   generate the specified probability distributions and count prescriptions.
2. **Historical rows:** supply the actual `bbh_catalogs_gwtc4_bp2p` originals
   with `--original-dir` to preserve the low-z realization used by the manuscript.
   Same posterior input, seeds, ensemble size and numerical stack reproduce the
   high-z sequence. This exact path is available locally; public access to the
   original low-z files is still an author decision.

The surviving original generator evolved during the initial 500-catalog work.
Regression confirms the same 500 posterior indices and `catalog_001` primary
masses, but freshly generated low catalog source columns do not all match the
saved original rows. The conditional-q values and subsequent RNG sequence differ;
the complete earlier implementation/RNG history is not frozen. We do
not claim an all-row byte replay of the old low-z files from the current code.
The public sampler preserves the current documented conditional-q inverse-CDF
algorithm and all mass/spin model assumptions; historical source files are
retained as the authority for exact manuscript rows.

The tests check the high-z `catalog_001.csv` against the historical file byte for
byte. The original input files are never rewritten. Full public regeneration
should be assessed at the distribution/model level unless historical originals
are supplied. A byte-level regression on one high catalog is explicitly not
proof of a full 500-file byte-level replay.

## Stable random streams

Root seed 20260918 selects 500 joint rows without replacement, then 500 uint64
source seeds. Root seed 20260920 supplies 500 extension seeds. Each catalog's
source stream is independent of scheduling and worker count. Changing the
ensemble size changes this mapping; use `--catalog-ids 1-3` with a 500-member
preset to reproduce that subset of the historical ensemble. Changing either
seed creates a new ensemble/extension and should be reported in any analysis.

`run.json` saves effective settings, selected IDs, Python and dependency versions,
posterior hash, cosmology, status, and any original directory. Parameter summaries
save the posterior index and actual child RNG seeds. For an append run, hashes
of original, high and combined files are also saved. These additions do not
change the historical 15-column source schema.

## Historical acceptance data

The `reference/` files are small copies of the existing parameter tables and
validation reports, not input event catalogs. The expected manuscript totals are:

- 500 realizations, each representing Tcat=1 observer year.
- Original rows: 15,317,438 (mean 30,634.876 per realization).
- High-z rows: 34,493,679 (mean 68,987.358).
- Combined: 49,811,117 (mean 99,622.234).
- Reference I_high/I_low=2.2516710958959476.
- 30 -> 50 boundary change in that ratio: about 3.97e-7.

These counts are expectations for the recorded historical files, not assertions
about every new stochastic realization. The new validator scans all rows,
checks physical ranges and mass relationships, and compares the complete append
prefix and suffix. Historical validation had a narrower scope, including byte
prefixes and redshift ranges; its reports are retained as historical evidence.
