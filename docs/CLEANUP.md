# Publication cleanup record

Prepared on 2026-11-17 from the existing workspace and the supplied manuscript.
No Git repository or remote was initialized and nothing was uploaded.

Public working code now lives in `bbh_population/`, with `configs/`, `tests/`,
`docs/`, `reference/` and one portable SLURM example. The old scripts are local
history, not alternate public defaults. The source package includes dependencies,
a reproducibility stack, a command-line interface and GitHub Actions unit tests.

Moved, rather than deleted:

- Ten former root Python scripts (generation, filtering, plotting, validation)
  to `_local_archive/20261117/scripts/`.
- Ten machine-specific SLURM submission scripts to
  `_local_archive/20261117/slurm/`.
- All former root `slurm-*.out` and `slurm-*.err` files to
  `_local_archive/20261117/logs/`.
- The exploratory `BP2P.ipynb` to `_local_archive/20261117/notebooks/`.
- The superseded whole-range redshift diagnostics and root bytecode cache to
  `_local_archive/20261117/results/` and `_local_archive/20261117/cache/`.

Recovery is by moving the named file back out of the archive. Old job scripts
refer to their historical locations and are reference material; use
`examples/generate.slurm` for new runs. The archive is ignored by Git.

Preserved in place and excluded from source distribution: all event catalog
directories, SNR/detector products, `o4a-astro/`, `data/`, `paper/`, legacy
`filter_catalog_for_dCS/` material and generated plot/validation directories.
There is no deletion of source data or catalog output. Small copies of the
historical parameter tables and validation summaries are in `reference/`.

Numerical scope: the public sampler retains the existing mass and spin formulas
and current sampling algorithm. The default GWTC-4 mode again uses the
posterior-specific power-law exponent. The manuscript mode retains the
conditional count prescription rather than replacing it with a new absolute
high-z count. The whole-range Madau-Dickinson and fresh-100 experiments are
archived to prevent accidental use as the manuscript preset.

Open author decisions: event-data deposition and code license/software citation
metadata. The supplied draft's public identifier can be added after publication.
