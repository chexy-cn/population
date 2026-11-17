# Verification of the publication cleanup

The cleanup was checked in the existing Conda GW numerical environment recorded
in `requirements-reproduction.txt`. This records completed checks, not forecasts.

- Synthetic-posterior tests cover variance/finite filtering and preservation of
  complete joint rows, distinct posterior selection, seeded replay, reference
  continuity, 20/30/50 integration convergence, empty Poisson realizations,
  overwrite refusal, and mutation detection in the appended suffix.
- The serial and two-worker synthetic pipelines produce identical catalog bytes.
- All 500 seeded posterior indices agree with the original 500-row parameter table.
- Optional local regression confirms the historical high-z catalog 001 byte for
  byte. Regenerated historical low-z columns are not all byte-identical; that
  boundary is documented in `REPRODUCIBILITY.md`.
- An end-to-end two-worker append using the actual original catalogs 001–003
  generated 306,919 high-z rows and 443,105 combined rows. SHA-256 comparisons
  of all three high-z and all three combined outputs agree with the existing
  manuscript files. The new validator checked all 15 fields in all those source
  rows and the exact original prefix and high-z suffix.
- The optional diagnostic plotting command generated all four figures and the
  combined-redshift figure was visually inspected.
- A Python wheel was built successfully without packaging input HDF5 files or
  event CSV catalogs. Dependency metadata and CLI entry points are included.
- The wheel was installed into an isolated temporary directory and its module
  CLI was exercised outside the source checkout.
- The GWTC-4-only CLI generated one actual 29,205-source catalog from the official
  input; all source fields passed validation.

No all-500-file regeneration was performed during cleanup. Existing event data
were not edited. The original historical validation results remain separately
available under `reference/` and should not be confused with the three-catalog
end-to-end cleanup regression.
