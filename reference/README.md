# Historical small reference products

These files are copies of the existing run's parameter tables and validation
reports. They describe the original 500-member ensemble used by the manuscript:

- `original_parameters.csv`: posterior indices, parameters, low counts, duration.
- `paper_extension_parameters.csv`: matching posterior parameters and appended counts.
- `historical_gwtc4_validation.json`: original low-catalog acceptance report.
- `historical_append_validation.json`: historical append acceptance report.

They contain no source event rows or posterior HDF5 data. They can be used to
audit the seeded hyperparameter selection and count totals. Historical reports
have their original scope; new validation is implemented in
`bbh_population/validation.py` and checks more fields and the append suffix too.

The full source catalogs remain local pending the authors' publication decision.
