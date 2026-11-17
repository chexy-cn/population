# External input and provenance

Official release: <https://doi.org/10.5281/zenodo.16911563>
(follow the GWTC-4 population release link at <https://gwosc.org/GWTC-4.0/>).
The release may package the HDF5 files in an archive; extract the BP2P file:

```text
BBHMassSpinRedshift_BrokenPowerLawTwoPeaks_GaussianComponentSpins_PowerLawRedshift.h5
```

Recorded local input:

- Size: 342,099,032 bytes.
- SHA-256: `80ba79716f5ba6e9fd8c8765cf1cc102dbf63f33afd5dfd49ee5f6797a1253ec`.
- Dataset: `posterior/hyperparameter_samples`.
- Column labels: `hyperparameters` attribute.
- Model labels: `model_names` attribute.
- Required posterior filter: finite mapped parameters; finite `variance < 1`.
- Valid rows: 3,962 for this file. Indices refer to the original HDF5 row numbers.

Verify with `sha256sum /path/to/the/file.h5`. A changed official release may use
a different filename/checksum; exact realization replay needs the recorded input.
Do not silently substitute a dominant-mode-only posterior, maximum-posterior
parameters, GWTC-3, or a different mass/spin model.

GWTC parameter names map to clearer output labels in
`bbh_population/sampling.py:HDF_PARAMETER_MAP`. For example `lamb -> kappa`,
`rate -> rate_0`, `lam_0/lam_1 -> lambda_0/lambda_1`, and
`mu_spin/sigma_spin/xi_spin -> mu_tilt/sigma_tilt/zeta_tilt`.

Third-party release contents remain local and are excluded from Git. Their
licenses and citations remain the responsibility of the data release; this
repository does not relicense those products.
