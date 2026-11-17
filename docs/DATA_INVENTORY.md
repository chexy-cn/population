# Local event data and publication options

Inventory date: 2026-11-17. Sizes below are actual file byte counts (not rounded
disk allocation). All existing catalog directories have been kept at their
original paths. They are excluded from Git by default.

| Dataset | Catalog files | Bytes | GiB (2^30 bytes) | Role |
| --- | ---: | ---: | ---: | --- |
| `bbh_catalogs_gwtc4_bp2p` | 500 | 4,378,389,349 | 4.08 | Original low-z sources and parameter table |
| `catalog_gwtc4_highz_extension_z1p5_30` | 500 | 9,851,977,518 | 9.18 | Added high-z rows |
| `catalog_gwtc4_extended_z30_appended` | 500 | 14,230,310,053 | 13.25 | **Manuscript intrinsic catalog** |
| `catalog_gwtc4_extended_z30_appended_snr_raw_gr_filtered` | 3,000 | 79,064,303,442 | 73.63 | Downstream detector-specific analysis products |

The manuscript intrinsic catalog has 49,811,117 source rows. Mean numbered CSV
size is 28,460,175 bytes (27.14 MiB); the largest is 97,956,945 bytes (93.42 MiB).
Directory totals include small parameter/metadata files. Per-file catalog counts
in detector directories count repeated detector-specific products, not distinct
intrinsic sources.

Streaming gzip level-6 measurements, without saving compressed copies:

| Combined CSV | Original bytes | Compressed bytes |
| --- | ---: | ---: |
| catalog_001.csv | 27,163,946 | 12,826,001 |
| catalog_376.csv (median-size file) | 26,755,277 | 12,632,578 |
| catalog_495.csv (largest file) | 97,956,945 | 46,221,529 |

These representative measurements suggest about 6.7 GB (6.25 GiB) for all
combined CSVs compressed separately. That is an estimate, not a completed archive
measurement. No catalog archive has been uploaded or selected for publication.

Possible author choices:

1. Publish only code, configuration, small parameter/validation tables and input
   provenance. Readers can regenerate the model; exact historical low-z rows
   remain unavailable unless separately supplied.
2. Deposit all combined CSVs in a research-data repository (e.g. an institutional
   archive or Zenodo), cite its versioned DOI from the code repository, and leave
   the Git repository lightweight. Combined files already contain both components;
   publishing separate high-z copies too would duplicate most bytes.
3. Publish a clearly labeled small example plus the full model-generation code.
   This does not replace access to all historical source rows.

The JSON inventory for all 21 local catalog/result directories is
`docs/catalog_inventory.json`. Rerun:

```bash
python scripts/catalog_inventory.py --output docs/catalog_inventory.json
```

The official posterior HDF5 is a further 342,099,032 bytes and is an external
input with its own provenance and distribution terms. The entire downloaded
release, other papers, and detector/SNR products are not part of the source package.
