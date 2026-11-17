# Publishing the manuscript catalog on Zenodo

The authors selected Zenodo for the full historical manuscript intrinsic catalog.
Account creation, dataset-license selection, upload and publication are pending.
The software license is MIT; that choice is separate from the event-data license.

## Account

Create a production account at <https://zenodo.org/signup/>. Zenodo supports
ORCID, GitHub or a local account. Complete the required profile/email verification.
Using GitHub to sign in does not automatically transfer ignored catalog files
from this repository to Zenodo.

## Dataset to deposit

Deposit `catalog_gwtc4_extended_z30_appended`: 500 complete intrinsic catalogs,
49,811,117 sources, about 14.23 GB uncompressed. Include its
`extension_parameters_summary.csv` and `extension_metadata.txt`, a dataset README
explaining the schema/model/seeds, and checksums. The README should link to a
specific GitHub code version and the GWTC-4 input provenance.

Pack numbered CSVs into one or several archives before uploading: Zenodo's
documented default is 50 GB and 100 files per record, so 500 separate CSV uploads
would exceed the file-count quota. Several self-contained archives are convenient
for retries. The estimated compressed size is about 6.7 GB, not yet measured on
the full ensemble. Do not include the SNR-filtered results or duplicate standalone
high-z files unless they are explicitly part of the intended dataset release.

## Web upload

1. Sign in, choose **New upload**, and set resource type **Dataset**.
2. Upload the archives, description and checksum files.
3. Set an informative title, for example *GWTC-4-informed BBH catalogs with a
   high-redshift extension accompanying Constraining the dynamical Chern-Simons
   gravity with future gravitational-wave detectors*.
4. Enter the actual dataset creators, affiliations and ORCIDs, description and
   version. Dataset creators are an author decision, not inferred automatically
   from software commit authors or the manuscript author list.
5. Choose the dataset license. CC BY 4.0 is a possible choice for attributed
   reuse; it has not yet been selected by the authors. Retain the GWTC-4 input
   provenance and its applicable notices.
6. Save the draft and reserve a new dataset DOI if desired. The dataset DOI is
   distinct from the manuscript DOI and only registers when the record is published.
7. Review filenames, model description, counts, license and metadata, then publish.
8. Put the resulting dataset DOI and versioned download link in this README and
   the paper's data-availability statement.

## Server upload alternative

The catalog lives on the computing server. Uploading from that server through
Zenodo's REST API avoids downloading several GB to a laptop and uploading it again.
The account owner creates an API token locally and supplies it to the upload
process through a private environment or credential file. Never put the token
in source code, a Git commit, a command URL, a log or a chat message.

API upload can populate a draft; final publication can still be performed after
review in the browser. No account token or draft record exists in this project
yet, and this document does not itself create or publish a record.

Official instructions:

- Account: <https://help.zenodo.org/docs/get-started/create-an-account/>
- Upload and metadata: <https://help.zenodo.org/docs/deposit/create-new-upload/>
- REST API: <https://developers.zenodo.org/>
