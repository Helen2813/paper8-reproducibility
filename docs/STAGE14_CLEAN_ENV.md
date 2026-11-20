# Stage 14: new R/Python environments, verified on author's Mac

**Reported status:** `PASS_FRESH_ENV_PORTABLE_DERIVED_CORE` on November 20,
2026. This status and numbers were supplied in the author's Terminal log;
the machine JSON receipt, exact micromamba solution, and installation logs
are not distributed in this review package. The local Stage 15 audit verifies
the retained report if present.

## Procedure and frozen versions

The Stage 14 v4 script created **two independent new micromamba prefixes**
under Downloads. Python 3.12.16 with exact pip versions
`numpy==2.5.4`, `pandas==3.1.1`, `scipy==1.18.1`, and
`scikit-learn==1.9.1`; R 4.4.3 with exact Bioconda packages
`sva==3.54.0`, `affy==1.84.0`, `Biobase==2.66.0`, and
`hgu133plus2cdf==2.18.0`. The installer used conda-forge, bioconda, and
PyPI, and checked import paths were inside the newly created prefixes.

After installation it re-ran `portable/portable_cli.py verify` and `all`
against nine SHA-256-verified private derived matrices and compared against
archived reference records.

**Results in supplied Terminal output:**
- ComBat max absolute difference: `1.77635683940025e-15` on
  9,577,656 stored entries.
- Frozen classifier max absolute score difference:
  `1.270095140171179e-13` for 37 test examples.
- Balanced accuracy: `0.9523809523809523`.

The micromamba log included the message `Error: ERROR: no packages specified`
while linking R/Bioconductor packages. The transaction still reported
success, the independent package version and library-path import check
passed, and both computational workers finished successfully. Retain the
unabridged installation logs as private provenance; this message should not
be silently erased from the record.

**Scope boundaries:** fresh dependency installation was verified on one
Mac (macOS Apple Silicon), not a second physical machine; fresh-environment
replay started from derived matrices, not GEO CEL files; null-model fitting
was not repeated. The historical Stage 14 v4 runner is a local transitional
script which expects a SHA-pinned previous Stage 13 staging ZIP and private
inputs; it is not presented here as a general public fresh-install command.

`portable/environment-combat.yml` and
`portable/requirements-classifier.txt` contain dependency specifications,
not a fully locked transitive dependency graph or signed install artifact.
Capture `micromamba list --json`, pip freeze, channel URLs and exact
package hashes before claiming byte-for-byte installation provenance.
