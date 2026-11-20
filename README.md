# Paper 8 — reproducibility code for an endometriosis classifier audit

Companion code for:

> Krikun E, Alkhateeb A. *Hidden Dependencies in Endometriosis Microarray
> Classifiers: Sample Reuse, Outcome-Informed Batch Correction, and Diagnostic
> Reassessment.* medRxiv preprint (2026).
> https://doi.org/10.64898/2026.11.15.26368137

**Scope:** reproducibility and methodological audit — **not** a clinically
validated diagnostic test. The preprint has not been peer reviewed.

## What is in this repository?

- `portable/` — executable, SHA-256-gated saved-matrix ComBat and frozen
  classifier reconstruction, with pinned scientific dependencies.
- `replays/` — author-local standalone replay scripts for raw-CEL RMA/MAS5,
  Huang-inspired own-copy omission, date controls, and computational audits.
  **Historical evidence scripts:** these are not a portable one-command
  raw-GEO-to-every-result workflow.
- `diagnostics/` — null-fit replication investigations, **not** the unrecovered
  historical v3.70/v3.71 null-model drivers.
- `docs/` — reported results, dependencies, reproducibility limits and
  how to supply authorized input files.
- `tests/`, `scripts/` — offline source-integrity and synthetic safety checks.

**This repository deliberately contains no CEL files, expression matrices,
clinical tables, per-GSM supplementary references, private logs, or archives.**
No source-specific licensed datasets have been redistributed.

## Static checks (no scientific dependencies or data required)

```bash
python3 -B scripts/public_code_preflight.py
python3 -B -m unittest discover -s tests -v
python3 -B portable/portable_cli.py verify
```

A PASS confirms only source integrity and test gates, **not** a fresh numerical
reproduction or independent privacy/licensing clearance.

## Numerical replay from authorized frozen inputs

The numerical core checks 9 saved-matrix resources and **2 additional**
reference tables by exact SHA-256. All 11 files must be supplied outside this
Git repository. Their approved provenance and hashes are recorded in
`portable/input_manifest.json`; none are included here.

For an author with an authorized full source project:

```bash
python3 -B portable/portable_cli.py export-inputs \
  --source-root /path/to/authorized/endo-research \
  --inputs-dir /path/outside/git/private_inputs
```

After separately copying the exact two verified Supplement S2/S3 reference
files as `reference_S2.csv` and `reference_S3.csv` **into that private input
directory**, run with the pinned R and Python environments:

```bash
python3 -B portable/portable_cli.py all \
  --inputs-dir /path/outside/git/private_inputs \
  --out-dir /path/outside/git/private_results \
  --rscript /path/to/new/R/bin/Rscript \
  --python /path/to/new/python/bin/python3
```

Use `combat`, `classifier`, or `all` as the subcommand. See
`docs/REPLAY_GUIDE.md` for installation and expected numeric outcomes. No
source project or downloaded data are modified.

## Evidence and honest limitations

Reported macOS replay outcomes (Nov 2026):

| Audit | Reported result | Important boundary |
|---|---|---|
| CEL → per-series RMA | 232 arrays, 12,684,600 entries, max difference 0 | Four series normalized separately |
| GSE51981 single-CEL MAS5 | 148 fresh R runs, 324 probes, 47,952 entries, max difference 0 | Specific frozen probe panel |
| Saved-matrix ComBat | 44,341 × 216 entries, max difference 1.78 × 10⁻¹⁵ | Test diagnoses intentionally included |
| Frozen elastic-net replay | 37 scores, max difference 1.27 × 10⁻¹³, BA 0.952381 | Recreates outcome-informed historical pipeline; not unbiased validation |
| New scientific environments | Passed on the author's Apple Silicon Mac | No second-machine verification |
| Null-model replay | Partial/HOLD | v3.70/v3.71 drivers not recovered; not 3,582 exact reconstructions |

These are the **author-reported Terminal outcomes**, not computations
performed by GitHub Actions. See `docs/REPLAY_EVIDENCE.md` and
`docs/RELEASE_LIMITATIONS.md` for further qualifications.

## Licensing and citation

`CITATION.cff` cites the associated scientific preprint; no software DOI has
been issued. The author has selected the **MIT License** for original software
and documentation she owns (see `LICENSE`). The license does **not** grant
rights to upstream third-party code or source-derived data, and it does not
authorize redistribution of private expression matrices, clinical metadata,
or restricted datasets. Any embedded upstream code must be identified and
handled under its own applicable terms before publication; see
`docs/CODE_PROVENANCE_REVIEW.md`.

**Release state:** MIT-selected code-only review candidate; third-party code
provenance and other manual publication gates remain pending. A tagged GitHub
release and Zenodo deposit should follow that review.
