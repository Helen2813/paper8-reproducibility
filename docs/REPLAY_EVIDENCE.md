# Paper 8 — computational replay evidence (reported results)

| Check | Scope | Reported result | Limits |
|---|---|---|---|
| Derived-table check | S1–S6/S11, 20 hashes | PASS on macOS | Arithmetic check, not model recomputation |
| Stage 6, ComBat | 44,341 rows x 216 arrays | 9,577,656 values, max absolute error 1.77635683940025e-15 | Archived RMA inputs; tests jointly used during outcome-informed ComBat |
| Stage 7, frozen classifier | 179 train, 37 test, 270 genes | 37/37 predictions match; maximum continuous score error 1.27009514017e-13; BA 0.952381 | Historical frozen author configuration has `max_iter=100`; convergence not established |
| Stage 8, null model diagnostics | 18 RMA-only pilot fits, 12 ComBat pilot fits | RMA-only BA 18/18 matched; ComBat BA 11/12 matched | All 12 ComBat continuous score sets differ; 1 changed classification; historical drivers v3.70/v3.71 unavailable; provenance of assignment-file hashes unresolved |
| Stage 8H, SAGA seed sensitivity | 12 configurations x 6 seeds | Seed spread comparable with some score differences | This does not prove historical random state or drivers |
| Stage 9A, Huang own-copy-out | 26 own-copy test observations, 37 omission fits | Pooled AUC 0.323308 vs 1.000 with copies; matched 78 scores and 597 placebo values | Independent reconstruction; not original Huang fitted model |
| Stage 9B, Huang seed and bootstrap | 5 RF seeds; 5,000 independent bootstrap draws | Five seeded AUC values match; independent bootstrap 95% CI 0.1053–0.5641 vs historical 0.1053–0.5714 | Conditional intervals, historical bootstrap RNG sequence not confirmed |
| Stage 10A, nested 324 probes | 148 samples, 25 folds, 527 OOF rows | Sample-level scores match exactly | Internal, not patient-independent; some date-blocked splits ineligible |
| Stage 10B, dates | Five 2009 dates | Pooled AUC 0.697756 and within-date summary reproduced | Original v4.23 OOF reference ZIP unavailable; individual score identity not independently checked |
| Stage 10C, date-only negative controls | 4 tasks | 3,326 numeric values matched exactly | Data from archived MAS5/precomputed sample metadata, not raw CEL |
| Stage 11A–C, raw CEL -> RMA | Four GEO series; 232 CEL | 12,684,600 values matched archived RMA, max absolute error 0 | Original per-series normalization is not equivalent to independent train/test processing within a batch |
| Stage 11E v2, raw CEL -> MAS5 | GSE51981, 148 CEL, 324 probes | 148 fresh R invocations, 0 reused checkpoints; 47,952 values matched exactly | This covers MAS5 GSE51981 324 probes, not every MAS5 source used in Huang |

These are **reported run outcomes**; exact external run artifacts have not all been incorporated into this staging archive. Verification on the author's machine is not a cold-environment reproduction on a second system.

## Stage 14 (November 20, 2026): clean dependency environment

The author supplied a Terminal log of `PASS_FRESH_ENV_PORTABLE_DERIVED_CORE`:
new independent micromamba R and Python environments passed exact-version
and library-path validation; full archived derived-matrix ComBat and 37
frozen classifier predictions matched the prior numerical references with
max absolute differences `1.77635683940025e-15` and
`1.270095140171179e-13`, respectively. This is a reported successful
run on the author's Mac. The underlying private JSON receipt and original
installation logs are **not** bundled. The release remains unable to claim
second-machine portability, fully automated raw-GEO replay or original null
reconstruction. The Stage 15 local audit can verify the saved JSON receipt.
