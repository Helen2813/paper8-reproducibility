# Manual publication clearance checklist

- [ ] Confirm authorship and permission for included scripts, especially any
      re-used upstream code fragments in `replays/` or `diagnostics/`.
- [x] Select MIT for original software owned by Elena Krikun; the repository
      now includes an author-approved `LICENSE`. This does not clear upstream
      third-party code or data.
- [ ] Review and approve citation metadata and preprint authorship.
- [ ] Confirm repository contains no expression matrices, per-patient metadata,
      raw CEL archives, proprietary code, tokens, or author-local reports.
- [ ] Keep private numerical input files outside every Git repository.
- [ ] Document that the classifier's outcome-informed ComBat is leakage-prone,
      not a valid prospectively unbiased diagnostic estimate.
- [ ] Disclose incomplete historical null-model driver provenance and partial
      numerical mismatch in null replays.
- [ ] Run local/CI static tests and review the exact `git diff --cached` files.
- [ ] Only then decide whether to issue a GitHub Release / Zenodo software DOI.

The automated static preflight is **not** a privacy, source-rights, or security
clearance. MIT has been selected by the author; these checks do not clear third-party rights or approve a tagged public release.
