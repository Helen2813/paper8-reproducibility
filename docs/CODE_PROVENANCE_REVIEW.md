# Code provenance and release review (manual HOLD)

Elena Krikun approved MIT for original software she owns in this repository.
The copyright notice does **not** assert ownership of scientific methods,
third-party publications, downloaded scripts, published metadata, or other
materials that may appear in the underlying private project.

Before a tagged public release, review `portable/`, `replays/`, and
`diagnostics/` for literal third-party code copied or adapted from upstream
repositories, papers, notebooks, or supplementary archives. Record each
source, its copyright statement and license; either obtain permission, retain
required notices in an appropriate place, or replace/remove the affected
portion. Research methods and independently written implementations should be
distinguished from copied code. In particular, check all historical replay
scripts against the source files used during the original study; the integrity
tests alone cannot establish authorship or rights.

The software MIT license does not apply to data that are not included (raw
CEL files, per-GSM and donor-level tables, saved expression matrices, and
private numerical references). The absence of such files from this Git tree
does not by itself establish rights to redistribute them in another venue.

**Status:** MIT selected for original code; upstream code provenance and
manual publication clearance remain open. The static gate intentionally
reports a publication HOLD even when all integrity tests pass.
