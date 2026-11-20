# Portable numerical core — saved matrices, not raw GEO

The `portable_cli.py` `all` mode replays outcome-informed ComBat and the frozen
270-gene elastic-net classifier, starting from preserved derived-matrix files.
It is **not** a prospective endometriosis diagnostic workflow.

This code-only public snapshot contains **no expression matrices and no
per-GSM reference tables**. The 9 scientific matrices/metadata/coefficient
files and 2 S2/S3 published comparison tables must be supplied outside the
repo, with exact SHA-256 values from `input_manifest.json`.

`python3 -B portable/portable_cli.py verify` runs offline static checks.

`python3 -B portable/portable_cli.py export-inputs --source-root SRC --inputs-dir DEST`
creates private verified copies of the nine scientific inputs under a new,
previously empty directory outside Git. It does not export the additional
reference tables; these must be supplied separately and checked by `all` or
`classifier`. The program will HOLD if either is missing or altered.

See `docs/REPLAY_GUIDE.md` for exact dependency versions and commands.
