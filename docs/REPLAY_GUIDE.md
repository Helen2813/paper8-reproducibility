# Running the saved-matrix replay in separate scientific environments

**Scope:** a locked numerical replay from archived, outcome-informed matrices;
not a fully independent clinical validation and not the original null driver.

## Dependencies

Direct versions successfully installed by the author on macOS Apple Silicon:

- Python 3.12.16, NumPy 2.5.4, pandas 3.1.1, SciPy 1.18.1, scikit-learn 1.9.1;
- R 4.4.3, `sva` 3.54.0 (other Bioconductor packages pinned in
  `portable/environment-combat.yml`).

The packages above are direct pins, not a signed full transitive lockfile;
availability can change. Stage 14 previously passed on the author's Mac,
not a second physical computer.

From a trusted `micromamba` installation (outside the Git worktree):

```bash
PY_PREFIX="$HOME/Downloads/P8_REPRO_PYTHON"
R_PREFIX="$HOME/Downloads/P8_REPRO_R"
micromamba create -y -p "$PY_PREFIX" --override-channels -c conda-forge python=3.12.16 pip
"$PY_PREFIX/bin/python3" -I -m pip install --only-binary=:all: --no-cache-dir \
  numpy==2.5.4 pandas==3.1.1 scipy==1.18.1 scikit-learn==1.9.1
micromamba create -y -p "$R_PREFIX" --override-channels -c conda-forge -c bioconda \
  r-base=4.4.3 bioconductor-sva=3.54.0 bioconductor-affy=1.84.0 \
  bioconductor-biobase=2.66.0 bioconductor-hgu133plus2cdf=2.18.0
```

The exact above block is a documented recipe; it was not validated verbatim
on a clean second machine. The Stage 14 installer used a more detailed
isolation and package check.

## Prepare authorized scientific and reference inputs (outside Git)

`portable/input_manifest.json` lists 9 raw-derived matrices and metadata with
original project-relative provenance and SHA-256, plus two additional S2/S3
reference files. The latter were included in an earlier private staging ZIP;
**they are excluded from this code-only release** pending redistribution review.

```bash
python3 -B portable/portable_cli.py export-inputs \
  --source-root "/path/to/authorized/endo-research" \
  --inputs-dir "/path/to/new/private_inputs"
```

Obtain the two authorized reference files and copy them into that same
**external** input directory as `reference_S2.csv` and `reference_S3.csv`.
The CLI verifies both reference hashes before classifier replay. File source
and format must be legitimate; do not disable hash gates to make a replay pass.

## Run

```bash
"$PY_PREFIX/bin/python3" -I portable/portable_cli.py all \
  --inputs-dir "/path/to/new/private_inputs" \
  --out-dir "/path/to/new/private_results" \
  --rscript "$R_PREFIX/bin/Rscript" \
  --python "$PY_PREFIX/bin/python3"
```

Expected author-reported numeric replay on the frozen inputs:
ComBat maximum absolute error `1.77635683940025e-15` (44,341 × 216 entries),
classifier score maximum absolute error `1.270095140171179e-13`
(37 test cases), BA `0.9523809523809523`.

**No original training/test data are included in the GitHub code package.**
A failure to obtain verified inputs is a documented reproducibility limitation,
not a reason to substitute synthetic test data for reported science.
