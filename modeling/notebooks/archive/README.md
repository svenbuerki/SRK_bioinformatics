# Archived Notebooks (Old Data Path)

These notebooks were archived on 2026-05-23 when the project migrated from the
**protein-level binary representation** (94 SRK proteins, 0/1 detection matrix)
to the **biological-allele dosage representation** (49 observed alleles with
explicit copy counts, collapsed via synonymy groups into ~27 biological alleles).

## What's here

| Notebook | Old data path |
|---|---|
| `00_load_data.ipynb` | Reads `data/revised/` (94-protein binary matrix) and writes `data/population.pkl` |
| `06_demography.ipynb` | Loads `data/population.pkl` for logistic growth / Ne work |
| `07_collapse_prediction.ipynb` | Loads `data/population.pkl` for collapse cascade |
| `08_real_analysis.ipynb` | Full pipeline on the old `data/population.pkl` |
| `08.1_real_analysis_pop27.ipynb` | Pop 27 subset using `data/filespopulation27/` |

## Reviving an archived notebook

These files have not been rewritten for their new location. If you re-run them
from `notebooks/archive/`, the `sys.path.insert(0, "../src")` and
`pd.read_csv("../data/...")` paths inside them will point one level too high.
Either:

1. Copy the notebook back to `notebooks/` and re-run; **or**
2. Edit the path strings inside to `../../src` and `../../data/...`.

The underlying data they read (`data/revised/`, `data/filespopulation27/`) is
still present and untouched.

## Why they were archived, not deleted

The protein-level model and the GFS / TP2 framing of crossing strategy were
developed against this representation. The notebooks are the reference
record of those analyses and the figures they produced. The biological
conclusions still stand; only the input representation changed.

## Replacements on the new data path

| Old | New |
|---|---|
| `00_load_data.ipynb` | `notebooks/00_load_data_Salleles.ipynb` |
| `08_real_analysis.ipynb` | `notebooks/08.2_real_analysis_Salleles.ipynb` |
| `06_demography.ipynb`, `07_collapse_prediction.ipynb`, `08.1_real_analysis_pop27.ipynb` | To be re-implemented on the new path (see `Genotype_Quality_Crossing_Plan.md`) |
