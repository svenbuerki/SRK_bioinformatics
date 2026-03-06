# SRK Allele Data — Provenance & Pipeline

## Raw Data Files

| File | Description | Source |
|---|---|---|
| `SRK_individual_genotypes.tsv` | Binary presence/absence matrix (119 individuals x 172 alleles). Each cell is 0/1. | Sequencing / genotyping pipeline |
| `SRK_individual_allele_table.tsv` | Long-format read-level data (648 rows). Each row = one sequencing read mapping an allele to an individual. | Sequencing / genotyping pipeline |
| `SRK_allele_accumulation_curves.pdf` | Allele accumulation (rarefaction) curves | Upstream analysis |
| `SRK_chisq_species_population_frequency_plots.pdf` | Chi-squared species/population frequency plots | Upstream analysis |

## Generated Files

| File | Produced by | Description |
|---|---|---|
| `imputed_population.pkl` | `notebooks/00_load_data.ipynb` | Pickle containing imputed genotypes, allele pool, name/ID mappings, per-individual metadata, population groups, core allele list |
| `imputed_genotypes.tsv` | `notebooks/00_load_data.ipynb` | Human-readable TSV of imputed tetraploid genotypes (one row per individual, 4 allele columns + integer IDs + population group) |

To regenerate the derived files, run `notebooks/00_load_data.ipynb` from top to bottom.

## Pipeline Summary

1. **Load** both TSVs; validate same 119 individuals, binary matrix values, and read-allele consistency.
2. **Build read counts** per (individual, allele) pair from the long-format table.
3. **Impute tetraploid genotypes** (4 allele copies each) using Strategy A:

| Detected alleles | Rule |
|---|---|
| 1 | Homozygous: (A,A,A,A) |
| 2 | Proportional to reads, fill 4 slots (largest remainder method) |
| 3 | Proportional to reads, fill 4 slots (largest remainder method) |
| 4 | Direct: one copy each |
| >4 | Top 4 by read count, proportional fill (Strategy A) |

   **Largest remainder method**: multiply each allele's read proportion by 4, floor, distribute remaining slots to largest fractional remainders, ties broken alphabetically.

4. **Validate** all genotypes: exactly 4 elements, canonical sorted form, allele IDs in pool.
5. **Export** pickle and TSV to this directory.

## Key Data Quality Notes

- Only 6/119 individuals had exactly 4 detected alleles (high confidence).
- ~126/172 alleles are singletons (detected in only 1 individual).
- ~35% of allele calls rest on a single sequencing read.
- Individuals span 5 library groups (Library001–005) plus 1 hybrid (BEA_hybrid).

## Downstream Analysis

`notebooks/06_real_analysis.ipynb` loads `imputed_population.pkl` and runs the full crossing/equilibrium/optimization pipeline. See `CLAUDE.md` at the project root for complete project documentation.
