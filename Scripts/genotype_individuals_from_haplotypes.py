#!/usr/bin/env python3

import csv
from collections import defaultdict

# ----------------------------
# User input
# ----------------------------
tsv_file = input("Enter haplotype key TSV file: ").strip()
out_file = tsv_file.replace(".tsv", "_genotypes.tsv")

# ----------------------------
# Parse haplotypes per individual
# ----------------------------
genotypes = defaultdict(set)

with open(tsv_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        seq_id = row["Original_sequence_ID"]
        hap = row["Protein_haplotype"]

        parts = seq_id.split("_")

        if len(parts) < 4:
            print(f"WARNING: unexpected ID format, skipping: {seq_id}")
            continue

        library = parts[1]
        barcode = parts[2]

        individual = f"{library}_{barcode}"

        genotypes[individual].add(hap)

# ----------------------------
# Write genotype table
# ----------------------------
with open(out_file, "w") as out:
    out.write("Individual\tHaplotypes\tN_haplotypes\n")

    for ind in sorted(genotypes):
        haps = sorted(genotypes[ind])
        out.write(
            f"{ind}\t{','.join(haps)}\t{len(haps)}\n"
        )

# ----------------------------
# Summary
# ----------------------------
print("\n✅ Genotyping complete")
print(f"Individuals genotyped: {len(genotypes)}")
print(f"Genotype table written to: {out_file}")

