#!/usr/bin/env python3

import csv
from collections import defaultdict

# ----------------------------
# User input
# ----------------------------
hap_key_file = input("Enter SRK_exon_haplotype_key.tsv: ").strip()
allele_key_file = input("Enter exon haplotype → allele key TSV: ").strip()

out_file = "SRK_individual_allele_table.tsv"

# ----------------------------
# Load exon_haplotype → allele
# ----------------------------
hap_to_allele = {}

with open(allele_key_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        hap_to_allele[row["Sequence_ID"]] = row["Allele"]

# ----------------------------
# Parse individuals and score alleles
# ----------------------------
individual_alleles = defaultdict(list)

with open(hap_key_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        seq_id = row["Original_sequence_ID"]
        exon_hap = row["Exon_haplotype"]

        if exon_hap not in hap_to_allele:
            # exon haplotype not assigned to an allele (e.g. filtered out)
            continue

        allele = hap_to_allele[exon_hap]

        # Extract individual: 2nd and 3rd elements
        parts = seq_id.split("_")
        if len(parts) < 3:
            print(f"WARNING: unexpected ID format, skipping: {seq_id}")
            continue

        individual = f"{parts[1]}_{parts[2]}"
        individual_alleles[individual].append(allele)

# ----------------------------
# Write output
# ----------------------------
with open(out_file, "w") as out:
    out.write("Individual\tAllele\n")
    for individual in sorted(individual_alleles):
        for allele in individual_alleles[individual]:
            out.write(f"{individual}\t{allele}\n")

# ----------------------------
# Summary
# ----------------------------
print("\n✅ Allele scoring complete")
print(f"Individuals scored: {len(individual_alleles)}")
print(f"Output file: {out_file}")

