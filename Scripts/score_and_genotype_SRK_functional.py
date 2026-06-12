#!/usr/bin/env python3

import csv
from collections import defaultdict

# ----------------------------
# User input
# ----------------------------

hap_key_file = input("Enter SRK_exon_haplotype_key.tsv: ").strip()
func_allele_key_file = input("Enter SRK_protein_functional_alleles.tsv: ").strip()

out_long = "SRK_individual_functional_alleles.tsv"
out_matrix = "SRK_individual_functional_genotypes.tsv"

# ----------------------------
# Load DNA haplotype → functional allele
# ----------------------------

hap_to_func = {}

with open(func_allele_key_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        hap_to_func[row["DNA_Haplotype_ID"]] = row["Functional_allele"]

print(f"Loaded {len(hap_to_func)} haplotype → functional allele mappings")

# ----------------------------
# Parse individuals and score alleles
# ----------------------------

individual_alleles = defaultdict(set)

with open(hap_key_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        seq_id = row["Original_sequence_ID"]
        exon_hap = row["Exon_haplotype"]

        if exon_hap not in hap_to_func:
            continue  # haplotype filtered out or non-functional

        func_allele = hap_to_func[exon_hap]

        # Extract individual ID (2nd + 3rd underscore-delimited fields)
        parts = seq_id.split("_")
        if len(parts) < 3:
            print(f"WARNING: unexpected ID format, skipping: {seq_id}")
            continue

        individual = f"{parts[1]}_{parts[2]}"
        individual_alleles[individual].add(func_allele)

# ----------------------------
# Write long-format allele table
# ----------------------------

with open(out_long, "w") as out:
    out.write("Individual\tFunctional_allele\n")
    for ind in sorted(individual_alleles):
        for allele in sorted(individual_alleles[ind]):
            out.write(f"{ind}\t{allele}\n")

# ----------------------------
# Build genotype matrix
# ----------------------------

individuals = sorted(individual_alleles.keys())
alleles = sorted({a for alleles in individual_alleles.values() for a in alleles})

with open(out_matrix, "w") as out:
    out.write("Individual\t" + "\t".join(alleles) + "\n")

    for ind in individuals:
        row = []
        for allele in alleles:
            row.append("1" if allele in individual_alleles[ind] else "0")
        out.write(ind + "\t" + "\t".join(row) + "\n")

# ----------------------------
# Summary
# ----------------------------

print("\n✅ Functional allele genotyping complete")
print(f"Individuals genotyped: {len(individuals)}")
print(f"Functional alleles: {len(alleles)}")
print(f"Long format: {out_long}")
print(f"Genotype matrix: {out_matrix}")

