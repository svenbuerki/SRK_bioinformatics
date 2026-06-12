#!/usr/bin/env python3

import csv
from collections import defaultdict, OrderedDict

# ----------------------------
# User input
# ----------------------------
hap_key_file = input("Enter SRK_exon_haplotype_key.tsv: ").strip()
allele_key_file = input("Enter exon haplotype → allele key TSV: ").strip()

long_out = "SRK_individual_allele_table.tsv"
wide_out = "SRK_individual_genotypes.tsv"

# ----------------------------
# Load exon_haplotype → allele
# ----------------------------
hap_to_allele = {}

with open(allele_key_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        hap_to_allele[row["Sequence_ID"]] = row["Allele"]

print(f"Loaded {len(hap_to_allele)} haplotype → allele mappings")

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
            continue

        allele = hap_to_allele[exon_hap]

        # Extract individual: 2nd and 3rd underscore-separated fields
        parts = seq_id.split("_")
        if len(parts) < 3:
            print(f"WARNING: unexpected ID format, skipping: {seq_id}")
            continue

        individual = f"{parts[1]}_{parts[2]}"
        individual_alleles[individual].append(allele)

# ----------------------------
# Write long-format table
# ----------------------------
with open(long_out, "w") as out:
    out.write("Individual\tAllele\n")
    for individual in sorted(individual_alleles):
        for allele in individual_alleles[individual]:
            out.write(f"{individual}\t{allele}\n")

# ----------------------------
# Build genotype matrix
# ----------------------------
# Get sorted list of all alleles
all_alleles = sorted(
    {allele for alleles in individual_alleles.values() for allele in alleles}
)

# Convert each individual’s alleles to presence/absence
genotype_matrix = OrderedDict()

for individual in sorted(individual_alleles):
    allele_set = set(individual_alleles[individual])
    genotype_matrix[individual] = {
        allele: 1 if allele in allele_set else 0
        for allele in all_alleles
    }

# ----------------------------
# Write wide-format genotype table
# ----------------------------
with open(wide_out, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")

    header = ["Individual"] + all_alleles
    writer.writerow(header)

    for individual, geno in genotype_matrix.items():
        row = [individual] + [geno[allele] for allele in all_alleles]
        writer.writerow(row)

# ----------------------------
# Summary
# ----------------------------
print("\n✅ Allele scoring and genotyping complete")
print(f"Individuals scored: {len(genotype_matrix)}")
print(f"Unique alleles:     {len(all_alleles)}")
print(f"Long table:         {long_out}")
print(f"Genotype matrix:    {wide_out}")

