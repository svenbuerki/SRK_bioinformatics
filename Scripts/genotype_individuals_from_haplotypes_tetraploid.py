#!/usr/bin/env python3

import csv
from collections import defaultdict

# ----------------------------
# User input
# ----------------------------
tsv_file = input("Enter haplotype key TSV file: ").strip()
out_file = tsv_file.replace(".tsv", "_tetraploid_genotypes.tsv")

# ----------------------------
# Parse haplotypes and subgenomes per individual
# ----------------------------
haplotypes_per_ind = defaultdict(set)
subgenomes_per_ind = defaultdict(list)

with open(tsv_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        seq_id = row["Original_sequence_ID"]
        hap_label = row["Protein_haplotype"]  # keep for Haplotypes column

        parts = seq_id.split("_")
        if len(parts) < 4:
            print(f"WARNING: unexpected ID format, skipping: {seq_id}")
            continue

        element1 = parts[0]  # first element, e.g., hap1
        library = parts[1]
        barcode = parts[2]
        tig = parts[3]       # fourth element, e.g., tig00000002

        individual = f"{library}_{barcode}"

        # Store haplotypes
        haplotypes_per_ind[individual].add(hap_label)

        # Store subgenome info using Original_sequence_ID elements
        subgenomes_per_ind[individual].append((element1, tig))

# ----------------------------
# Compute Tetraploid genotype (letters + numbers)
# ----------------------------
tetraploid_genotype = {}
for ind, element_tig_list in subgenomes_per_ind.items():
    # Group haplotypes by tig (subgenome)
    tig_dict = defaultdict(list)
    for element1, tig in element_tig_list:
        tig_dict[tig].append(element1)

    # Assign letters to tigs sorted alphabetically
    letters = {}
    for i, tig in enumerate(sorted(tig_dict.keys())):
        letters[tig] = chr(65 + i)  # 65 = 'A'

    genotype_parts = []
    for tig in sorted(tig_dict.keys()):
        hap_list = sorted(tig_dict[tig])
        if len(hap_list) == 1:
            genotype_parts.append(letters[tig])
        else:
            for j, hap in enumerate(hap_list, start=1):
                genotype_parts.append(f"{letters[tig]}{j}")

    tetraploid_genotype[ind] = "".join(genotype_parts)

# ----------------------------
# Write output table
# ----------------------------
with open(out_file, "w") as out:
    out.write("Individual\tHaplotypes\tSubgenomes\tTetraploid_genotype\n")
    for ind in sorted(haplotypes_per_ind):
        haps = sorted(haplotypes_per_ind[ind])
        # Build Subgenomes column as element1:tig from Original_sequence_ID
        subgenome_entries = sorted(f"{element1}:{tig}" for element1, tig in subgenomes_per_ind[ind])
        out.write(
            f"{ind}\t{','.join(haps)}\t{','.join(subgenome_entries)}\t{tetraploid_genotype[ind]}\n"
        )

# ----------------------------
# Summary
# ----------------------------
print("\n✅ Tetraploid genotyping complete")
print(f"Individuals genotyped: {len(haplotypes_per_ind)}")
print(f"Output written to: {out_file}")

