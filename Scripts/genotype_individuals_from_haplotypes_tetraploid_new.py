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
subgenomes_per_ind = defaultdict(lambda: defaultdict(set))  # ind -> tig -> set of Protein_haplotypes

with open(tsv_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        seq_id = row["Original_sequence_ID"]
        prot_hap = row["Protein_haplotype"]

        parts = seq_id.split("_")
        if len(parts) < 4:
            print(f"WARNING: unexpected ID format, skipping: {seq_id}")
            continue

        element1 = parts[0]  # hap1, hap2, ...
        library = parts[1]
        barcode = parts[2]
        tig = parts[3]

        individual = f"{library}_{barcode}"

        haplotypes_per_ind[individual].add(prot_hap)
        subgenomes_per_ind[individual][tig].add(prot_hap)

# ----------------------------
# Compute Tetraploid genotype
# ----------------------------
tetraploid_genotype = {}
for ind, tig_dict in subgenomes_per_ind.items():
    # Assign letters to tigs sorted alphabetically
    letters = {tig: chr(65 + i) for i, tig in enumerate(sorted(tig_dict.keys()))}
    
    genotype_parts = []
    for tig in sorted(tig_dict.keys()):
        unique_prots = sorted(tig_dict[tig])
        if len(unique_prots) == 1:
            # Only one protein haplotype in this tig
            genotype_parts.append(letters[tig])
        else:
            # Multiple protein haplotypes: add numbers
            for i, _ in enumerate(unique_prots, start=1):
                genotype_parts.append(f"{letters[tig]}{i}")
    
    tetraploid_genotype[ind] = "".join(genotype_parts)

# ----------------------------
# Write output
# ----------------------------
with open(out_file, "w") as out:
    out.write("Individual\tHaplotypes\tSubgenomes\tTetraploid_genotype\n")
    for ind in sorted(haplotypes_per_ind):
        haps = sorted(haplotypes_per_ind[ind])
        subgenome_entries = []
        for tig, prots in sorted(subgenomes_per_ind[ind].items()):
            for p in sorted(prots):
                subgenome_entries.append(f"{p}:{tig}")  # keep protein hap + tig for reference
        out.write(f"{ind}\t{','.join(haps)}\t{','.join(subgenome_entries)}\t{tetraploid_genotype[ind]}\n")

print("\n✅ Tetraploid genotyping complete")
print(f"Individuals genotyped: {len(haplotypes_per_ind)}")
print(f"Output written to: {out_file}")

