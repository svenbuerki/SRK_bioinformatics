#!/usr/bin/env python3

from Bio import Phylo

# ----------------------------
# User input
# ----------------------------

tree_file = input("Enter IQ-TREE .contree file: ").strip()
bootstrap_cutoff = float(input("Enter bootstrap cutoff (e.g. 85): ").strip())

out_file = tree_file.replace(".contree", "_allele_key.tsv")

# ----------------------------
# Load tree
# ----------------------------

tree = Phylo.read(tree_file, "newick")

alleles = {}
assigned = set()
allele_counter = 0

# ----------------------------
# Identify TERMINAL supported clades only
# ----------------------------

for clade in tree.find_clades(order="postorder"):

    if clade.confidence is None:
        continue

    if clade.confidence < bootstrap_cutoff:
        continue

    tips = [t.name for t in clade.get_terminals()]

    # KEY FIX:
    # Skip if ANY child clade also meets bootstrap cutoff

    child_supported = False

    for child in clade.clades:
        if child.confidence is not None and child.confidence >= bootstrap_cutoff:
            child_supported = True
            break

    if child_supported:
        continue

    # Assign allele

    allele_counter += 1
    allele_name = f"Allele_{allele_counter:03d}"

    for tip in tips:
        alleles[tip] = allele_name
        assigned.add(tip)

# ----------------------------
# Assign singletons
# ----------------------------

for tip in tree.get_terminals():

    if tip.name not in alleles:

        allele_counter += 1
        alleles[tip.name] = f"Allele_{allele_counter:03d}"

# ----------------------------
# Write output
# ----------------------------

with open(out_file, "w") as out:

    out.write("Sequence_ID\tAllele\n")

    for seq in sorted(alleles):
        out.write(f"{seq}\t{alleles[seq]}\n")

# ----------------------------
# Summary
# ----------------------------

print("\nAllele inference complete")
print("Bootstrap cutoff:", bootstrap_cutoff)
print("Alleles defined:", len(set(alleles.values())))
print("Output:", out_file)
