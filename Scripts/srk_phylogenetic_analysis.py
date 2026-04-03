#!/usr/bin/env python3

import csv
import os
import shutil
import subprocess
import sys
from Bio import Phylo, SeqIO

# ----------------------------
# Files
# ----------------------------
metadata_file   = "sampling_metadata.csv"
genotype_file   = "SRK_individual_genotypes.tsv"
alignment_file  = "SRK_functional_proteins_aligned.fasta"
filtered_alignment = "SRK_ingroup_proteins_aligned.fasta"
iqtree_prefix   = "SRK_phylogeny"

# IQ-TREE generates several output files sharing the same prefix
IQTREE_EXTENSIONS = ["treefile", "iqtree", "log", "model.gz", "contree", "splits.nex", "bionj", "mldist"]
iqtree_outputs = [f"{iqtree_prefix}.{ext}" for ext in IQTREE_EXTENSIONS]

class_file         = "SRK_protein_classes.tsv"
ind_class_file     = "SRK_individual_class_genotype.tsv"
rooted_treefile    = f"{iqtree_prefix}_midpoint_rooted.treefile"

OUTPUT_FILES = [filtered_alignment, class_file, ind_class_file,
                rooted_treefile] + iqtree_outputs

# ----------------------------
# Locate IQ-TREE3 executable
# ----------------------------
iqtree_exe = None
for candidate in ["./iqtree3", "./iqtree", "iqtree3", "iqtree"]:
    if candidate.startswith("./"):
        if os.path.isfile(candidate[2:]):
            iqtree_exe = candidate
            break
    elif shutil.which(candidate):
        iqtree_exe = candidate
        break

if iqtree_exe is None:
    user_path = input("IQ-TREE3 not found. Enter path to executable (e.g. ./iqtree3): ").strip()
    if not os.path.isfile(user_path.lstrip("./")):
        print(f"ERROR: File not found: {user_path}")
        sys.exit(1)
    iqtree_exe = user_path
print(f"Using IQ-TREE executable: {iqtree_exe}")

# ----------------------------
# Overwrite check
# ----------------------------
existing = [f for f in OUTPUT_FILES if os.path.exists(f)]
if existing:
    print("\nWARNING: The following output files already exist and will be overwritten:")
    for f in existing:
        print(f"  {f}")
    confirm = input("\nProceed and overwrite? (y/n): ").strip().lower()
    if confirm != 'y':
        print("Aborted. No files were modified.")
        sys.exit(0)

# ----------------------------
# Step 1: Load ingroup individuals from metadata
# ----------------------------
print("Loading metadata...")
ingroup_individuals = set()

with open(metadata_file, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row.get("Ingroup", "").strip() == "1":
            library = row["Library"].strip()
            barcode = row["barcode"].strip()
            # Normalise to Library001_barcode01 format regardless of whether
            # the metadata stores bare numbers (1, 2) or full strings (Library001, barcode01)
            if library.lstrip("-").isdigit():
                library = f"Library{int(library):03d}"
            if barcode.lstrip("-").isdigit():
                barcode = f"barcode{int(barcode):02d}"
            individual_id = f"{library}_{barcode}"
            ingroup_individuals.add(individual_id)

print(f"Ingroup individuals found in metadata: {len(ingroup_individuals)}")

if not ingroup_individuals:
    print("ERROR: No ingroup individuals found.")
    print("  Check that the 'Ingroup' column contains 1 for ingroup samples.")
    sys.exit(1)

# ----------------------------
# Step 2: Load genotype matrix and collect proteins present in ingroup
# ----------------------------
print("\nLoading genotype matrix...")
ingroup_proteins = set()
genotype_individuals = set()
matched_individuals = set()
ingroup_individual_proteins = {}   # individual -> set of proteins (for Step 5)

with open(genotype_file, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter="\t")
    protein_columns = [col for col in (reader.fieldnames or []) if col != "Individual"]

    for row in reader:
        individual = row["Individual"]
        genotype_individuals.add(individual)
        if individual in ingroup_individuals:
            matched_individuals.add(individual)
            proteins = {p for p in protein_columns if row.get(p, "0").strip() == "1"}
            ingroup_individual_proteins[individual] = proteins
            ingroup_proteins.update(proteins)

print(f"Individuals in genotype file:        {len(genotype_individuals)}")
print(f"Ingroup individuals matched:         {len(matched_individuals)}")
print(f"Proteins present in ingroup:         {len(ingroup_proteins)}")

# Diagnostic if some individuals from metadata are absent from the genotype file.
# This is expected for individuals excluded by the abundance or max_alleles filter
# in define_functional_proteins.py. Only investigate if the number seems too high.
unmatched = ingroup_individuals - genotype_individuals
if unmatched:
    print(f"\nNOTE: {len(unmatched)} ingroup individuals from metadata absent from genotype file.")
    print(f"  These were likely excluded by the abundance or max_alleles filter upstream.")
    print(f"  Absent individuals: {sorted(unmatched)[:5]}{' ...' if len(unmatched) > 5 else ''}")
    if len(unmatched) > len(ingroup_individuals) * 0.5:
        print(f"  WARNING: >50% of ingroup individuals are missing — check ID format.")
        print(f"  Example metadata IDs:  {sorted(ingroup_individuals)[:3]}")
        print(f"  Example genotype IDs:  {sorted(genotype_individuals)[:3]}")

if not ingroup_proteins:
    print("\nERROR: No proteins found for ingroup individuals. Cannot build tree.")
    sys.exit(1)

# ----------------------------
# Step 3: Filter alignment to ingroup proteins
# ----------------------------
print("\nFiltering alignment to ingroup proteins...")
all_records = list(SeqIO.parse(alignment_file, "fasta"))
filtered_records = [rec for rec in all_records if rec.id in ingroup_proteins]

print(f"Sequences in full alignment:         {len(all_records)}")
print(f"Sequences retained for ingroup:      {len(filtered_records)}")

# Warn if any expected proteins are absent from the alignment
missing_from_alignment = ingroup_proteins - {rec.id for rec in all_records}
if missing_from_alignment:
    print(f"\nWARNING: {len(missing_from_alignment)} ingroup proteins not found in alignment:")
    for p in sorted(missing_from_alignment):
        print(f"  {p}")

if not filtered_records:
    print("\nERROR: No sequences matched ingroup proteins in the alignment.")
    sys.exit(1)

if len(filtered_records) < 4:
    print(f"\nWARNING: Only {len(filtered_records)} sequences — tree may not be informative.")

SeqIO.write(filtered_records, filtered_alignment, "fasta")
print(f"\nFiltered alignment written: {filtered_alignment}")

# ----------------------------
# Step 4: Run IQ-TREE3
# ----------------------------
print("\nRunning IQ-TREE3...")

cmd = [
    iqtree_exe,
    "-s", filtered_alignment,
    "-m", "TEST",
    "-mset", "LG,WAG,JTT",
    "-mfreq", "FU",
    "-mrate", "G4,I,R4",
    "-bb", "1000",
    "-alrt", "1000",
    "-redo",
    "-nt", "AUTO",
    "-pre", iqtree_prefix
]

print("Command:", " ".join(cmd))
print()
subprocess.run(cmd, check=True)

# ----------------------------
# Step 5: Classify proteins into SRK Class I and Class II
# ----------------------------
# Strategy: find the longest internal branch in the unrooted tree.
# In SRK phylogenies the two classes are separated by an exceptionally long
# branch due to ancient balancing selection. The subtree connected by that
# branch = Class II (larger, recessive); everything else = Class I (smaller,
# dominant). This avoids the rooting problem and handles the Class I grade.

print("\nClassifying proteins into SRK classes...")
treefile = f"{iqtree_prefix}.treefile"
tree = Phylo.read(treefile, "newick")

# Find the non-terminal clade connected by the longest branch
max_branch_len = 0
best_clade = None

for clade in tree.find_clades(order="level"):
    if clade == tree.root or clade.is_terminal():
        continue
    if clade.branch_length and clade.branch_length > max_branch_len:
        max_branch_len = clade.branch_length
        best_clade = clade

if best_clade is None:
    print("WARNING: Could not identify a splitting branch. Check tree file.")
    class_I, class_II = set(), {t.name for t in tree.get_terminals()}
else:
    all_leaves  = {t.name for t in tree.get_terminals()}
    group_A     = {t.name for t in best_clade.get_terminals()}
    group_B     = all_leaves - group_A
    # Class I = smaller group (dominant), Class II = larger group (recessive)
    if len(group_A) <= len(group_B):
        class_I, class_II = group_A, group_B
    else:
        class_I, class_II = group_B, group_A

print(f"Longest internal branch length:  {max_branch_len:.4f}")
print(f"Class I  (dominant,  smaller):   {len(class_I)} proteins")
print(f"Class II (recessive, larger):    {len(class_II)} proteins")
print()
print("NOTE: Validate this split in FigTree using the midpoint-rooted tree.")
print("      For higher confidence, add published SRK reference sequences with")
print("      known class assignments to anchor the tree.")

# --- Write protein class assignments ---
with open(class_file, "w") as f:
    f.write("Protein\tClass\n")
    for p in sorted(class_I):
        f.write(f"{p}\tClass_I\n")
    for p in sorted(class_II):
        f.write(f"{p}\tClass_II\n")

# --- Write midpoint-rooted tree for visualisation in FigTree/iTOL ---
tree_viz = Phylo.read(treefile, "newick")
tree_viz.root_at_midpoint()
Phylo.write(tree_viz, rooted_treefile, "newick")

# --- Per-individual class genotype ---
# For each ingroup individual, list which Class I and Class II proteins they
# carry and derive a simple dominance prediction for crossing decisions.
with open(ind_class_file, "w") as f:
    f.write("Individual\t"
            "n_Class_I\tClass_I_proteins\t"
            "n_Class_II\tClass_II_proteins\t"
            "Dominance_prediction\n")
    for individual in sorted(ingroup_individual_proteins):
        proteins = ingroup_individual_proteins[individual]
        c1 = sorted(proteins & class_I)
        c2 = sorted(proteins & class_II)

        # Dominance rules:
        # - Carries Class I allele(s): expressed as dominant over Class II
        # - Carries only Class II: codominant with other Class II carriers
        # - No functional proteins: potentially self-compatible
        if c1 and c2:
            prediction = "Class_I_dominant"
        elif c1:
            prediction = "Class_I_only_dominant"
        elif c2:
            prediction = "Class_II_only_codominant"
        else:
            prediction = "no_functional_SRK"

        f.write(f"{individual}\t"
                f"{len(c1)}\t{';'.join(c1)}\t"
                f"{len(c2)}\t{';'.join(c2)}\t"
                f"{prediction}\n")

# ----------------------------
# Summary
# ----------------------------
print("\nFINAL SUMMARY")
print(f"\nIngroup individuals:     {len(matched_individuals)}")
print(f"Ingroup proteins:        {len(ingroup_proteins)}")
print(f"Sequences in tree:       {len(filtered_records)}")
print(f"\nClass I  (dominant):     {len(class_I)} proteins")
print(f"Class II (recessive):    {len(class_II)} proteins")
print(f"\nFiltered alignment:      {filtered_alignment}")
print(f"Tree file (unrooted):    {iqtree_prefix}.treefile")
print(f"Tree file (midpt rooted):{rooted_treefile}")
print(f"IQ-TREE log:             {iqtree_prefix}.iqtree")
print(f"Protein classes:         {class_file}")
print(f"Individual class geno:   {ind_class_file}")
