#!/usr/bin/env python3

import csv
import os
import sys
from collections import defaultdict, OrderedDict


# ----------------------------
# Output files & overwrite check
# ----------------------------

long_out = "SRK_individual_protein_table.tsv"
wide_out = "SRK_individual_genotypes.tsv"

OUTPUT_FILES = [long_out, wide_out]

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
# User input
# ----------------------------

protein_key_file = input(
    "Enter SRK_functional_protein_key.tsv: "
).strip()



# ----------------------------
# Load protein assignments
# ----------------------------

individual_proteins = defaultdict(list)


with open(protein_key_file, newline="", encoding="utf-8-sig") as f:

    reader = csv.DictReader(f, delimiter="\t")

    # Validate expected columns before processing rows
    expected = {"Original_sequence_ID", "Protein"}
    if not expected.issubset(set(reader.fieldnames or [])):
        print(f"ERROR: Expected columns {expected}")
        print(f"       Found columns:    {reader.fieldnames}")
        print("Check that you provided the correct file (SRK_functional_protein_key.tsv).")
        sys.exit(1)

    for row in reader:

        seq_id = row["Original_sequence_ID"]

        protein = row["Protein"]


        parts = seq_id.split("_")


        # FIXED extraction
        if len(parts) < 2:

            print("WARNING: cannot parse:", seq_id)

            continue


        individual = f"{parts[0]}_{parts[1]}"


        individual_proteins[individual].append(protein)



print("Individuals found:", len(individual_proteins))



# ----------------------------
# Write long table
# ----------------------------

with open(long_out, "w") as out:

    out.write("Individual\tProtein\n")

    for individual in sorted(individual_proteins):

        for protein in individual_proteins[individual]:

            out.write(f"{individual}\t{protein}\n")



# ----------------------------
# Build genotype matrix
# ----------------------------

all_proteins = sorted({

    protein

    for proteins in individual_proteins.values()

    for protein in proteins

})


genotype_matrix = OrderedDict()


for individual in sorted(individual_proteins):

    protein_set = set(individual_proteins[individual])

    genotype_matrix[individual] = {

        protein: 1 if protein in protein_set else 0

        for protein in all_proteins

    }



# ----------------------------
# Write matrix
# ----------------------------

with open(wide_out, "w", newline="") as out:

    writer = csv.writer(out, delimiter="\t")

    header = ["Individual"] + all_proteins

    writer.writerow(header)


    for individual, geno in genotype_matrix.items():

        row = [individual] + [

            geno[p]

            for p in all_proteins

        ]

        writer.writerow(row)



# ----------------------------
# Summary
# ----------------------------

print("\nFINAL SUMMARY")
print(f"\nIndividuals genotyped: {len(genotype_matrix)}")
print(f"Functional proteins: {len(all_proteins)}")
print(f"\nOutputs:")
print(f"  {long_out}")
print(f"  {wide_out}")
