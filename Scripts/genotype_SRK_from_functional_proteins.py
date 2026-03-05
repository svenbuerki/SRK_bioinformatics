#!/usr/bin/env python3

import csv
from collections import defaultdict, OrderedDict


# ----------------------------
# User input
# ----------------------------

protein_key_file = input(
    "Enter SRK_functional_protein_key.tsv: "
).strip()


long_out = "SRK_individual_protein_table.tsv"

wide_out = "SRK_individual_genotypes.tsv"



# ----------------------------
# Load protein assignments
# ----------------------------

individual_proteins = defaultdict(list)


with open(protein_key_file, newline="") as f:

    reader = csv.DictReader(f, delimiter="\t")

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

print("\n✅ Functional protein genotyping complete")

print("Individuals:", len(genotype_matrix))

print("Functional alleles:", len(all_proteins))

print("Outputs:")

print(long_out)

print(wide_out)
