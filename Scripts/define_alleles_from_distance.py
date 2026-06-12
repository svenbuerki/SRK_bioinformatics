#!/usr/bin/env python3

import sys
from Bio import AlignIO
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


# ----------------------------
# User settings
# ----------------------------

FASTA = "SRK_exon_haplotypes.aligned.fasta"

DIST_THRESHOLD = 0.01   # 3% distance cutoff

OUT_TABLE = "SRK_distance_defined_alleles.tsv"

OUT_FASTA = "SRK_distance_defined_alleles.fasta"


# ----------------------------
# Read alignment
# ----------------------------

alignment = AlignIO.read(FASTA, "fasta")

names = [record.id for record in alignment]

seqs = [str(record.seq).upper() for record in alignment]

n = len(seqs)


# ----------------------------
# Compute distance matrix
# ----------------------------

def genetic_distance(seq1, seq2):

    differences = 0
    valid = 0

    for a, b in zip(seq1, seq2):

        if a == "-" or b == "-":
            continue

        valid += 1

        if a != b:
            differences += 1

    return differences / valid


dist_matrix = np.zeros((n, n))

for i in range(n):

    for j in range(i+1, n):

        d = genetic_distance(seqs[i], seqs[j])

        dist_matrix[i,j] = d
        dist_matrix[j,i] = d


# ----------------------------
# Hierarchical clustering
# ----------------------------

condensed = squareform(dist_matrix)

Z = linkage(condensed, method="average")

clusters = fcluster(Z, DIST_THRESHOLD, criterion="distance")


# ----------------------------
# Assign allele names
# ----------------------------

allele_names = {}

for i, cluster in enumerate(clusters):

    allele = f"Allele_{cluster:03d}"

    allele_names[names[i]] = allele


# ----------------------------
# Write table
# ----------------------------

with open(OUT_TABLE, "w") as out:

    out.write("Sequence\tAllele\n")

    for name in names:

        out.write(f"{name}\t{allele_names[name]}\n")


# ----------------------------
# Write representative fasta
# ----------------------------

seen = {}

with open(OUT_FASTA, "w") as out:

    for record in alignment:

        allele = allele_names[record.id]

        if allele not in seen:

            seen[allele] = True

            out.write(f">{allele}\n{record.seq}\n")


# ----------------------------
# Summary
# ----------------------------

print()

print("Allele definition complete")

print("Sequences:", n)

print("Alleles:", len(set(clusters)))

print("Output table:", OUT_TABLE)

print("Output fasta:", OUT_FASTA)
