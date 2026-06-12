#!/usr/bin/env python3

from Bio import SeqIO
import csv
from collections import defaultdict
import math

# ----------------------------
# User input
# ----------------------------
genotype_tsv = input("Enter genotype TSV (from genotype_individuals_from_haplotypes.py): ").strip()
hap_fasta = input("Enter protein haplotype FASTA (aligned): ").strip()

out_dist = genotype_tsv.replace(".tsv", "_protein_distance.tsv")

# ----------------------------
# Load haplotype sequences
# ----------------------------
hap_seqs = {}

for rec in SeqIO.parse(hap_fasta, "fasta"):
    hap_seqs[rec.id] = str(rec.seq)

hap_ids = list(hap_seqs.keys())

# Check alignment consistency
seq_len = len(next(iter(hap_seqs.values())))
for s in hap_seqs.values():
    if len(s) != seq_len:
        raise ValueError("Haplotype FASTA is not aligned")

# ----------------------------
# Pairwise haplotype distances
# ----------------------------
def hap_distance(seq1, seq2):
    mismatches = 0
    compared = 0

    for a, b in zip(seq1, seq2):
        if a in "-X" or b in "-X":
            continue
        compared += 1
        if a != b:
            mismatches += 1

    if compared == 0:
        return math.nan

    return mismatches / compared

hap_dist = {}

for i in range(len(hap_ids)):
    for j in range(i, len(hap_ids)):
        h1 = hap_ids[i]
        h2 = hap_ids[j]
        d = hap_distance(hap_seqs[h1], hap_seqs[h2])
        hap_dist[(h1, h2)] = d
        hap_dist[(h2, h1)] = d

# ----------------------------
# Load individual genotypes
# ----------------------------
individuals = defaultdict(list)

with open(genotype_tsv) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        ind = row["Individual"]
        haps = row["Haplotypes"].split(",")
        individuals[ind] = haps

ind_ids = sorted(individuals)

# ----------------------------
# Individual–individual distances
# ----------------------------
def individual_distance(hapsA, hapsB):
    mins = []

    for hA in hapsA:
        dists = [
            hap_dist[(hA, hB)]
            for hB in hapsB
            if not math.isnan(hap_dist[(hA, hB)])
        ]
        if dists:
            mins.append(min(dists))

    if not mins:
        return math.nan

    return sum(mins) / len(mins)

# ----------------------------
# Write distance matrix
# ----------------------------
with open(out_dist, "w") as out:
    out.write("Individual\t" + "\t".join(ind_ids) + "\n")

    for i in ind_ids:
        row = []
        for j in ind_ids:
            d = individual_distance(individuals[i], individuals[j])
            row.append(f"{d:.4f}" if not math.isnan(d) else "NA")
        out.write(i + "\t" + "\t".join(row) + "\n")

# ----------------------------
# Summary
# ----------------------------
print("\n✅ Protein-based genetic distance complete")
print(f"Individuals: {len(ind_ids)}")
print(f"Haplotypes: {len(hap_ids)}")
print(f"Distance matrix written to: {out_dist}")

