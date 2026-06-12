#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
import sys

# ----------------------------
# User input
# ----------------------------
fasta_in = input("Enter exon DNA FASTA file: ").strip()

hap_fasta_out = "SRK_exon_haplotypes.fasta"
key_out = "SRK_exon_haplotype_key.tsv"

# ----------------------------
# Containers
# ----------------------------
haplotypes = OrderedDict()   # dna_sequence -> hap_name
sequence_to_hap = {}         # original_id -> hap_name

hap_counter = 0
discarded_stop = 0
discarded_frame = 0

# ----------------------------
# Process sequences
# ----------------------------
for record in SeqIO.parse(fasta_in, "fasta"):
    dna_seq = str(record.seq).upper()

    # Length must be divisible by 3
    if len(dna_seq) % 3 != 0:
        discarded_frame += 1
        continue

    protein = str(Seq(dna_seq).translate())

    # Discard if internal stop codon
    if "*" in protein[:-1]:
        discarded_stop += 1
        continue

    # Remove terminal stop if present
    protein = protein.rstrip("*")

    # Collapse by DNA sequence (not protein)
    if dna_seq not in haplotypes:
        hap_counter += 1
        hap_name = f"SRK_exon_haplotype_{hap_counter:03d}"
        haplotypes[dna_seq] = hap_name

    sequence_to_hap[record.id] = haplotypes[dna_seq]

# ----------------------------
# Write haplotype FASTA (DNA)
# ----------------------------
hap_records = []

for dna_seq, hap_name in haplotypes.items():
    hap_records.append(
        SeqRecord(
            Seq(dna_seq),
            id=hap_name,
            description=""
        )
    )

SeqIO.write(hap_records, hap_fasta_out, "fasta")

# ----------------------------
# Write key table
# ----------------------------
with open(key_out, "w") as out:
    out.write("Original_sequence_ID\tExon_haplotype\n")
    for seq_id, hap in sequence_to_hap.items():
        out.write(f"{seq_id}\t{hap}\n")

# ----------------------------
# Summary
# ----------------------------
print("\n✅ Exon haplotype collapsing complete")
print(f"Functional exon haplotypes: {hap_counter}")
print(f"Discarded (frameshift):     {discarded_frame}")
print(f"Discarded (stop codons):    {discarded_stop}")
print(f"Haplotype FASTA written to: {hap_fasta_out}")
print(f"Key table written to:       {key_out}")

