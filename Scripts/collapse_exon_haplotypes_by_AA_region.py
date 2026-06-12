#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict

# ----------------------------
# User input
# ----------------------------
fasta_in = input("Enter exon haplotype FASTA (DNA): ").strip()
aa_len = int(input("Enter number of amino acids to consider (e.g. 400): ").strip())

nt_len = aa_len * 3

fasta_out = fasta_in.replace(".fasta", f"_AA1_{aa_len}_haplotypes.fasta")
key_out = fasta_in.replace(".fasta", f"_AA1_{aa_len}_haplotype_key.tsv")

# ----------------------------
# Collapse haplotypes
# ----------------------------
haplotypes = OrderedDict()   # subseq -> hap_name
seq_to_hap = {}              # original_id -> hap_name
hap_counter = 0

for record in SeqIO.parse(fasta_in, "fasta"):
    seq = str(record.seq)

    if len(seq) < nt_len:
        print(f"WARNING: {record.id} shorter than {nt_len} nt, skipping")
        continue

    sub_seq = seq[:nt_len]

    if sub_seq not in haplotypes:
        hap_counter += 1
        hap_name = f"SRK_exon_haplotype_{hap_counter:03d}"
        haplotypes[sub_seq] = hap_name

    seq_to_hap[record.id] = haplotypes[sub_seq]

# ----------------------------
# Write haplotype FASTA
# ----------------------------
hap_records = [
    SeqRecord(Seq(seq), id=hap, description="")
    for seq, hap in haplotypes.items()
]

SeqIO.write(hap_records, fasta_out, "fasta")

# ----------------------------
# Write key table
# ----------------------------
with open(key_out, "w") as out:
    out.write("Original_sequence_ID\tExon_haplotype\n")
    for seq_id, hap in seq_to_hap.items():
        out.write(f"{seq_id}\t{hap}\n")

# ----------------------------
# Summary
# ----------------------------
print("\n✅ Exon haplotype collapsing complete")
print(f"Amino acids used: {aa_len}")
print(f"Nucleotides used: {nt_len}")
print(f"Unique exon haplotypes: {hap_counter}")
print(f"Haplotype FASTA: {fasta_out}")
print(f"Key table: {key_out}")

