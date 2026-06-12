#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict

# ----------------------------
# User input
# ----------------------------
aa_fasta = input("Enter aligned protein FASTA file: ").strip()

hap_fasta_out = "SRK_protein_haplotypes.fasta"
key_out = "SRK_protein_haplotype_key.tsv"

# ----------------------------
# Collapse identical sequences
# ----------------------------
haplotypes = OrderedDict()   # sequence -> haplotype name
sequence_to_hap = {}         # original ID -> haplotype name

hap_counter = 0

for record in SeqIO.parse(aa_fasta, "fasta"):
    seq_str = str(record.seq)

    if seq_str not in haplotypes:
        hap_counter += 1
        hap_name = f"SRK_prot_haplotype_{hap_counter:03d}"
        haplotypes[seq_str] = hap_name

    sequence_to_hap[record.id] = haplotypes[seq_str]

# ----------------------------
# Write haplotype FASTA
# ----------------------------
hap_records = []

for seq_str, hap_name in haplotypes.items():
    hap_records.append(
        SeqRecord(
            Seq(seq_str),
            id=hap_name,
            description=""
        )
    )

SeqIO.write(hap_records, hap_fasta_out, "fasta")

# ----------------------------
# Write key table
# ----------------------------
with open(key_out, "w") as out:
    out.write("Original_sequence_ID\tProtein_haplotype\n")
    for seq_id, hap in sequence_to_hap.items():
        out.write(f"{seq_id}\t{hap}\n")

# ----------------------------
# Summary
# ----------------------------
print("\n✅ Protein haplotype collapsing complete")
print(f"Unique protein haplotypes: {hap_counter}")
print(f"Haplotype FASTA written to: {hap_fasta_out}")
print(f"Key table written to: {key_out}")

