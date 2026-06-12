#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
import glob
import os
import sys

# ----------------------------
# Parameters
# ----------------------------
pattern = "*_exons_backfilled.fasta"

hap_fasta_out = "SRK_exon_haplotypes.fasta"
key_out = "SRK_exon_haplotype_key.tsv"

# ----------------------------
# Find input FASTA files
# ----------------------------
fasta_files = sorted(glob.glob(pattern))

print("🔍 Searching for files matching:", pattern)
for f in fasta_files:
    print("  →", f)

if not fasta_files:
    print("\n❌ No *_exons_backfilled.fasta files found in current directory.")
    print("   Current directory:", os.getcwd())
    sys.exit(1)

# ----------------------------
# Merge and collapse sequences
# ----------------------------
haplotypes = OrderedDict()    # sequence -> haplotype name
sequence_to_hap = {}          # original ID -> haplotype name
hap_counter = 0

for fasta in fasta_files:
    print(f"\n📥 Processing {fasta}")
    for record in SeqIO.parse(fasta, "fasta"):
        seq_str = str(record.seq)

        if seq_str not in haplotypes:
            hap_counter += 1
            hap_name = f"SRK_exon_haplotype_{hap_counter:03d}"
            haplotypes[seq_str] = hap_name

        sequence_to_hap[record.id] = haplotypes[seq_str]

# ----------------------------
# Write haplotype FASTA
# ----------------------------
hap_records = [
    SeqRecord(Seq(seq), id=hap, description="")
    for seq, hap in haplotypes.items()
]

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
print(f"Input FASTA files: {len(fasta_files)}")
print(f"Unique exon haplotypes: {hap_counter}")
print(f"Haplotype FASTA: {hap_fasta_out}")
print(f"Haplotype key:   {key_out}")

