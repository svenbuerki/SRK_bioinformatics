#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
import glob
import os
import sys

# ----------------------------
# User input
# ----------------------------
frame = int(input("Enter translation frame (1, 2, or 3): ").strip())
if frame not in (1, 2, 3):
    sys.exit("❌ Frame must be 1, 2, or 3")

frame -= 1  # convert to 0-based

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
    print("\n❌ No *_exons_backfilled.fasta files found.")
    print("   Current directory:", os.getcwd())
    sys.exit(1)

# ----------------------------
# Merge, filter, collapse
# ----------------------------
haplotypes = OrderedDict()   # DNA sequence -> haplotype name
sequence_to_hap = {}
hap_counter = 0
skipped = 0

for fasta in fasta_files:
    print(f"\n📥 Processing {fasta}")
    for record in SeqIO.parse(fasta, "fasta"):

        # 1. Remove alignment gaps
        dna_str = str(record.seq).replace("-", "").upper()

        # 2. Apply frame
        dna_str = dna_str[frame:]

        # 3. Trim to multiple of 3
        dna_str = dna_str[:len(dna_str) - (len(dna_str) % 3)]

        if len(dna_str) < 3:
            skipped += 1
            continue

        dna_seq = Seq(dna_str)

        # 4. Translate for QC only
        try:
            prot = str(dna_seq.translate(to_stop=False))
        except Exception:
            skipped += 1
            continue

        # 5. Reject if internal stop codons
        if "*" in prot[:-1]:
            skipped += 1
            continue

        # 6. Collapse DNA haplotypes
        if dna_str not in haplotypes:
            hap_counter += 1
            hap_name = f"SRK_exon_haplotype_{hap_counter:03d}"
            haplotypes[dna_str] = hap_name

        sequence_to_hap[record.id] = haplotypes[dna_str]

# ----------------------------
# Write haplotype FASTA (DNA)
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
print("\n✅ Functional exon haplotype collapsing complete")
print(f"Input FASTA files:        {len(fasta_files)}")
print(f"Functional haplotypes:    {hap_counter}")
print(f"Skipped (non-functional): {skipped}")
print(f"Haplotype FASTA:          {hap_fasta_out}")
print(f"Haplotype key:            {key_out}")

