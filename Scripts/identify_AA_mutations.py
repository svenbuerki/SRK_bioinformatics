#!/usr/bin/env python3
from Bio import SeqIO
import pandas as pd

# -----------------------------
# User inputs
# -----------------------------
msa_file = input("Enter aligned protein FASTA file: ").strip()
output_file = input("Enter output CSV file name: ").strip()

# -----------------------------
# Read sequences
# -----------------------------
records = list(SeqIO.parse(msa_file, "fasta"))
hap_ids = [r.id for r in records]
seqs = [str(r.seq) for r in records]
seq_length = len(seqs[0])

# -----------------------------
# Identify variable positions
# -----------------------------
data = []

for pos in range(seq_length):
    aa_at_pos = [seq[pos] for seq in seqs if seq[pos] != 'X']  # Exclude missing data
    if len(set(aa_at_pos)) > 1:  # More than one unique AA, it's variable
        for i, seq in enumerate(seqs):
            aa = seq[pos]
            if aa != 'X':
                data.append({
                    "haplotype": hap_ids[i],
                    "position": pos + 1,  # 1-based position
                    "AA": aa
                })

# -----------------------------
# Save output
# -----------------------------
df = pd.DataFrame(data)
df.to_csv(output_file, index=False)
print(f"Variable AA positions written to: {output_file}")

