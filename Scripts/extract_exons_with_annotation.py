#!/usr/bin/env python3
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ----------------------------
# User inputs
# ----------------------------
msa_file = input("Enter MSA FASTA path (aligned sequences, including canonical reference): ").strip()
canonical_id = input("Enter canonical reference ID (exact header in MSA): ").strip()
augustus_csv = input("Enter AUGUSTUS CSV with exon positions: ").strip()
strand = input("Enter strand of annotation for canonical reference (+/-): ").strip()
min_exon_len = input("Enter minimum exon length (default 0): ").strip()

min_exon_len = int(min_exon_len) if min_exon_len else 0

# ----------------------------
# Load MSA
# ----------------------------
msa = SeqIO.to_dict(SeqIO.parse(msa_file, "fasta"))

if canonical_id not in msa:
    raise ValueError("Canonical reference ID not found in MSA")

canonical_seq = msa[canonical_id].seq

# ----------------------------
# Read AUGUSTUS CSV
# ----------------------------
with open(augustus_csv, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    exons = [row for row in reader if row['V3'] == 'exon']

# Filter by minimum exon length and sort by start
exons = [e for e in exons if int(e['V5']) - int(e['V4']) + 1 >= min_exon_len]
exons.sort(key=lambda x: int(x['V4']))

# ----------------------------
# Build genomic -> MSA mapping
# ----------------------------
genomic_to_msa = {}
seq_index = 0
for i, nt in enumerate(canonical_seq):
    if nt != '-':
        genomic_to_msa[seq_index + 1] = i  # 1-based genomic positions
        seq_index += 1

# ----------------------------
# Extract exon sequences
# ----------------------------
exon_records = []

for seq_name, seq_record in msa.items():
    seq_aln = seq_record.seq
    exon_seq = Seq("")

    for e in exons:
        start_gen = int(e['V4'])
        end_gen = int(e['V5'])
        try:
            msa_start = genomic_to_msa[start_gen]
            msa_end = genomic_to_msa[end_gen]
        except KeyError:
            continue  # skip if exon position not in mapping

        exon_seq += seq_aln[msa_start:msa_end+1]

    if strand == '-':
        exon_seq = exon_seq.reverse_complement()

    exon_records.append(SeqRecord(exon_seq, id=seq_name, description=""))

# ----------------------------
# Write exon-only FASTA
# ----------------------------
out_file = msa_file.replace(".fasta", "_exons.fasta")
SeqIO.write(exon_records, out_file, "fasta")
print(f"Exon DNA FASTA written to: {out_file}")

