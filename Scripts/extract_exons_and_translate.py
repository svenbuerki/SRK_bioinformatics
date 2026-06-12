#!/usr/bin/env python3
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# -----------------------------
# User input
# -----------------------------
haplo_fasta = input("Enter haplotype FASTA path (MSA, including canonical reference): ").strip()
ref_id = input("Enter canonical reference ID (exactly as in MSA headers): ").strip()
augustus_csv = input("Enter AUGUSTUS CSV file with exon positions: ").strip()
strand_input = input("Is the canonical annotation on negative strand? [y/N]: ").strip().lower()
reverse_complement = strand_input == 'y'

min_exon_len = input("Enter minimum exon length to keep (default 0): ").strip()
min_exon_len = int(min_exon_len) if min_exon_len else 0

# -----------------------------
# Load MSA
# -----------------------------
msa_records = list(SeqIO.parse(haplo_fasta, "fasta"))
ref_record = next((r for r in msa_records if r.id == ref_id), None)
if not ref_record:
    raise ValueError(f"Reference {ref_id} not found in MSA FASTA.")

# -----------------------------
# Map genomic positions to MSA indices
# -----------------------------
genomic_to_msa = {}
g_pos = 1  # assuming 1-based positions in AUGUSTUS
for i, nt in enumerate(ref_record.seq):
    if nt != '-':
        genomic_to_msa[g_pos] = i
        g_pos += 1

# -----------------------------
# Read AUGUSTUS CSV and collect exon coordinates
# -----------------------------
exon_coords = []
with open(augustus_csv, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['V3'].lower() == 'exon':
            start = int(row['V4'])
            end = int(row['V5'])
            exon_coords.append((start, end))

if not exon_coords:
    raise ValueError("No exons found in AUGUSTUS CSV.")

# -----------------------------
# Extract exons for each sequence in MSA
# -----------------------------
dna_records = []
aa_records = []

for rec in msa_records:
    if rec.id == ref_id:
        continue  # skip canonical reference

    full_seq = str(rec.seq)

    exon_seq_parts = []
    for start, end in exon_coords:
        # map genomic to MSA, ensure end nucleotide included
        try:
            msa_start = genomic_to_msa[start]
            msa_end = max([i for g, i in genomic_to_msa.items() if g <= end])
        except KeyError:
            continue  # skip if position not found in MSA

        part = full_seq[msa_start:msa_end + 1]
        exon_seq_parts.append(part)

    exon_seq = ''.join(exon_seq_parts).replace('-', '')  # remove gaps

    if reverse_complement:
        exon_seq = str(Seq(exon_seq).reverse_complement())

    if len(exon_seq) < min_exon_len:
        continue

    dna_records.append(SeqRecord(Seq(exon_seq), id=rec.id, description=""))

    # trim to multiple of 3 for translation
    trimmed_len = len(exon_seq) - (len(exon_seq) % 3)
    aa_seq = str(Seq(exon_seq[:trimmed_len]).translate(to_stop=False))
    aa_records.append(SeqRecord(Seq(aa_seq), id=rec.id, description=""))

# -----------------------------
# Write output
# -----------------------------
out_prefix = haplo_fasta.rsplit('.', 1)[0]
dna_outfile = f"{out_prefix}_aligned_exons.fasta"
aa_outfile = f"{out_prefix}_aligned_exons_aa.fasta"

SeqIO.write(dna_records, dna_outfile, "fasta")
SeqIO.write(aa_records, aa_outfile, "fasta")

print(f"✅ Extraction complete.")
print(f"DNA exons written to: {dna_outfile}")
print(f"Translated AA sequences written to: {aa_outfile}")

