#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq

# ----------------------------
# User input
# ----------------------------
msa_fasta = input("Enter MSA FASTA file: ").strip()
ref_id = input("Enter canonical reference ID (exact FASTA header): ").strip()
WINDOW = 25  # bp to backfill from reference

# ----------------------------
# Load MSA
# ----------------------------
records = list(SeqIO.parse(msa_fasta, "fasta"))
seqs = {r.id: list(str(r.seq)) for r in records}

if ref_id not in seqs:
    raise ValueError("Reference ID not found in MSA")

ref_seq = seqs[ref_id]
aln_len = len(ref_seq)

# ----------------------------
# Helper: find real sequence bounds
# ----------------------------
def real_bounds(seq):
    start = next((i for i, x in enumerate(seq) if x != '-'), None)
    end = next((i for i in range(len(seq)-1, -1, -1) if seq[i] != '-'), None)
    return start, end

# ----------------------------
# Process sequences
# ----------------------------
new_records = []

for sid, seq in seqs.items():
    new_seq = seq.copy()
    start, end = real_bounds(seq)

    if start is None:
        continue  # skip empty sequence

    # ---- Backfill first 25 bp ----
    for i in range(min(WINDOW, aln_len)):
        if new_seq[i] == '-':
            new_seq[i] = ref_seq[i]

    # ---- Backfill last 25 bp ----
    for i in range(max(0, aln_len - WINDOW), aln_len):
        if new_seq[i] == '-':
            new_seq[i] = ref_seq[i]

    # ---- Replace extra leading gaps with Ns ----
    for i in range(0, start):
        if new_seq[i] == '-':
            new_seq[i] = 'N'

    # ---- Replace extra trailing gaps with Ns ----
    for i in range(end + 1, aln_len):
        if new_seq[i] == '-':
            new_seq[i] = 'N'

    new_records.append(
        SeqIO.SeqRecord(
            Seq("".join(new_seq)),
            id=sid,
            description=""
        )
    )

# ----------------------------
# Write output
# ----------------------------
out_fasta = msa_fasta.replace(".fasta", "_backfilled.fasta")
SeqIO.write(new_records, out_fasta, "fasta")

print("\n✅ Done!")
print(f"Backfilled alignment written to: {out_fasta}")

