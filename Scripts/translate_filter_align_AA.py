#!/usr/bin/env python3
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

# ----------------------------
# User input
# ----------------------------
dna_fasta = input("Enter aligned DNA FASTA file: ").strip()
frame = int(input("Enter translation frame (1, 2, or 3): ").strip())

if frame not in (1, 2, 3):
    raise ValueError("Frame must be 1, 2, or 3")

# ----------------------------
# Output files
# ----------------------------
aa_raw = dna_fasta.replace(".fasta", f"_frame{frame}_AA_raw.fasta")
aa_filtered = dna_fasta.replace(".fasta", f"_frame{frame}_AA_filtered.fasta")
log_file = dna_fasta.replace(".fasta", f"_frame{frame}_stopcodon_log.tsv")
aa_aligned = dna_fasta.replace(".fasta", f"_frame{frame}_AA_filtered_aligned.fasta")

GENETIC_CODE = 1  # plant nuclear code

# ----------------------------
# Translate DNA → AA
# ----------------------------
aa_records = []

for record in SeqIO.parse(dna_fasta, "fasta"):
    ungapped = str(record.seq).replace("-", "").upper()

    cds = ungapped[frame - 1 :]
    cds = cds[: len(cds) - (len(cds) % 3)]

    if len(cds) == 0:
        aa_seq = ""
    else:
        aa_seq = str(Seq(cds).translate(table=GENETIC_CODE, to_stop=False))

    aa_records.append(
        SeqIO.SeqRecord(
            Seq(aa_seq),
            id=record.id,
            description=f"frame={frame}"
        )
    )

SeqIO.write(aa_records, aa_raw, "fasta")

# ----------------------------
# Filter internal stop codons
# ----------------------------
filtered_records = []

with open(log_file, "w") as log:
    log.write("Sequence_ID\tAA_length\tStop_count\tStop_positions\tStatus\n")

    for rec in aa_records:
        seq = str(rec.seq)
        stop_positions = [i + 1 for i, aa in enumerate(seq[:-1]) if aa == "*"]
        stop_count = len(stop_positions)

        if stop_count == 0:
            status = "OK"
            filtered_records.append(rec)
        else:
            status = "REMOVED"

        log.write(
            f"{rec.id}\t{len(seq)}\t{stop_count}\t"
            f"{','.join(map(str, stop_positions))}\t{status}\n"
        )

SeqIO.write(filtered_records, aa_filtered, "fasta")

# ----------------------------
# Align AA sequences with MAFFT
# ----------------------------
try:
    subprocess.run(
        ["mafft", "--auto", aa_filtered],
        check=True,
        stdout=open(aa_aligned, "w"),
        stderr=subprocess.DEVNULL
    )
except FileNotFoundError:
    print("\n⚠ MAFFT not found. Skipping alignment.")
    aa_aligned = None

# ----------------------------
# Summary
# ----------------------------
print("\n✅ Translation and filtering complete")
print(f"Raw AA FASTA:        {aa_raw}")
print(f"Filtered AA FASTA:   {aa_filtered}")
print(f"Stop codon log:      {log_file}")

if aa_aligned:
    print(f"Aligned AA FASTA:    {aa_aligned}")

