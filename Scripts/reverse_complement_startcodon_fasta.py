#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ----------------------------
# User input
# ----------------------------
in_fasta = input("Enter input DNA FASTA file: ").strip()

out_fasta = in_fasta.replace(".fasta", "_revcomp_ATGfixed.fasta")
log_file = in_fasta.replace(".fasta", "_startcodon_fix.log")

# ----------------------------
# Process sequences
# ----------------------------
records_out = []

with open(log_file, "w") as log:
    log.write("Sequence_ID\tAction\tNew_start\n")

    for record in SeqIO.parse(in_fasta, "fasta"):
        seq_rc = record.seq.reverse_complement()
        seq_str = str(seq_rc)

        action = "OK"

        if seq_str.startswith("ATG"):
            pass
        elif seq_str.startswith("TG"):
            seq_str = "A" + seq_str
            action = "Added_A"
        elif seq_str.startswith("G"):
            seq_str = "AT" + seq_str
            action = "Added_AT"
        else:
            action = "No_ATG_found"

        log.write(f"{record.id}\t{action}\t{seq_str[:3]}\n")

        records_out.append(
            SeqRecord(
                Seq(seq_str),
                id=record.id,
                description=record.description
            )
        )

# ----------------------------
# Write output
# ----------------------------
SeqIO.write(records_out, out_fasta, "fasta")

print("\n✅ Reverse complement + start codon correction complete")
print(f"Corrected FASTA: {out_fasta}")
print(f"Log file:        {log_file}")

