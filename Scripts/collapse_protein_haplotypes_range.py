#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict

# ----------------------------
# User input
# ----------------------------
aa_fasta = input("Enter aligned protein FASTA file: ").strip()

start_pos = int(input("Enter start AA position (1-based): ").strip())
end_pos = int(input("Enter end AA position: ").strip())

hap_fasta_out = f"SRK_protein_haplotypes_AA{start_pos}_{end_pos}.fasta"
key_out = f"SRK_protein_haplotype_key_AA{start_pos}_{end_pos}.tsv"

# ----------------------------
# Collapse identical haplotypes
#   - X is treated as missing data
#   - sequences differing only by X are identical
# ----------------------------
haplotypes = OrderedDict()      # masked_sequence -> haplotype name
representative_seq = {}         # masked_sequence -> original sequence
sequence_to_hap = {}            # original ID -> haplotype name

hap_counter = 0

for record in SeqIO.parse(aa_fasta, "fasta"):
    full_seq = str(record.seq)

    # Slice alignment (1-based → 0-based)
    sliced_seq = full_seq[start_pos - 1:end_pos]

    # Mask X as wildcard for comparison
    masked_seq = sliced_seq.replace("X", "*")

    if masked_seq not in haplotypes:
        hap_counter += 1
        hap_name = f"SRK_prot_haplotype_{hap_counter:03d}"
        haplotypes[masked_seq] = hap_name
        representative_seq[masked_seq] = sliced_seq

    sequence_to_hap[record.id] = haplotypes[masked_seq]

# ----------------------------
# Write haplotype FASTA
# ----------------------------
hap_records = []

for masked_seq, hap_name in haplotypes.items():
    hap_records.append(
        SeqRecord(
            Seq(representative_seq[masked_seq]),
            id=hap_name,
            description=f"AA_range={start_pos}-{end_pos}"
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
print(f"AA range analysed: {start_pos}-{end_pos}")
print(f"Unique protein haplotypes (X treated as missing): {hap_counter}")
print(f"Haplotype FASTA written to: {hap_fasta_out}")
print(f"Haplotype key written to: {key_out}")

