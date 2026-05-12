#!/usr/bin/env python3
"""
pad_representatives.py — Step 22a-pre helper

Pads the LEPA allele-representative FASTA (`SRK_protein_allele_representatives.fasta`,
output of Step 10a) to a uniform sequence length so that downstream `mafft --add`
produces a clean combined alignment in Step 22a.

Why padding is needed
---------------------
`mafft --add reference.fasta query.fasta` requires that all sequences in `query.fasta`
share the same length (because `--add` treats the query as a pre-aligned profile).
The Step 10a representatives are mostly the same length (typically the alignment
width minus any trailing gaps in the representative chosen for each allele), but a
small number can be one or two residues shorter when the chosen representative had
trailing alignment gaps stripped. Padding to the longest length with `-` (gap)
characters at the C-terminus brings every sequence to the same width without
altering biological content — the gaps are simply unaligned terminal positions
that `mafft --add` will absorb into the combined alignment.

Inputs
------
INPUT_FASTA  = "SRK_protein_allele_representatives.fasta"  (Step 10a output)

Output
------
OUTPUT_FASTA = "SRK_protein_allele_representatives_padded.fasta"
               — same headers and sequences, all right-padded to the maximum
                 length with '-'.

When to run
-----------
Run this script every time `define_SRK_alleles_from_distance.py` (Step 10a) is
re-executed with a new N_ALLELES value, or whenever the input dataset changes
(new library, re-run after a filter change). The output is the immediate input
to the `mafft --add` step that rebuilds `SRK_combined_alignment.fasta`.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

INPUT_FASTA = "SRK_protein_allele_representatives.fasta"
OUTPUT_FASTA = "SRK_protein_allele_representatives_padded.fasta"
PAD_CHAR = "-"

records = list(SeqIO.parse(INPUT_FASTA, "fasta"))
max_len = max(len(r.seq) for r in records)

padded = []
n_padded = 0
for r in records:
    s = str(r.seq)
    if len(s) < max_len:
        s = s + PAD_CHAR * (max_len - len(s))
        n_padded += 1
    padded.append(SeqRecord(Seq(s), id=r.id, description=r.description.split(None, 1)[1] if " " in r.description else ""))

SeqIO.write(padded, OUTPUT_FASTA, "fasta")

print(f"Input  : {INPUT_FASTA}  ({len(records)} sequences)")
print(f"Output : {OUTPUT_FASTA}")
print(f"Target length: {max_len} aa")
print(f"Padded {n_padded} sequence(s) with '{PAD_CHAR}' at C-terminus")
