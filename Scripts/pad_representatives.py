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
OUTPUT_FASTA = "Tables/Phase4/step22a_protein_allele_representatives_padded.fasta"
               — same headers and sequences, all right-padded to the maximum
                 length with '-'.

When to run
-----------
Run this script every time `define_SRK_alleles_from_distance.py` (Step 10a) is
re-executed with a new N_ALLELES value, or whenever the input dataset changes
(new library, re-run after a filter change). The output is the immediate input
to the `mafft --add` step that rebuilds `SRK_combined_alignment.fasta`.
"""

import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

INPUT_FASTA   = "Tables/Phase2/step10a_protein_allele_representatives.fasta"
ALLELE_TABLE  = "Tables/Phase2/step11_individual_allele_table.tsv"
OUTPUT_FASTA  = "Tables/Phase4/step22a_protein_allele_representatives_padded.fasta"
PAD_CHAR = "-"

# Restrict to alleles observed in at least one LEPA-ingroup individual (Step 11).
# Step 10a clusters from the full Step 9 protein pool, which includes outgroup
# species (L. montanum, L. freemontii, L. philonitron) — those carry their own
# allele clusters that would contaminate the LEPA variability landscape.
ingroup_alleles = set()
with open(ALLELE_TABLE, encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        a = row.get("Allele", "").strip()
        if a:
            ingroup_alleles.add(a)

records_all = list(SeqIO.parse(INPUT_FASTA, "fasta"))
records = [r for r in records_all if r.id in ingroup_alleles]
n_dropped = len(records_all) - len(records)
if n_dropped > 0:
    dropped = sorted({r.id for r in records_all} - ingroup_alleles)
    print(f"Dropped {n_dropped} outgroup-only allele(s): {', '.join(dropped)}")

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

print(f"Input  : {INPUT_FASTA}  ({len(records_all)} sequences)")
print(f"Output : {OUTPUT_FASTA}  ({len(records)} LEPA-ingroup sequences)")
print(f"Target length: {max_len} aa")
print(f"Padded {n_padded} sequence(s) with '{PAD_CHAR}' at C-terminus")
