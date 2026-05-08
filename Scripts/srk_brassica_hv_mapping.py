#!/usr/bin/env python3
"""
srk_brassica_hv_mapping.py

Build a Brassica-residue → LEPA-alignment-column mapping by reading the
mafft --add --keeplength alignment that contains B. rapa SRK9 alongside the
LEPA allele representatives. Outputs SRK_brassica_hv_mapping.tsv with the
12 SCR9-contact residues from Ma et al. 2016 (Cell Research) and the canonical
hvI / hvII / hvIII region spans inferred from those residues.

Reference:
  Ma R, Han Z, Hu Z, Lin G, Gong X, Zhang H, Nasrallah JB, Chai J. 2016.
  Structural basis for specific self-incompatibility response in Brassica.
  Cell Research 26: 1320-1329. (PDB 5GYY; B. rapa eSRK9-SCR9 complex)
"""

from Bio import SeqIO
import csv
import sys

ALN_EXPANDED = "SRK_protein_allele_representatives_with_brassica_full.fasta"
BR_RECORD_ID = "BrSRK9"
OUT_TSV      = "SRK_brassica_hv_mapping.tsv"

# Ma et al. 2016: SCR9-contact residues in B. rapa SRK9 (full-protein
# numbering, signal peptide residues 1-29 included). Each tuple is
# (residue_number, expected_aa, hv_region).
MA2016_CONTACTS = [
    (211, "V", "hvI"),
    (267, "F", "hvI"),
    (287, "P", "hvI"),
    (290, "F", "hvI"),
    (277, "L", "hvII"),
    (278, "V", "hvII"),
    (282, "P", "hvII"),
    (331, "H", "hvII"),
    (332, "T", "hvII"),
    (333, "R", "hvII"),
    (325, "E", "hvIII"),
    (330, "D", "hvIII"),
]

# Read the expanded mafft --add alignment (LEPA + Brassica, full length).
# We then derive a mapping from each *expanded* column back to its original
# LEPA-alignment column by skipping columns where every LEPA sequence is a
# gap (Brassica insertion columns).
records = list(SeqIO.parse(ALN_EXPANDED, "fasta"))
brassica_row = None
lepa_rows    = []
for rec in records:
    if rec.id == BR_RECORD_ID:
        brassica_row = str(rec.seq).upper()
    else:
        lepa_rows.append(str(rec.seq).upper())
if brassica_row is None:
    sys.exit(f"ERROR: '{BR_RECORD_ID}' not found in {ALN_EXPANDED}")
exp_len = len(brassica_row)

# Identify Brassica-insertion columns (every LEPA seq has '-' there)
expanded_to_lepa = [None] * exp_len   # None = Brassica-insertion column
lepa_col_idx = 0
for c in range(exp_len):
    is_insertion = all(row[c] == "-" for row in lepa_rows)
    if is_insertion:
        expanded_to_lepa[c] = None    # column does not exist in original LEPA alignment
    else:
        lepa_col_idx += 1             # 1-based original LEPA column
        expanded_to_lepa[c] = lepa_col_idx
n_lepa_cols = lepa_col_idx

# Build Brassica residue_number -> expanded_column (1-based) lookup
residue_to_expanded_col = {}
res_idx = 0
for c, aa in enumerate(brassica_row):
    if aa != "-":
        res_idx += 1
        residue_to_expanded_col[res_idx] = c + 1     # 1-based expanded column

n_brassica_residues = res_idx
print(f"BrSRK9 row: {n_brassica_residues} non-gap residues "
      f"in {exp_len}-column expanded alignment")
print(f"Original LEPA alignment columns recovered: {n_lepa_cols} "
      f"(matches expected 843)")

# Verify each contact residue and translate to LEPA-alignment column
rows = []
for resnum, expected_aa, hv in MA2016_CONTACTS:
    if resnum not in residue_to_expanded_col:
        print(f"  WARNING: residue {resnum} ({expected_aa}, {hv}) "
              f"missing from Brassica row")
        continue
    exp_col = residue_to_expanded_col[resnum]
    actual_aa = brassica_row[exp_col - 1]
    lepa_col = expanded_to_lepa[exp_col - 1]
    ok = "OK" if actual_aa == expected_aa else f"MISMATCH (saw {actual_aa})"
    rows.append({
        "Brassica_residue":      resnum,
        "Brassica_aa":           expected_aa,
        "HV_region":             hv,
        "Expanded_aln_col":      exp_col,
        "LEPA_aln_col":          lepa_col if lepa_col is not None else "NA(insertion)",
        "Verification":          ok,
    })

print(f"\n{'Residue':>8s}  {'AA':>3s}  {'HV':>5s}  {'ExpCol':>7s}  "
      f"{'LEPACol':>8s}  Verify")
for r in rows:
    print(f"{r['Brassica_residue']:>8d}  {r['Brassica_aa']:>3s}  "
          f"{r['HV_region']:>5s}  {r['Expanded_aln_col']:>7d}  "
          f"{str(r['LEPA_aln_col']):>8s}  {r['Verification']}")

with open(OUT_TSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(f"\nWritten {OUT_TSV}")

# Compute per-HV-region span (min/max LEPA alignment column) for plot bands
print("\n-- HV region spans in LEPA alignment coordinates --")
hv_spans = {}
for r in rows:
    if isinstance(r["LEPA_aln_col"], int):
        hv_spans.setdefault(r["HV_region"], []).append(r["LEPA_aln_col"])
for hv in ["hvI", "hvII", "hvIII"]:
    cols = sorted(hv_spans.get(hv, []))
    if cols:
        print(f"  {hv}: LEPA alignment columns {cols[0]}-{cols[-1]} "
              f"({len(cols)} contact residues)")
