#!/usr/bin/env python3
"""
genotype_SRK_from_alleles.py

Builds individual-level genotype matrices based on distance-defined S-allele
bins rather than raw functional protein sequences.

Data flow
---------
  SRK_functional_protein_key.tsv     Original_sequence_ID → Protein
        ↓  (parse individual from sequence ID)
  SRK_protein_allele_assignments.tsv  Protein → Allele
        ↓  (join + deduplicate)
  Per-individual allele set
        ↓
  SRK_individual_allele_table.tsv     long format  (Individual, Protein, Allele)
  SRK_individual_allele_genotypes.tsv wide matrix  (individuals × alleles, binary)

Reference sequences (IDs not starting with "Library") are excluded.

Usage
-----
    python genotype_SRK_from_alleles.py
"""

import csv
import os
import sys
from collections import defaultdict, OrderedDict

# ─────────────────────────────────────────────────────────────────────────────
# User settings
# ─────────────────────────────────────────────────────────────────────────────

PROTEIN_KEY_FILE  = "SRK_functional_protein_key.tsv"
ALLELE_TSV_FILE   = "SRK_protein_allele_assignments.tsv"

# Optional metadata file with an 'Ingroup' column (1 = ingroup, 0 = outgroup).
# Set to None to skip metadata filtering entirely.
METADATA_FILE     = "sampling_metadata.csv"

# If True (and METADATA_FILE is set), restrict genotyping to ingroup samples
# (Ingroup == 1). Set to False to include all Library* individuals regardless
# of ingroup status.
INGROUP_ONLY      = True

LONG_OUT  = "SRK_individual_allele_table.tsv"
WIDE_OUT  = "SRK_individual_allele_genotypes.tsv"

# Only process individuals whose parsed ID starts with this prefix.
# Reference sequences (e.g. SRK_BEA_*) are excluded.
INDIVIDUAL_PREFIX = "Library"

# ─────────────────────────────────────────────────────────────────────────────
# Overwrite check
# ─────────────────────────────────────────────────────────────────────────────

existing = [f for f in [LONG_OUT, WIDE_OUT] if os.path.exists(f)]
if existing:
    print("\nWARNING: The following output files already exist and will be overwritten:")
    for f in existing:
        print(f"  {f}")
    confirm = input("\nProceed and overwrite? (y/n): ").strip().lower()
    if confirm != "y":
        print("Aborted. No files were modified.")
        sys.exit(0)

# ─────────────────────────────────────────────────────────────────────────────
# Load metadata and build ingroup allow-list
# ─────────────────────────────────────────────────────────────────────────────

ingroup_ids = None   # None = no filtering applied

if METADATA_FILE is not None:
    if not os.path.exists(METADATA_FILE):
        print(f"WARNING: METADATA_FILE '{METADATA_FILE}' not found — skipping metadata filter.")
    else:
        ingroup_ids = set()
        n_ingroup = n_outgroup = 0
        with open(METADATA_FILE, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            if "SampleID" not in (reader.fieldnames or []) or \
               "Ingroup" not in (reader.fieldnames or []):
                print(f"ERROR: {METADATA_FILE} must have columns 'SampleID' and 'Ingroup'")
                sys.exit(1)
            for row in reader:
                if row["Ingroup"].strip() == "1":
                    ingroup_ids.add(row["SampleID"].strip())
                    n_ingroup += 1
                else:
                    n_outgroup += 1
        print(f"Metadata loaded: {n_ingroup} ingroup, {n_outgroup} outgroup samples")
        if INGROUP_ONLY:
            print(f"INGROUP_ONLY = True  →  restricting to {n_ingroup} ingroup samples")
        else:
            print(f"INGROUP_ONLY = False →  ingroup filter disabled, all Library* samples retained")
            ingroup_ids = None   # disable filter

# ─────────────────────────────────────────────────────────────────────────────
# Load allele assignments: Protein → Allele
# ─────────────────────────────────────────────────────────────────────────────

protein_to_allele = {}

with open(ALLELE_TSV_FILE, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter="\t")
    expected = {"Protein", "Allele"}
    if not expected.issubset(set(reader.fieldnames or [])):
        print(f"ERROR: {ALLELE_TSV_FILE} must have columns {expected}")
        print(f"       Found: {reader.fieldnames}")
        sys.exit(1)
    for row in reader:
        protein_to_allele[row["Protein"]] = row["Allele"]

print(f"Loaded {len(protein_to_allele)} protein → allele mappings")

# ─────────────────────────────────────────────────────────────────────────────
# Load functional protein key and map to individuals
# Individual is parsed as  parts[0]_parts[1]  from the sequence ID
# (e.g. Library001_barcode01 from Library001_barcode01_hap1_rep1_tig00000001)
# ─────────────────────────────────────────────────────────────────────────────

# individual → list of (protein, allele) tuples (before deduplication)
individual_data = defaultdict(list)
n_skipped_ref   = 0
n_skipped_nomap = 0

with open(PROTEIN_KEY_FILE, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter="\t")
    expected = {"Original_sequence_ID", "Protein"}
    if not expected.issubset(set(reader.fieldnames or [])):
        print(f"ERROR: {PROTEIN_KEY_FILE} must have columns {expected}")
        print(f"       Found: {reader.fieldnames}")
        sys.exit(1)

    for row in reader:
        seq_id  = row["Original_sequence_ID"]
        protein = row["Protein"]

        parts = seq_id.split("_")
        if len(parts) < 2:
            print(f"WARNING: cannot parse individual from: {seq_id}")
            continue

        individual = f"{parts[0]}_{parts[1]}"

        # Skip reference sequences
        if not individual.startswith(INDIVIDUAL_PREFIX):
            n_skipped_ref += 1
            continue

        # Skip outgroup samples if ingroup filter is active
        if ingroup_ids is not None and individual not in ingroup_ids:
            continue

        # Map protein to allele
        allele = protein_to_allele.get(protein)
        if allele is None:
            n_skipped_nomap += 1
            continue

        individual_data[individual].append((protein, allele))

print(f"Reference sequences excluded : {n_skipped_ref}")
print(f"Proteins with no allele match: {n_skipped_nomap}")
print(f"Individuals found            : {len(individual_data)}")

# ─────────────────────────────────────────────────────────────────────────────
# Deduplicate: keep unique (protein, allele) pairs per individual,
# then collapse to unique alleles for the genotype matrix
# ─────────────────────────────────────────────────────────────────────────────

# individual → sorted list of unique (protein, allele) pairs
individual_pairs = {
    ind: sorted(set(pairs), key=lambda x: (x[1], x[0]))
    for ind, pairs in individual_data.items()
}

# individual → set of unique alleles
individual_alleles = {
    ind: sorted({allele for _, allele in pairs})
    for ind, pairs in individual_pairs.items()
}

# ─────────────────────────────────────────────────────────────────────────────
# Write long table  (Individual | Protein | Allele)
# ─────────────────────────────────────────────────────────────────────────────

with open(LONG_OUT, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Individual", "Protein", "Allele"])
    for individual in sorted(individual_pairs):
        for protein, allele in individual_pairs[individual]:
            writer.writerow([individual, protein, allele])

print(f"\nLong table written to: {LONG_OUT}")

# ─────────────────────────────────────────────────────────────────────────────
# Build and write wide binary genotype matrix  (individuals × alleles)
# ─────────────────────────────────────────────────────────────────────────────

all_alleles = sorted({
    allele
    for alleles in individual_alleles.values()
    for allele in alleles
})

genotype_matrix = OrderedDict()
for individual in sorted(individual_alleles):
    allele_set = set(individual_alleles[individual])
    genotype_matrix[individual] = {a: 1 if a in allele_set else 0
                                   for a in all_alleles}

with open(WIDE_OUT, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Individual"] + all_alleles)
    for individual, geno in genotype_matrix.items():
        writer.writerow([individual] + [geno[a] for a in all_alleles])

print(f"Genotype matrix written to : {WIDE_OUT}")

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

allele_counts = [len(v) for v in individual_alleles.values()]

print("\n── FINAL SUMMARY ──")
print(f"Individuals genotyped : {len(genotype_matrix)}")
print(f"S-alleles in matrix   : {len(all_alleles)}")
print(f"Min alleles/individual: {min(allele_counts)}")
print(f"Max alleles/individual: {max(allele_counts)}")
print(f"Mean alleles/individual: {sum(allele_counts)/len(allele_counts):.2f}")

print("\nAllele carrier counts (individuals per allele):")
for allele in all_alleles:
    count = sum(1 for geno in genotype_matrix.values() if geno[allele] == 1)
    print(f"  {allele}: {count} individuals")
