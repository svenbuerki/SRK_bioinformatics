#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

# ----------------------------
# User input
# ----------------------------

dna_fasta = "SRK_exon_haplotypes.aligned.fasta"
protein_fasta = "SRK_exon_haplotypes.translated.fasta"
allele_key_out = "SRK_protein_functional_alleles.tsv"

# Translation settings
GENETIC_CODE = 1   # Standard genetic code
STOP_SYMBOL = "*"

# ----------------------------
# Load DNA sequences
# ----------------------------

records = list(SeqIO.parse(dna_fasta, "fasta"))

protein_seqs = {}
dna_to_protein = {}

# ----------------------------
# Translate DNA → protein
# ----------------------------

for rec in records:
    seq_id = rec.id
    
    # Remove gaps from alignment before translation
    dna_seq = str(rec.seq).replace("-", "").upper()
    
    # Ensure length is multiple of 3
    remainder = len(dna_seq) % 3
    if remainder != 0:
        dna_seq = dna_seq[:len(dna_seq) - remainder]
    
    # Translate
    protein = str(Seq(dna_seq).translate(table=GENETIC_CODE))
    
    # Optional: remove terminal stop codon
    if protein.endswith(STOP_SYMBOL):
        protein = protein[:-1]
    
    protein_seqs[seq_id] = protein
    dna_to_protein[seq_id] = protein

# ----------------------------
# Define functional alleles
# ----------------------------

protein_to_allele = {}
allele_id_counter = 1

for protein in sorted(set(protein_seqs.values())):
    allele_name = f"FuncAllele_{allele_id_counter}"
    protein_to_allele[protein] = allele_name
    allele_id_counter += 1

# ----------------------------
# Write protein FASTA
# ----------------------------

with open(protein_fasta, "w") as out:
    for seq_id, protein in protein_seqs.items():
        allele = protein_to_allele[protein]
        out.write(f">{seq_id}|{allele}\n")
        out.write(f"{protein}\n")

# ----------------------------
# Write allele key TSV
# ----------------------------

with open(allele_key_out, "w") as out:
    out.write("DNA_Haplotype_ID\tProtein_sequence\tFunctional_allele\n")
    
    for seq_id in sorted(dna_to_protein):
        protein = dna_to_protein[seq_id]
        allele = protein_to_allele[protein]
        out.write(f"{seq_id}\t{protein}\t{allele}\n")

# ----------------------------
# Summary
# ----------------------------

print("\n✅ Functional allele definition complete")
print(f"DNA haplotypes: {len(dna_to_protein)}")
print(f"Functional protein alleles: {len(set(protein_to_allele.values()))}")
print(f"Protein FASTA: {protein_fasta}")
print(f"Allele key TSV: {allele_key_out}")

