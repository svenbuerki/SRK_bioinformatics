#!/usr/bin/env python3
# Enhanced version with abundance filtering tracking

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter, defaultdict
import glob
import os
import subprocess
import sys
import re

# ----------------------------
# Files
# ----------------------------
# Two input naming conventions are accepted:
#   - Legacy:  *_exons_backfilled.fasta              (libraries pre-Step-7b)
#   - Current: *_exons_backfilled_Nfilt.fasta        (libraries with Step 7b applied)
# When both files exist for the same library, the _Nfilt version is used
# and the legacy file is shadowed (skipped).
patterns = ("*_exons_backfilled.fasta", "*_exons_backfilled_Nfilt.fasta")
raw_fasta = "SRK_proteins_raw.fasta"
aligned_fasta = "SRK_proteins_aligned.fasta"
final_fasta = "SRK_functional_proteins.fasta"
repr_aligned_fasta = "SRK_functional_proteins_aligned.fasta"
key_file = "SRK_functional_protein_key.tsv"
report_file = "SRK_individual_status_report.tsv"
sc_file = "SRK_self_compatible_candidates.txt"
too_many_file = "SRK_too_many_alleles.txt"

OUTPUT_FILES = [raw_fasta, aligned_fasta, final_fasta, repr_aligned_fasta,
                key_file, report_file, sc_file, too_many_file]

# ----------------------------
# Overwrite check
# ----------------------------
existing = [f for f in OUTPUT_FILES if os.path.exists(f)]
if existing:
    print("\nWARNING: The following output files already exist and will be overwritten:")
    for f in existing:
        print(f"  {f}")
    confirm = input("\nProceed and overwrite? (y/n): ").strip().lower()
    if confirm != 'y':
        print("Aborted. No files were modified.")
        sys.exit(0)

# ----------------------------
# Parameters
# ----------------------------
frame = int(input("Translation frame (1,2,3): ")) - 1
min_length = 100
min_count = int(input("Minimum sequences per protein (recommended 2-5): "))
max_alleles = int(input("Maximum functional proteins per individual (e.g. 4 for allotetraploid): "))

# Add option for enhanced validation
use_enhanced_validation = input("Use enhanced SRK domain validation? (y/n, default=n): ").strip().lower()
use_enhanced_validation = use_enhanced_validation == 'y'

# ----------------------------
# DEFINE FUNCTIONS FIRST
# ----------------------------
def validate_srk_features(protein_seq):
    """Enhanced SRK domain validation - OPTIONAL"""
    if not use_enhanced_validation:
        return True  # Skip enhanced validation, use original logic
    
    # Check for minimum cysteine content (SRK should have ~12)
    if protein_seq.count('C') < 6:
        return False
    
    # Check for basic kinase motifs
    kinase_motifs = [
        r'[LIVMFYC]G.{0,20}GxGxxG',  # ATP-binding
        r'[LIVMFYC]K.{6,20}D[LIVMFY]K',  # Catalytic
    ]
    
    motif_count = sum(1 for motif in kinase_motifs if re.search(motif, protein_seq))
    return motif_count >= 1

def generate_individual_report(individual_status, individuals_in_final_key, max_alleles):
    """Generate comprehensive individual status report with abundance filtering tracking"""
    
    report_data = []
    
    for individual, status in individual_status.items():
        total_seq = status['total_sequences']
        functional = status['functional_proteins']
        
        # Check if individual made it to final genotype data
        in_final_data = individual in individuals_in_final_key
        proteins_in_final = individuals_in_final_key.get(individual, 0)
        
        # Classify individual status
        if total_seq == 0:
            classification = "no_sequences"
        elif functional == 0:
            classification = "no_functional_proteins"
        elif not in_final_data and individual_status[individual].get('too_many_alleles', False):
            classification = "too_many_alleles"
        elif not in_final_data:
            classification = "dropped_abundance_filter"
        elif functional < total_seq * 0.5:
            classification = "low_functional_rate"
        else:
            classification = "normal"
        
        # Most common failure reasons
        failure_counts = Counter(status['failed_reasons'])
        top_failures = failure_counts.most_common(3)
        
        report_data.append({
            'Individual': individual,
            'Total_sequences': total_seq,
            'Functional_proteins': functional,
            'Functional_rate': functional/total_seq if total_seq > 0 else 0,
            'Proteins_in_final_data': proteins_in_final,  # NEW COLUMN
            'In_final_genotype_data': 'Yes' if in_final_data else 'No',  # NEW COLUMN
            'Classification': classification,
            'Top_failure_reasons': ';'.join([f"{reason}({count})" for reason, count in top_failures])
        })
    
    return report_data

# ----------------------------
# Enhanced tracking dictionaries
# ----------------------------
individual_status = defaultdict(lambda: {
    'total_sequences': 0,
    'functional_proteins': 0,
    'failed_reasons': [],
    'sequences_details': []
})

# ----------------------------
# Step 1: Translate and store mapping (ORIGINAL LOGIC + tracking)
# ----------------------------
print("\nTranslating sequences with enhanced tracking...")

raw_records = []
original_to_protein = {}

# Build input file list from both naming conventions, preferring _Nfilt
# over _backfilled when both exist for the same library.
nfilt_files = set(glob.glob(patterns[1]))
backfilled_files = set(glob.glob(patterns[0]))
shadowed = {f.replace("_Nfilt.fasta", ".fasta") for f in nfilt_files}
input_files = sorted(
    (nfilt_files | (backfilled_files - shadowed)) - set(OUTPUT_FILES)
)

if shadowed & backfilled_files:
    print("\nLegacy _backfilled.fasta files shadowed by _Nfilt versions (will be skipped):")
    for f in sorted(shadowed & backfilled_files):
        print(f"  {f}")

# Exclude known output files in case they accidentally match the input pattern
for fasta in input_files:
    print(f"Processing {fasta}...")
    
    for record in SeqIO.parse(fasta, "fasta"):
        # Extract individual ID from sequence name
        individual_id = "_".join(record.id.split('_')[:2])  # Adjust based on your naming
        
        individual_status[individual_id]['total_sequences'] += 1
        
        # ORIGINAL LOGIC FROM YOUR WORKING SCRIPT
        dna = str(record.seq).replace("-", "").upper()
        dna = dna[frame:]
        dna = dna[:len(dna)-(len(dna)%3)]
        
        failure_reasons = []
        prot = ""
        
        if len(dna) < 3:
            failure_reasons.append("too_short_dna")
            continue
        
        prot = str(Seq(dna).translate()).rstrip("*")
        
        # ORIGINAL VALIDATION CHECKS
        if len(prot) < min_length:
            failure_reasons.append("too_short_protein")
        
        if "X" in prot:
            failure_reasons.append("ambiguous_aa")
        
        if "*" in prot:
            failure_reasons.append("premature_stop")
        
        # OPTIONAL enhanced validation (only if requested)
        if use_enhanced_validation and not validate_srk_features(prot):
            failure_reasons.append("missing_srk_domains")
        
        # Store details for reporting
        seq_detail = {
            'sequence_id': record.id,
            'length_dna': len(dna),
            'length_protein': len(prot),
            'failure_reasons': failure_reasons
        }
        individual_status[individual_id]['sequences_details'].append(seq_detail)
        
        # ORIGINAL SUCCESS CRITERIA (same as your working script)
        if not failure_reasons:
            raw_records.append(
                SeqRecord(Seq(prot), id=record.id, description="")
            )
            original_to_protein[record.id] = prot
            individual_status[individual_id]['functional_proteins'] += 1
        else:
            individual_status[individual_id]['failed_reasons'].extend(failure_reasons)

print("Raw functional proteins:", len(raw_records))
SeqIO.write(raw_records, raw_fasta, "fasta")

# ----------------------------
# Step 2: Align (ORIGINAL)
# ----------------------------
print("\nRunning MAFFT alignment...")
subprocess.run(
    f"mafft --auto --amino {raw_fasta} > {aligned_fasta}",
    shell=True,
    check=True
)

# ----------------------------
# Step 3: Count aligned proteins (ORIGINAL)
# ----------------------------
aligned_records = list(SeqIO.parse(aligned_fasta, "fasta"))
aligned_seq_dict = {}
aligned_seqs = []

for rec in aligned_records:
    ungapped = str(rec.seq).replace("-", "")
    aligned_seq_dict[rec.id] = ungapped
    aligned_seqs.append(ungapped)

counts = Counter(aligned_seqs)
print("Unique aligned proteins:", len(counts))

# ----------------------------
# Step 4: Apply abundance filter (ORIGINAL)
# ----------------------------
filtered = {
    seq: count
    for seq, count in counts.items()
    if count >= min_count
}
print("After abundance filter:", len(filtered))

# ----------------------------
# Step 5: Assign protein names (ORIGINAL)
# ----------------------------
protein_names = {}
for i, seq in enumerate(filtered):
    protein_names[seq] = f"SRK_protein_{i+1:03d}"

# ----------------------------
# Step 6: Write final FASTA (ORIGINAL)
# ----------------------------
final_records = []
for seq, name in protein_names.items():
    count = filtered[seq]
    final_records.append(
        SeqRecord(
            Seq(seq),
            id=name,
            description=f"count={count}"
        )
    )

SeqIO.write(final_records, final_fasta, "fasta")

# ----------------------------
# Step 6b: Align representative sequences (one per protein) for phylogenetic analysis
# ----------------------------
print("\nAligning representative protein sequences for phylogenetic analysis...")
subprocess.run(
    f"mafft --auto --amino {final_fasta} > {repr_aligned_fasta}",
    shell=True,
    check=True
)
print(f"Representative alignment written: {repr_aligned_fasta}")

# ----------------------------
# Step 7: Build per-individual protein sets and filter by max_alleles
# ----------------------------

# Build per-individual set of distinct protein names that passed all filters
individual_to_proteins = defaultdict(set)
for seq_id, prot in original_to_protein.items():
    if prot in protein_names:
        individual_id = "_".join(seq_id.split('_')[:2])
        individual_to_proteins[individual_id].add(protein_names[prot])

# Identify overloaded individuals
overloaded_individuals = {
    ind for ind, proteins in individual_to_proteins.items()
    if len(proteins) > max_alleles
}

if overloaded_individuals:
    print(f"\nFiltering {len(overloaded_individuals)} individuals with >{max_alleles} distinct proteins...")
    for ind in sorted(overloaded_individuals):
        n = len(individual_to_proteins[ind])
        print(f"  {ind}: {n} proteins (excluded)")
        individual_status[ind]['too_many_alleles'] = True

# Write key excluding overloaded individuals
with open(key_file, "w") as out:
    out.write("Original_sequence_ID\tProtein\tCount\n")
    for seq_id, prot in original_to_protein.items():
        if prot in protein_names:
            individual_id = "_".join(seq_id.split('_')[:2])
            if individual_id not in overloaded_individuals:
                name = protein_names[prot]
                count = filtered[prot]
                out.write(f"{seq_id}\t{name}\t{count}\n")

# ----------------------------
# NEW: Track which individuals made it to final data
# ----------------------------
individuals_in_final_key = defaultdict(int)

# Count how many proteins each individual has in the final key file
for seq_id, prot in original_to_protein.items():
    if prot in protein_names:
        individual_id = "_".join(seq_id.split('_')[:2])
        if individual_id not in overloaded_individuals:
            individuals_in_final_key[individual_id] += 1

# ----------------------------
# NEW: Generate enhanced reports with abundance tracking
# ----------------------------
report_data = generate_individual_report(individual_status, individuals_in_final_key, max_alleles)

# Save individual status report with new columns
with open(report_file, "w") as f:
    f.write("Individual\tTotal_sequences\tFunctional_proteins\tFunctional_rate\t"
            "Proteins_in_final_data\tIn_final_genotype_data\tClassification\tTop_failure_reasons\n")
    for row in report_data:
        f.write(f"{row['Individual']}\t{row['Total_sequences']}\t{row['Functional_proteins']}\t"
                f"{row['Functional_rate']:.3f}\t{row['Proteins_in_final_data']}\t"
                f"{row['In_final_genotype_data']}\t{row['Classification']}\t{row['Top_failure_reasons']}\n")

# Identify potentially self-compatible individuals (updated categories)
sc_candidates = [row['Individual'] for row in report_data
                if row['Classification'] in ['no_functional_proteins', 'low_functional_rate', 'dropped_abundance_filter']]

with open(sc_file, "w") as f:
    f.write("# Individuals with no or low functional SRK proteins\n")
    f.write("# These may be self-compatible\n")
    for individual in sc_candidates:
        f.write(f"{individual}\n")

# Write separate file for individuals excluded by allele count filter
too_many = [row['Individual'] for row in report_data
            if row['Classification'] == 'too_many_alleles']

with open(too_many_file, "w") as f:
    f.write(f"# Individuals excluded: >{max_alleles} distinct functional proteins\n")
    f.write("# Likely reflects assembly artefacts, contamination, or mixed samples\n")
    for individual in too_many:
        n = len(individual_to_proteins.get(individual, set()))
        f.write(f"{individual}\t{n}_proteins\n")

# ----------------------------
# Final summary with abundance filtering stats
# ----------------------------
print("\nFINAL SUMMARY")
print(f"\nFinal functional proteins: {len(protein_names)}")
print(f"Protein FASTA: {final_fasta}")
print(f"Aligned representatives (phylogenetics): {repr_aligned_fasta}")
print(f"Protein key: {key_file}")

classification_counts = Counter([row['Classification'] for row in report_data])
print("\nIndividual classifications:")
for classification, count in sorted(classification_counts.items(), key=lambda x: -x[1]):
    print(f"  {classification}: {count}")

dropped_abundance = [row for row in report_data if row['Classification'] == 'dropped_abundance_filter']
print(f"\nIndividuals dropped by abundance filter (min_count={min_count}): {len(dropped_abundance)}")
if dropped_abundance:
    print("\nDropped individuals:")
    for row in dropped_abundance:
        print(f"  {row['Individual']}: {row['Functional_proteins']} functional proteins → 0 in final data")

too_many = [row['Individual'] for row in report_data if row['Classification'] == 'too_many_alleles']
print(f"\nIndividuals excluded by allele count filter (max_alleles={max_alleles}): {len(too_many)}")
if too_many:
    print("\nExcluded individuals:")
    for ind in too_many:
        n = len(individual_to_proteins.get(ind, set()))
        print(f"  {ind}: {n} distinct proteins (exceeds max {max_alleles})")

print(f"\nPotential self-compatible individuals: {len(sc_candidates)}")
print(f"\nDetailed report: {report_file}")
print(f"SC candidates: {sc_file}")
print(f"Too-many-alleles exclusions: {too_many_file}")

