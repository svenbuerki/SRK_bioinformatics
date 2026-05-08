#!/usr/bin/env python3
"""
srk_fetch_reference_alleles.py

Fetch multi-allele SRK protein references from NCBI for cross-species
variability comparison:
  - Brassica rapa SRK S-haplotypes
  - Brassica oleracea SRK S-haplotypes
  - Arabidopsis lyrata SRK S-haplotypes

Strategy: NCBI Entrez esearch on protein database, filter by length and
deduplicate by source (one representative per S-haplotype where possible).

Outputs:
  brassica_SRK_alleles.fasta            (B. rapa + B. oleracea, S-haplotype-deduplicated)
  arabidopsis_lyrata_SRK_alleles.fasta  (A. lyrata, S-haplotype-deduplicated)
  reference_SRK_fetch_log.tsv           accession, species, length, S-haplotype, title
"""

from Bio import Entrez, SeqIO
import time
import re
import csv
import sys
import os

Entrez.email = "svenbuerki@boisestate.edu"

# Target lengths for full-length SRK (ectodomain + TM + kinase ≈ 800-870 aa).
# We accept 700+ to allow some partial sequences that still have full S-domain.
MIN_LENGTH = 700
MAX_LENGTH = 900

# Maximum sequences to fetch per query (we'll deduplicate after)
MAX_PER_QUERY = 200

QUERIES = {
    "Brassica_rapa":      'Brassica rapa[Organism] AND ("S-receptor kinase"[Title] OR "SRK"[Title] OR "S-locus receptor"[Title])',
    "Brassica_oleracea":  'Brassica oleracea[Organism] AND ("S-receptor kinase"[Title] OR "SRK"[Title] OR "S-locus receptor"[Title] OR "S-domain"[Title] OR "self-incompatibility"[Title])',
    "Arabidopsis_lyrata": 'Arabidopsis lyrata[Organism] AND ("S-receptor kinase"[Title] OR "S-locus receptor"[Title] OR "SRK"[Title] OR "self-incompatibility"[Title])',
    "Arabidopsis_halleri": 'Arabidopsis halleri[Organism] AND ("S-receptor kinase"[Title] OR "S-locus receptor"[Title] OR "SRK"[Title] OR "self-incompatibility"[Title])',
}

# Per-species output buckets
# A. lyrata + A. halleri are pooled because they share the same SRK system
# and most published S-haplotypes come from a mix of the two species
SPECIES_TO_FILE = {
    "Brassica_rapa":       "brassica_SRK_alleles.fasta",
    "Brassica_oleracea":   "brassica_SRK_alleles.fasta",       # combined Brassica
    "Arabidopsis_lyrata":  "arabidopsis_SRK_alleles.fasta",
    "Arabidopsis_halleri": "arabidopsis_SRK_alleles.fasta",    # combined Arabidopsis
}

# Max representatives per species (to avoid bias from over-sampled species)
MAX_PER_SPECIES = 20

LOG_TSV = "reference_SRK_fetch_log.tsv"


def s_haplotype_label(title: str) -> str:
    """
    Extract a putative S-haplotype label (e.g. 'S9', 'BoS-29', 'BrS9') from the
    NCBI title. Returns 'unknown' if no S-haplotype tag is parseable.
    """
    # Common patterns: 'S-9', 'S9', 'S 9', 'BrS9', 'haplotype 9', 'allele S9'
    m = re.search(r"\bS[- ]?(\d+)\b", title)
    if m:
        return f"S{m.group(1)}"
    m = re.search(r"haplotype[- ](\w+)", title, re.IGNORECASE)
    if m:
        return f"H_{m.group(1)}"
    return "unknown"


def fetch_one_species(query: str, species: str):
    print(f"\n=== {species} ===")
    print(f"Query: {query}")
    handle = Entrez.esearch(db="protein", term=query, retmax=MAX_PER_QUERY)
    rec = Entrez.read(handle)
    handle.close()
    ids = rec["IdList"]
    print(f"  esearch returned {len(ids)} IDs")
    if not ids:
        return []

    # Fetch summaries to get title + length without paying full FASTA cost yet
    handle = Entrez.esummary(db="protein", id=",".join(ids))
    summaries = Entrez.read(handle)
    handle.close()

    candidates = []
    for s in summaries:
        title = s.get("Title", "")
        length = int(s.get("Length", 0))
        acc = s.get("AccessionVersion", s.get("Caption", ""))
        if not (MIN_LENGTH <= length <= MAX_LENGTH):
            continue
        candidates.append({
            "accession": acc,
            "length":    length,
            "title":     title,
            "species":   species,
            "haplotype": s_haplotype_label(title),
        })

    print(f"  {len(candidates)} candidates after length filter "
          f"({MIN_LENGTH}-{MAX_LENGTH} aa)")

    # Deduplicate by S-haplotype: keep the longest entry per haplotype
    by_h = {}
    for c in candidates:
        key = c["haplotype"]
        if key == "unknown":
            # keep all unknowns separately (use accession as tiebreaker)
            key = f"unknown_{c['accession']}"
        if key not in by_h or c["length"] > by_h[key]["length"]:
            by_h[key] = c
    deduped = list(by_h.values())

    # Sort: known haplotypes first (by S-number), then unknowns
    def sort_key(c):
        h = c["haplotype"]
        if h.startswith("S") and h[1:].isdigit():
            return (0, int(h[1:]))
        return (1, h)
    deduped.sort(key=sort_key)

    # Keep at most MAX_PER_SPECIES
    deduped = deduped[:MAX_PER_SPECIES]
    print(f"  {len(deduped)} representatives after S-haplotype dedup "
          f"(cap = {MAX_PER_SPECIES})")
    for c in deduped:
        print(f"    {c['accession']:18s} {c['length']:>4d} aa  "
              f"{c['haplotype']:>10s}  {c['title'][:80]}")

    return deduped


def fetch_fasta(accessions):
    """Fetch FASTA records for a list of accessions."""
    if not accessions:
        return []
    handle = Entrez.efetch(db="protein", id=",".join(accessions),
                           rettype="fasta", retmode="text")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    return records


def main():
    all_log = []
    species_to_records = {}
    for species, query in QUERIES.items():
        try:
            cands = fetch_one_species(query, species)
        except Exception as e:
            print(f"  ERROR fetching {species}: {e}", file=sys.stderr)
            continue
        time.sleep(0.4)
        recs = fetch_fasta([c["accession"] for c in cands])
        # Re-tag each record header with a clean identifier
        relabeled = []
        for c, r in zip(cands, recs):
            new_id = f"{species}_{c['haplotype']}_{c['accession']}"
            r.id = new_id
            r.name = new_id
            r.description = (f"{species} {c['haplotype']} "
                             f"{c['accession']} length={len(r.seq)}")
            relabeled.append(r)
            all_log.append({**c, "fetched_length": len(r.seq)})
        species_to_records[species] = relabeled
        time.sleep(0.4)

    # Write per-file FASTAs (Brassica species combined into one file)
    bucketed = {}
    for species, recs in species_to_records.items():
        outfile = SPECIES_TO_FILE[species]
        bucketed.setdefault(outfile, []).extend(recs)
    for outfile, recs in bucketed.items():
        SeqIO.write(recs, outfile, "fasta")
        print(f"\nWritten {outfile}: {len(recs)} sequences")

    with open(LOG_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["species", "haplotype", "accession",
                                       "length", "fetched_length", "title"])
        w.writeheader()
        for r in all_log:
            w.writerow(r)
    print(f"Written {LOG_TSV}: {len(all_log)} entries")


if __name__ == "__main__":
    main()
