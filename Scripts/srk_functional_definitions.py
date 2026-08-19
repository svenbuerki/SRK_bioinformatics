#!/usr/bin/env python3
"""Functional-SRK definition — cross-genera HV union + polymorphic linker
positions.

Motivation (2026-08-18): LEPA has drift-collapsed most of its within-species
SRK variability, so LEPA-only Shannon entropy scans (Step 22a) call only
1 HV run and miss most of the residues that carry specificity information.
Using the outgroup taxa (Brassica + Arabidopsis) as guiding force recovers
5 cross-genera HV runs (197–430) and 4 intervening linker zones. All 12 Ma
2016 SCR-contact residues fall inside the resulting HV union + linker set.

Under this expanded functional-site definition the per-pair identical
("synonymy") rate in LEPA drops from 42.7 % (LEPA-HV only, current pipeline)
to ≈ 4.8 % (cross-genera HV + polymorphic linkers) — the current pipeline
is massively over-collapsing distinct alleles into synonymy groups.

Inputs (all from Step 22a / Step 22b outputs):
    Tables/Phase5/step26a_HV_regions_per_species.tsv
    Tables/Phase5/step26a_combined_alignment.fasta
    Tables/Phase5/step26a_brassica_hv_mapping.tsv        (Ma 2016 contacts)
    Tables/Phase5/step26a_LEPA_HV_positions.tsv          (for comparison)
    Tables/Phase2/step11_individual_allele_genotypes.tsv (for collapse)

Outputs:
    Tables/Phase5/step26a_functional_site_positions.tsv
        Long-format table: one row per site in the functional set with
        annotation of which species HV it belongs to, whether it is a
        Ma 2016 contact, and whether it is polymorphic within LEPA / in
        the full alignment.

    Tables/Phase5/step26b_functional_synonymy_groups.csv
        One row per LEPA allele, with its "functional group" identifier.
        Groups collapse alleles that are identical across ALL functional
        sites (connected components of the identical-pair graph).

    Tables/Phase5/step26b_functional_allele_distances.tsv
        Long-format pairwise p-distance between LEPA alleles computed on
        the functional site set only. Parallel to
        step26b_HV_allele_distances.tsv (which uses the HV-only set).

    Tables/Phase2/step11_functional_allele_genotypes.tsv
        The Step 11 individual × allele wide count matrix with columns
        collapsed by functional synonymy group. Feeds directly into the
        functional accumulation-curve rerun (see
        SRK_functional_allele_accumulation.R).

The HV-only pipeline outputs are left completely untouched — the new
files are added side-by-side for direct comparison.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import networkx as nx

# ------------------------------------------------------------ paths
BASE = Path("Tables")
IN_HV_RUNS      = BASE / "Phase5" / "step26a_HV_regions_per_species.tsv"
IN_ALN          = BASE / "Phase5" / "step26a_combined_alignment.fasta"
IN_MA_CONTACTS  = BASE / "Phase5" / "step26a_brassica_hv_mapping.tsv"
IN_LEPA_HV      = BASE / "Phase5" / "step26a_LEPA_HV_positions.tsv"
IN_ALLELE_GENO  = BASE / "Phase2" / "step11_individual_allele_genotypes.tsv"

OUT_SITES       = BASE / "Phase5" / "step26a_functional_site_positions.tsv"
OUT_GROUPS      = BASE / "Phase5" / "step26b_functional_synonymy_groups.csv"
OUT_DISTS       = BASE / "Phase5" / "step26b_functional_allele_distances.tsv"
OUT_COLLAPSED   = BASE / "Phase2" / "step11_functional_allele_genotypes.tsv"
OUT_NET_PDF     = Path("figures") / "Phase5" / "step26b_functional_synonymy_network.pdf"
OUT_NET_PNG     = Path("figures") / "Phase5" / "step26b_functional_synonymy_network.png"

DOMAIN_REGION = (31, 430)                    # 1-based cols, per Step 22b
GAP_CHARS = {"-", "X", "?", "B", "Z", "J"}   # skip these when calling polymorphism

# ------------------------------------------------------------ helpers
def load_species_hv_runs(path: Path) -> dict[str, set[int]]:
    df = pd.read_csv(path, sep="\t")
    out = {}
    for species, sub in df.groupby("Species"):
        cols = set()
        for _, r in sub.iterrows():
            cols.update(range(int(r["Start_lepa1"]), int(r["End_lepa1"]) + 1))
        out[species] = cols
    return out

def runs_from_cols(cols: list[int]) -> list[tuple[int, int]]:
    cols = sorted(cols)
    if not cols:
        return []
    runs = []
    start = prev = cols[0]
    for c in cols[1:]:
        if c == prev + 1:
            prev = c
        else:
            runs.append((start, prev)); start = c; prev = c
    runs.append((start, prev))
    return runs

def polymorphic(col_data) -> bool:
    aa = set(col_data) - GAP_CHARS
    return len(aa) >= 2

def hamming_pdist(mat_sub: np.ndarray, allele_ids: list[str]) -> pd.DataFrame:
    """Pairwise p-distance across the given sub-alignment (rows = alleles,
    cols = functional sites). Ignores gap positions in each pair."""
    n = mat_sub.shape[0]
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            mask = ~np.isin(mat_sub[i], list(GAP_CHARS)) & ~np.isin(mat_sub[j], list(GAP_CHARS))
            if mask.sum() == 0:
                d = np.nan
            else:
                d = float((mat_sub[i][mask] != mat_sub[j][mask]).mean())
            rows.append({"Allele_a": allele_ids[i], "Allele_b": allele_ids[j],
                         "p_distance": d,
                         "n_sites_compared": int(mask.sum())})
    return pd.DataFrame(rows)

def draw_functional_network(allele_ids: list[str],
                            id_edges: list[tuple[str, str]],
                            allele_to_group: dict[str, str],
                            group_size: dict[str, int],
                            allele_freq: dict[str, int],
                            n_alleles: int, n_groups: int,
                            out_pdf: Path, out_png: Path) -> None:
    """Functional-site synonymy network — TIDY LAYOUT.

    Multi-allele groups arranged along the top row, one cluster per cell, with
    the group label ABOVE each cluster. Nodes within a cluster arranged in a
    small circle. Singletons placed in a compact grid across the bottom of
    the plot. Nodes coloured by functional group; sized by carrier frequency.
    """
    G = nx.Graph()
    for a in allele_ids:
        G.add_node(a)
    for a, b in id_edges:
        G.add_edge(a, b)

    # Multi-allele groups sorted by size (largest first) for stable colour assignment
    multi_groups = sorted(
        [g for g, sz in group_size.items() if sz > 1],
        key=lambda g: (-group_size[g], g),
    )
    palette = plt.get_cmap("tab10").colors + plt.get_cmap("Set2").colors
    group_colour = {g: palette[i % len(palette)] for i, g in enumerate(multi_groups)}
    SINGLETON_COL = "#B8B8B8"

    # ---- Custom tidy layout ----
    pos: dict[str, tuple[float, float]] = {}
    group_cell_x: dict[str, float] = {}
    group_cell_y: dict[str, float] = {}

    # Multi-allele groups: one row across top
    y_multi = 0.72
    for i, g in enumerate(multi_groups):
        cx = (i + 0.5) / len(multi_groups)
        alleles = sorted([a for a in G.nodes() if allele_to_group[a] == g])
        n = len(alleles)
        r = 0.04 + 0.006 * n                 # radius scales gently with n
        for j, a in enumerate(alleles):
            angle = 2 * np.pi * j / n - np.pi / 2   # start at 12 o'clock
            pos[a] = (cx + r * np.cos(angle), y_multi + r * np.sin(angle))
        group_cell_x[g] = cx
        group_cell_y[g] = y_multi + r + 0.03      # for label placement

    # Singletons: compact grid at bottom
    singletons = sorted(a for a in G.nodes()
                        if group_size[allele_to_group[a]] == 1)
    n_single = len(singletons)
    n_cols = 7
    row_h = 0.08
    grid_top_y = 0.42
    for i, a in enumerate(singletons):
        row = i // n_cols
        col = i % n_cols
        pos[a] = ((col + 0.5) / n_cols, grid_top_y - row * row_h)

    # Node colours + sizes
    node_colors = [group_colour.get(allele_to_group[a], SINGLETON_COL)
                   for a in G.nodes()]
    max_freq = max(allele_freq.values()) if allele_freq else 1
    node_sizes = [220 + 420 * (allele_freq.get(a, 0) / max_freq)
                  for a in G.nodes()]

    # ---- Draw ----
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(14, 8.5))
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 0.98)

    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#555555",
                           width=1.6, alpha=0.75)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, edgecolors="black", linewidths=1)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8, font_weight="bold",
                            labels={a: a.replace("Allele_", "") for a in G.nodes()})

    # Group labels above each multi-allele cluster
    for g in multi_groups:
        ax.text(group_cell_x[g], group_cell_y[g],
                f"{g}  (n = {group_size[g]})",
                ha="center", va="bottom", fontsize=11, fontweight="bold",
                color=group_colour[g])

    # Section label above the singleton grid
    ax.text(0.5, grid_top_y + row_h * 0.65,
            f"Singleton alleles  (n = {n_single})",
            ha="center", va="bottom", fontsize=11, fontweight="bold",
            color="#555555")
    # Faint horizontal divider between the two sections
    ax.plot([0.05, 0.95], [grid_top_y + row_h * 0.55] * 2,
            color="#D0D0D0", linewidth=0.8, zorder=0)

    ax.set_axis_off()
    n_edges = G.number_of_edges()
    fig.suptitle(
        f"Functional synonymy network — {n_alleles} LEPA alleles → {n_groups} functional groups\n"
        f"({len(multi_groups)} multi-allele groups + {n_single} singletons; "
        f"{n_edges} identity edges at functional sites)",
        fontsize=13, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def connected_components(n_nodes: int, edges: list[tuple[int, int]]) -> list[list[int]]:
    parent = list(range(n_nodes))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for a, b in edges:
        union(a, b)
    groups: dict[int, list[int]] = {}
    for i in range(n_nodes):
        r = find(i)
        groups.setdefault(r, []).append(i)
    return list(groups.values())

# ------------------------------------------------------------ main
def main() -> None:
    warnings.filterwarnings("ignore")
    OUT_SITES.parent.mkdir(parents=True, exist_ok=True)
    OUT_COLLAPSED.parent.mkdir(parents=True, exist_ok=True)

    # ---- 1. Cross-genera HV union + linker zones ----
    per_species_hv = load_species_hv_runs(IN_HV_RUNS)
    lepa_hv = per_species_hv.get("LEPA", set())
    brass_hv = per_species_hv.get("Brassica", set())
    arab_hv = per_species_hv.get("Arabidopsis", set())
    hv_union = lepa_hv | brass_hv | arab_hv
    print(f"HV per species — LEPA: {len(lepa_hv)}, "
          f"Brassica: {len(brass_hv)}, Arabidopsis: {len(arab_hv)}")
    print(f"HV union across genera: {len(hv_union)} cols "
          f"(+{len(hv_union - lepa_hv)} beyond LEPA-only)")

    union_runs = runs_from_cols(list(hv_union))
    print(f"HV union runs ({len(union_runs)}):")
    for i, (s, e) in enumerate(union_runs, 1):
        print(f"  Run {i}: cols {s}-{e}  ({e - s + 1} cols)")
    linker_zones = [(union_runs[i][1] + 1, union_runs[i + 1][0] - 1)
                    for i in range(len(union_runs) - 1)]
    print(f"Linker zones between HV runs ({len(linker_zones)}):")
    for i, (s, e) in enumerate(linker_zones, 1):
        print(f"  Linker {i}: cols {s}-{e}  ({e - s + 1} cols)")

    # ---- 2. Load alignment; find polymorphic linker positions ----
    all_seqs = list(SeqIO.parse(IN_ALN, "fasta"))
    lepa_seqs = [s for s in all_seqs if s.id.startswith("Allele_")]
    print(f"Alignment: {len(all_seqs)} total seqs "
          f"({len(lepa_seqs)} LEPA alleles), length {len(all_seqs[0].seq)}")

    mat_all = np.array([list(str(s.seq)) for s in all_seqs])
    mat_lepa = np.array([list(str(s.seq)) for s in lepa_seqs])

    linker_cols = set()
    for s, e in linker_zones:
        linker_cols.update(range(s, e + 1))
    linker_poly_all = {c for c in linker_cols if polymorphic(mat_all[:, c - 1])}
    print(f"Linker cols total: {len(linker_cols)}; "
          f"polymorphic in full alignment: {len(linker_poly_all)}")

    functional_sites = hv_union | linker_poly_all
    print(f"FUNCTIONAL SITE SET: {len(functional_sites)} cols "
          f"(HV union {len(hv_union)} + polymorphic linker {len(linker_poly_all - hv_union)})")

    # ---- 3. Write functional-site table ----
    ma_df = pd.read_csv(IN_MA_CONTACTS, sep="\t")
    ma_cols = set(ma_df["LEPA_aln_col"].astype(int).tolist())

    rows = []
    for c in sorted(functional_sites):
        rows.append({
            "LEPA_aln_col_1based": c,
            "Site_type": "HV_union" if c in hv_union else "linker_polymorphic",
            "In_LEPA_HV":       c in lepa_hv,
            "In_Brassica_HV":   c in brass_hv,
            "In_Arabidopsis_HV": c in arab_hv,
            "In_linker_zone":   c in linker_cols,
            "Ma_2016_contact":  c in ma_cols,
            "Polymorphic_LEPA": polymorphic(mat_lepa[:, c - 1]),
            "Polymorphic_all":  polymorphic(mat_all[:, c - 1]),
        })
    sites_df = pd.DataFrame(rows)
    sites_df.to_csv(OUT_SITES, sep="\t", index=False)
    print(f"→ wrote {OUT_SITES} ({len(sites_df)} rows)")

    # ---- 4. Functional synonymy groups from LEPA alleles ----
    functional_cols = sorted(functional_sites)
    mat_lepa_func = mat_lepa[:, [c - 1 for c in functional_cols]]
    allele_ids = [s.id for s in lepa_seqs]
    n_alleles = len(allele_ids)

    # Compute distance matrix
    dist_df = hamming_pdist(mat_lepa_func, allele_ids)
    dist_df.to_csv(OUT_DISTS, sep="\t", index=False)
    print(f"→ wrote {OUT_DISTS} ({len(dist_df)} pairs)")

    # Build identical-pair edges
    id_to_idx = {a: i for i, a in enumerate(allele_ids)}
    edges = [(id_to_idx[r["Allele_a"]], id_to_idx[r["Allele_b"]])
             for _, r in dist_df.iterrows() if r["p_distance"] == 0]
    components = connected_components(n_alleles, edges)
    components.sort(key=lambda comp: -len(comp))

    group_rows = []
    for gid, comp in enumerate(components, 1):
        group_id = f"FG{gid:03d}"
        for idx in sorted(comp):
            group_rows.append({
                "Allele": allele_ids[idx],
                "Functional_group": group_id,
                "Group_size": len(comp),
                "Group_type": "synonymy" if len(comp) > 1 else "singleton",
            })
    groups_df = pd.DataFrame(group_rows).sort_values(["Functional_group", "Allele"])
    groups_df.to_csv(OUT_GROUPS, index=False)

    n_groups = len(components)
    n_multi = sum(1 for c in components if len(c) > 1)
    n_single = n_groups - n_multi
    n_collapsed = sum(len(c) for c in components if len(c) > 1)
    print(f"Functional synonymy: {n_alleles} alleles → {n_groups} groups "
          f"({n_multi} multi-allele, {n_single} singletons; "
          f"{n_collapsed} alleles collapsed into multi-allele groups)")
    print(f"→ wrote {OUT_GROUPS}")

    # ---- 5. Collapse the individual × allele matrix by functional group ----
    geno = pd.read_csv(IN_ALLELE_GENO, sep="\t", encoding="utf-8-sig")
    allele_to_group = dict(zip(groups_df["Allele"], groups_df["Functional_group"]))

    # Map every allele column to its functional group; sum within groups
    allele_cols_in_matrix = [c for c in geno.columns if c != "Individual"]
    missing = [a for a in allele_cols_in_matrix if a not in allele_to_group]
    if missing:
        print(f"WARNING: {len(missing)} alleles in step11 lack a functional group "
              f"(first few: {missing[:5]}). They will be kept as their own singletons.")
        # Assign each missing to a unique group
        next_gid = len(components) + 1
        for a in missing:
            allele_to_group[a] = f"FG{next_gid:03d}"
            next_gid += 1

    # Build collapsed matrix
    long = geno.melt(id_vars="Individual", var_name="Allele", value_name="Count")
    long["Functional_group"] = long["Allele"].map(allele_to_group)
    collapsed = (long.groupby(["Individual", "Functional_group"], sort=False)["Count"]
                     .sum()
                     .unstack("Functional_group", fill_value=0)
                     .reset_index())
    # Reorder columns: Individual first, then functional groups sorted alphabetically
    other_cols = sorted([c for c in collapsed.columns if c != "Individual"])
    collapsed = collapsed[["Individual"] + other_cols]
    collapsed.to_csv(OUT_COLLAPSED, sep="\t", index=False)
    print(f"→ wrote {OUT_COLLAPSED} "
          f"({collapsed.shape[0]} individuals × {collapsed.shape[1] - 1} functional groups)")

    # ---- 6. Functional synonymy network figure ----
    # Allele carrier frequency (n individuals with count > 0) — for node sizing
    allele_freq = {a: int((geno[a] > 0).sum()) for a in allele_cols_in_matrix
                   if a in geno.columns}
    id_edges = [(r["Allele_a"], r["Allele_b"]) for _, r in dist_df.iterrows()
                if r["p_distance"] == 0]
    group_size = {g: int(sz) for g, sz in
                  groups_df.groupby("Functional_group")["Allele"].nunique().items()}
    draw_functional_network(
        allele_ids=allele_ids,
        id_edges=id_edges,
        allele_to_group=allele_to_group,
        group_size=group_size,
        allele_freq=allele_freq,
        n_alleles=n_alleles,
        n_groups=n_groups,
        out_pdf=OUT_NET_PDF,
        out_png=OUT_NET_PNG,
    )
    print(f"→ wrote {OUT_NET_PDF}")
    print(f"→ wrote {OUT_NET_PNG}")


if __name__ == "__main__":
    main()
