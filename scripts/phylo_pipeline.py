#!/usr/bin/env python3
"""
phylo_pipeline.py — Pure-Python Molecular Evolutionary Genetics Pipeline
========================================================================
Implements JC69 and Kimura 2-Parameter (K2P) evolutionary distance models
from their first-principles mathematical definitions, followed by Saitou-Nei (1987)
Neighbor-Joining (NJ) tree reconstruction and Felsenstein (1985) bootstrap resampling
— all executed without proprietary phylogenetics GUI software.

Supported Datasets:
  1. col1a1   : 8 vertebrate comparative species (Homo sapiens, Pan troglodytes,
                Mus musculus, Rattus norvegicus, Bos taurus, Gallus gallus,
                Danio rerio, Xenopus tropicalis).
  2. animals  : 9 animal & outgroup taxa alignment from MEGA format (.meg).
  3. tro_seq  : 30 zoonotic & sandfly vector epidemiological surveillance isolates.

Dependencies: biopython, numpy, matplotlib
Usage:
    python scripts/phylo_pipeline.py --dataset col1a1 --model k2p --bootstrap 100
    python scripts/phylo_pipeline.py --dataset animals --model k2p --bootstrap 100
    python scripts/phylo_pipeline.py --dataset tro_seq --model k2p --bootstrap 100

Outputs:
    results/{dataset}_analysis/{dataset}_nj_tree_{model}.nwk
    results/{dataset}_analysis/{dataset}_nj_tree_{model}.png

Author : Suleyman Hajizadeh
Organization: SuleimanMEGAMolecularEvolutionaryTree
"""

import os
import math
import time
import logging
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger("phylo_pipeline")

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION & PATHS
# ─────────────────────────────────────────────────────────────────────────────

DEFAULT_EMAIL = "suleyman.hacizade1@gmail.com"

# COL1A1 RefSeq accessions across 8 vertebrate species
ACCESSIONS = {
    "Homo_sapiens":        "NM_000088",
    "Pan_troglodytes":     "XM_016944896",
    "Mus_musculus":        "NM_007742",
    "Rattus_norvegicus":   "NM_053304",
    "Bos_taurus":          "NM_174520",
    "Gallus_gallus":       "NM_204775",
    "Danio_rerio":         "NM_001004631",
    "Xenopus_tropicalis":  "NM_001016853",
}

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data", "raw")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
COL1A1_LOCAL_DIR = os.path.join(DATA_DIR, "col1a1_sequences")
ANIMALS_MEG_PATH = os.path.join(DATA_DIR, "animal_alignments", "animals.meg")
TRO_SEQ_PATH = os.path.join(DATA_DIR, "col1a1_sequences", "TRO_Seq(30).fas")


# ─────────────────────────────────────────────────────────────────────────────
# 1. DATA INGESTION & PARSING
# ─────────────────────────────────────────────────────────────────────────────

def parse_meg_file(filepath: str) -> dict:
    """
    Parse a MEGA format (.meg) nucleotide alignment file.

    Parameters
    ----------
    filepath : str
        Path to the .meg file.

    Returns
    -------
    dict
        Mapping of taxon_name -> aligned nucleotide sequence string.
        Taxa with 100% gap/missing characters are omitted.
    """
    sequences = {}
    current_taxon = None
    current_seq_parts = []

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                header = line[1:].strip()
                if header.upper() == "MEGA":
                    continue
                if current_taxon:
                    seq_str = "".join(current_seq_parts).replace(" ", "").upper()
                    if any(c in "ACGTURYSWKMBDHV" for c in seq_str):
                        sequences[current_taxon] = seq_str
                    current_seq_parts = []
                current_taxon = header
            elif line.startswith("!"):
                continue
            else:
                if current_taxon:
                    current_seq_parts.append(line.replace(" ", ""))

        if current_taxon:
            seq_str = "".join(current_seq_parts).replace(" ", "").upper()
            if any(c in "ACGTURYSWKMBDHV" for c in seq_str):
                sequences[current_taxon] = seq_str

    return sequences


def parse_fasta_file(filepath: str) -> dict:
    """
    Parse a FASTA file into a mapping of identifier -> sequence string.
    """
    records = {}
    for rec in SeqIO.parse(filepath, "fasta"):
        records[rec.id] = str(rec.seq).upper()
    return records


def fetch_sequences(accessions: dict, local_dir: str = None, rettype: str = "fasta") -> dict:
    """
    Retrieve nucleotide sequences from local cache or NCBI Entrez.

    Parameters
    ----------
    accessions : dict
        Mapping of species_name -> RefSeq accession string.
    local_dir : str, optional
        Directory containing local FASTA files for offline/cached execution.
    rettype : str
        Entrez return type ('fasta' or 'gb').

    Returns
    -------
    dict
        Mapping of species_name -> sequence string.
    """
    if local_dir is None:
        local_dir = COL1A1_LOCAL_DIR

    records = {}
    for species, acc in accessions.items():
        found_local = False
        if local_dir and os.path.isdir(local_dir):
            for fname in os.listdir(local_dir):
                if species.lower() in fname.lower() and (fname.endswith(".fasta") or fname.endswith(".fas") or fname.endswith(".fa")):
                    fpath = os.path.join(local_dir, fname)
                    try:
                        record = SeqIO.read(fpath, "fasta")
                        records[species] = str(record.seq).upper()
                        logger.info(f"Loaded {species} from local cache: {fname}")
                        found_local = True
                        break
                    except Exception as parse_err:
                        logger.debug(f"Could not parse local file {fname}: {parse_err}")

        if found_local:
            continue

        logger.info(f"Fetching {species} ({acc}) from NCBI Entrez …")
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype=rettype, retmode="text")
            record = SeqIO.read(handle, rettype)
            handle.close()
            records[species] = str(record.seq).upper()
            time.sleep(0.4)
        except Exception as exc:
            logger.warning(f"Could not fetch {acc}: {exc}")
    return records


def load_dataset(name: str) -> dict:
    """
    Load sequence dataset based on keyword identifier.
    """
    if name == "col1a1":
        return fetch_sequences(ACCESSIONS)
    elif name == "animals":
        if not os.path.exists(ANIMALS_MEG_PATH):
            raise FileNotFoundError(f"Animals alignment not found at: {ANIMALS_MEG_PATH}")
        return parse_meg_file(ANIMALS_MEG_PATH)
    elif name == "tro_seq":
        if not os.path.exists(TRO_SEQ_PATH):
            raise FileNotFoundError(f"TRO_Seq fasta not found at: {TRO_SEQ_PATH}")
        return parse_fasta_file(TRO_SEQ_PATH)
    else:
        raise ValueError(f"Unknown dataset '{name}'. Choose 'col1a1', 'animals', or 'tro_seq'.")


# ─────────────────────────────────────────────────────────────────────────────
# 2. ALIGNMENT & COUNTS
# ─────────────────────────────────────────────────────────────────────────────

def align_pair(seq_a: str, seq_b: str) -> tuple:
    """
    Global pairwise alignment of two nucleotide sequences using Biopython.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(seq_a, seq_b)
    best = next(iter(alignments))
    aligned_a = str(best).split("\n")[0]
    aligned_b = str(best).split("\n")[2]
    return aligned_a, aligned_b


def compute_alignment_counts(aligned_a: str, aligned_b: str) -> dict:
    """
    Classify aligned nucleotide columns into matches, transitions (P), transversions (Q), and gaps.
    """
    purines = set("AG")
    pyrimidines = set("CT")
    n_match = n_P = n_Q = n_gap = 0

    for a, b in zip(aligned_a, aligned_b):
        if a not in "ACGT" or b not in "ACGT":
            n_gap += 1
            continue
        if a == b:
            n_match += 1
        elif (a in purines and b in purines) or (a in pyrimidines and b in pyrimidines):
            n_P += 1
        else:
            n_Q += 1

    n_sites = n_match + n_P + n_Q
    p = n_P / n_sites if n_sites > 0 else 0.0
    q = n_Q / n_sites if n_sites > 0 else 0.0
    return dict(n_sites=n_sites, n_P=n_P, n_Q=n_Q, p=p, q=q)


# ─────────────────────────────────────────────────────────────────────────────
# 3. MATHEMATICAL EVOLUTIONARY DISTANCE MODELS
# ─────────────────────────────────────────────────────────────────────────────

def jc69_distance(p_total: float) -> float:
    """
    Jukes-Cantor 1969 (JC69) evolutionary distance:
        d = -3/4 * ln(1 - 4/3 * p)
    """
    arg = 1.0 - (4.0 / 3.0) * p_total
    if arg <= 0:
        return np.inf
    return -0.75 * math.log(arg)


def k2p_distance(p: float, q: float) -> float:
    """
    Kimura 2-Parameter (K2P, 1980) evolutionary distance:
        d = -1/2 * ln(1 - 2P - Q) - 1/4 * ln(1 - 2Q)
    """
    w1 = 1.0 - 2.0 * p - q
    w2 = 1.0 - 2.0 * q
    if w1 <= 0 or w2 <= 0:
        return np.inf
    return -0.5 * math.log(w1) - 0.25 * math.log(w2)


def build_distance_matrix(sequences: dict, model: str = "k2p") -> tuple:
    """
    Compute pairwise evolutionary distance matrix and cache aligned sequences.

    Returns
    -------
    (labels, matrix, aligned_pairs)
    """
    labels = list(sequences.keys())
    n = len(labels)
    mat = np.zeros((n, n))
    aligned_pairs = {}

    for i in range(n):
        for j in range(i + 1, n):
            sa = sequences[labels[i]]
            sb = sequences[labels[j]]

            # If pre-aligned (equal length and containing alignment characters)
            if len(sa) == len(sb) and ("-" in sa or "-" in sb):
                aa, ab = sa, sb
            else:
                min_len = min(len(sa), len(sb))
                aa, ab = align_pair(sa[:min_len], sb[:min_len])

            aligned_pairs[(i, j)] = (aa, ab)
            counts = compute_alignment_counts(aa, ab)
            if model == "jc69":
                d = jc69_distance(counts["p"] + counts["q"])
            else:
                d = k2p_distance(counts["p"], counts["q"])

            if not np.isfinite(d):
                d = 2.5  # Numerical ceiling for deeply saturated/divergent comparisons
            mat[i, j] = mat[j, i] = d

    return labels, mat, aligned_pairs


# ─────────────────────────────────────────────────────────────────────────────
# 4. NEIGHBOR-JOINING ALGORITHM (Saitou & Nei, 1987)
# ─────────────────────────────────────────────────────────────────────────────

def neighbor_joining(labels: list, dist_matrix: np.ndarray) -> list:
    """
    Reconstruct a phylogenetic tree using Saitou & Nei (1987) Neighbor-Joining.
    """
    labels = list(labels)
    D = dist_matrix.copy().astype(float)
    joins = []

    while len(labels) > 2:
        n = len(labels)
        row_sums = D.sum(axis=1)

        # Compute Q-matrix
        Q = np.full((n, n), np.inf)
        for i in range(n):
            for j in range(i + 1, n):
                Q[i, j] = (n - 2) * D[i, j] - row_sums[i] - row_sums[j]
                Q[j, i] = Q[i, j]

        np.fill_diagonal(Q, np.inf)
        idx = np.unravel_index(np.argmin(Q), Q.shape)
        i, j = int(idx[0]), int(idx[1])

        d_ij = D[i, j]
        d_iu = d_ij / 2.0 + (row_sums[i] - row_sums[j]) / (2.0 * (n - 2)) if n > 2 else d_ij / 2.0
        d_ju = d_ij - d_iu

        joins.append((labels[i], labels[j], max(d_iu, 0), max(d_ju, 0)))

        new_label = f"({labels[i]},{labels[j]})"
        new_dists = []
        for k in range(n):
            if k != i and k != j:
                d_ku = (D[i, k] + D[j, k] - d_ij) / 2.0
                new_dists.append(d_ku)

        remaining = [k for k in range(n) if k != i and k != j]
        m = len(remaining)
        new_D = np.zeros((m + 1, m + 1))
        for a_idx, ri in enumerate(remaining):
            for b_idx, rj in enumerate(remaining):
                new_D[a_idx, b_idx] = D[ri, rj]
            new_D[a_idx, m] = new_D[m, a_idx] = new_dists[a_idx]

        D = new_D
        labels = [labels[k] for k in remaining] + [new_label]

    if len(labels) == 2:
        joins.append((labels[0], labels[1], max(D[0, 1] / 2.0, 0), max(D[0, 1] / 2.0, 0)))

    return joins


# ─────────────────────────────────────────────────────────────────────────────
# 5. BOOTSTRAP RESAMPLING & SPLIT SUPPORT (Felsenstein, 1985)
# ─────────────────────────────────────────────────────────────────────────────

def canonical_split(clade, all_taxa: list):
    """
    Normalize an unrooted phylogenetic bipartition (split) against a reference taxon.
    """
    ref_taxon = all_taxa[0]
    total_taxa = set(all_taxa)
    clade_set = set(clade)
    if ref_taxon in clade_set:
        return frozenset(total_taxa - clade_set)
    return frozenset(clade_set)


def get_splits_from_joins(joins: list, all_taxa: list) -> set:
    """
    Extract all internal bipartitions (splits) formed by NJ join steps.
    """
    leaf_map = {t: {t} for t in all_taxa}
    splits = set()
    for a, b, _, _ in joins[:-1]:
        leaves_a = leaf_map.get(a, {a})
        leaves_b = leaf_map.get(b, {b})
        merged = leaves_a | leaves_b
        parent_key = f"({a},{b})"
        leaf_map[parent_key] = merged
        if 2 <= len(merged) <= len(all_taxa) - 2:
            splits.add(canonical_split(merged, all_taxa))
    return splits


def bootstrap_neighbor_joining(labels: list, aligned_pairs: dict,
                               model: str = "k2p", n_boot: int = 100) -> dict:
    """
    Perform Felsenstein non-parametric bootstrap resampling and compute clade support percentages.
    """
    n = len(labels)
    split_counts = {}

    for _ in range(n_boot):
        D_b = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                sa, sb = aligned_pairs[(i, j)]
                L = len(sa)
                if L == 0:
                    d = 0.0
                else:
                    cols = np.random.choice(L, size=L, replace=True)
                    res_a = "".join(sa[k] for k in cols)
                    res_b = "".join(sb[k] for k in cols)
                    counts = compute_alignment_counts(res_a, res_b)
                    if model == "jc69":
                        d = jc69_distance(counts["p"] + counts["q"])
                    else:
                        d = k2p_distance(counts["p"], counts["q"])
                    if not np.isfinite(d):
                        d = 2.5
                D_b[i, j] = D_b[j, i] = d

        joins_b = neighbor_joining(labels, D_b)
        splits_b = get_splits_from_joins(joins_b, labels)
        for sp in splits_b:
            split_counts[sp] = split_counts.get(sp, 0) + 1

    return {sp: (count / n_boot) * 100.0 for sp, count in split_counts.items()}


def format_newick(joins: list, all_taxa: list = None, split_support: dict = None) -> str:
    """
    Format NJ joins into standard Newick tree notation with branch lengths and internal bootstrap support.
    """
    if not joins:
        return ";"
    nodes = {}
    leaf_map = {t: {t} for t in all_taxa} if all_taxa else {}

    for idx, (a, b, da, db) in enumerate(joins):
        child_a = f"{nodes[a]}:{da:.5f}" if a in nodes else f"{a}:{da:.5f}"
        child_b = f"{nodes[b]}:{db:.5f}" if b in nodes else f"{b}:{db:.5f}"
        parent_key = f"({a},{b})"

        support_str = ""
        if all_taxa and split_support and idx < len(joins) - 1:
            leaves_a = leaf_map.get(a, {a})
            leaves_b = leaf_map.get(b, {b})
            merged = leaves_a | leaves_b
            leaf_map[parent_key] = merged
            if 2 <= len(merged) <= len(all_taxa) - 2:
                can_sp = canonical_split(merged, all_taxa)
                val = split_support.get(can_sp, None)
                if val is not None:
                    support_str = f"{int(round(val))}"
        elif all_taxa:
            leaves_a = leaf_map.get(a, {a})
            leaves_b = leaf_map.get(b, {b})
            leaf_map[parent_key] = leaves_a | leaves_b

        nodes[parent_key] = f"({child_a},{child_b}){support_str}"

    last_key = f"({joins[-1][0]},{joins[-1][1]})"
    return nodes.get(last_key, "") + ";"


# ─────────────────────────────────────────────────────────────────────────────
# 6. PUBLICATION-QUALITY TREE VISUALISATION
# ─────────────────────────────────────────────────────────────────────────────

def plot_nj_tree(joins: list, labels_orig: list, output_path: str,
                 dist_matrix: np.ndarray, model: str = "k2p",
                 dataset_title: str = "COL1A1", split_support: dict = None) -> None:
    """
    Plot the Neighbor-Joining tree cladogram and distance matrix heatmap in publication dark mode.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, max(7, len(labels_orig) * 0.45)),
                             gridspec_kw={"width_ratios": [2.2, 1.2]})
    fig.patch.set_facecolor("#070b14")
    for ax in axes:
        ax.set_facecolor("#0f172a")

    ax = axes[0]
    y_positions = {sp: i for i, sp in enumerate(labels_orig)}
    leaf_map = {sp: {sp} for sp in labels_orig}
    colors = plt.cm.tab10(np.linspace(0, 1, len(labels_orig)))
    color_map = {sp: c for sp, c in zip(labels_orig, colors)}

    current_y = len(labels_orig)
    for idx, join in enumerate(joins):
        a, b, d_a, d_b = join
        ya = y_positions.get(a, current_y)
        yb = y_positions.get(b, current_y + 1)
        mid_y = (ya + yb) / 2.0

        x_a_start = 0 if a in labels_orig else 0.02
        ax.plot([x_a_start, d_a], [ya, ya],
                color=color_map.get(a, "#38bdf8"), linewidth=2.4, solid_capstyle="round")
        ax.plot([x_a_start, d_b], [yb, yb],
                color=color_map.get(b, "#38bdf8"), linewidth=2.4, solid_capstyle="round")

        max_d = max(d_a, d_b)
        ax.plot([max_d, max_d], [ya, yb],
                color="#00f2fe", linewidth=1.5, linestyle="--", alpha=0.8)

        # Annotate bootstrap support if available
        new_label = f"({a},{b})"
        leaves_a = leaf_map.get(a, {a})
        leaves_b = leaf_map.get(b, {b})
        merged = leaves_a | leaves_b
        leaf_map[new_label] = merged

        if split_support and idx < len(joins) - 1:
            if 2 <= len(merged) <= len(labels_orig) - 2:
                can_sp = canonical_split(merged, labels_orig)
                val = split_support.get(can_sp, None)
                if val is not None:
                    ax.text(max_d + 0.002, mid_y, f"{int(round(val))}%",
                            color="#38bdf8", fontsize=8, fontweight="bold",
                            va="center", ha="left")

        y_positions[new_label] = mid_y
        color_map[new_label] = "#00f2fe"
        current_y += 1

    for sp, y in [(s, y_positions[s]) for s in labels_orig]:
        ax.text(0.001, y, sp.replace("_", " "), color=color_map.get(sp, "white"),
                fontsize=8.5, va="center", fontstyle="italic")

    ax.set_title(f"{dataset_title} Neighbor-Joining Phylogeny\n({model.upper()} Distance Model • Felsenstein Bootstrap)",
                 color="white", fontsize=13, pad=12, fontweight="bold")
    ax.set_xlabel(f"{model.upper()} evolutionary distance (substitutions/site)", color="#94a3b8")
    ax.tick_params(colors="#94a3b8")
    for spine in ax.spines.values():
        spine.set_edgecolor((0.22, 0.74, 0.97, 0.3))

    # Heatmap
    ax2 = axes[1]
    n = len(labels_orig)
    clean_labels = [s.replace("_", " ") for s in labels_orig]
    finite_vals = dist_matrix[np.isfinite(dist_matrix)]
    vmax = np.nanmax(finite_vals) if len(finite_vals) > 0 else 1.0

    im = ax2.imshow(dist_matrix, cmap="magma_r", aspect="auto", vmin=0, vmax=vmax)
    ax2.set_xticks(range(n))
    ax2.set_yticks(range(n))
    ax2.set_xticklabels(clean_labels, rotation=45, ha="right",
                        color="#94a3b8", fontsize=max(6, 9 - n // 4), fontstyle="italic")
    ax2.set_yticklabels(clean_labels, color="#94a3b8",
                        fontsize=max(6, 9 - n // 4), fontstyle="italic")

    if n <= 12:
        for i in range(n):
            for j in range(n):
                val = dist_matrix[i, j]
                ax2.text(j, i, f"{val:.3f}" if np.isfinite(val) else "∞",
                         ha="center", va="center", fontsize=7,
                         color="#070b14" if val > vmax * 0.55 else "white")

    cbar = fig.colorbar(im, ax=ax2, pad=0.04)
    cbar.set_label(f"{model.upper()} Pairwise Distance", color="#94a3b8")
    cbar.ax.yaxis.set_tick_params(color="#94a3b8")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="#94a3b8")
    ax2.set_title(f"Pairwise {model.upper()} Distance Matrix", color="white", fontsize=11, pad=10)
    for spine in ax2.spines.values():
        spine.set_edgecolor((0.22, 0.74, 0.97, 0.3))

    fig.suptitle("Suleiman MEGA — Molecular Evolutionary Genetics Analysis",
                 color="#f1f5f9", fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    logger.info(f"Phylogeny plot saved to: {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 7. MAIN CONTROLLER
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Suleiman MEGA — Pure-Python Molecular Evolutionary Genetics Pipeline"
    )
    parser.add_argument(
        "--dataset",
        choices=["col1a1", "animals", "tro_seq"],
        default="col1a1",
        help="Dataset to analyze (col1a1, animals, or tro_seq; default: col1a1)"
    )
    parser.add_argument(
        "--model",
        choices=["k2p", "jc69"],
        default="k2p",
        help="Evolutionary distance model (default: k2p)"
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=100,
        help="Number of bootstrap replicates (default: 100, 0 to disable)"
    )
    parser.add_argument(
        "--email",
        default=DEFAULT_EMAIL,
        help=f"NCBI Entrez email address (default: {DEFAULT_EMAIL})"
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: results/{dataset}_analysis)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    Entrez.email = args.email

    dataset_name = args.dataset
    if args.output_dir is None:
        output_dir = os.path.join(RESULTS_DIR, f"{dataset_name}_analysis")
    else:
        output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    logger.info("=" * 65)
    logger.info(f"Suleiman MEGA Pipeline: Dataset={dataset_name.upper()} | Model={args.model.upper()}")
    logger.info("=" * 65)

    # 1. Ingest Dataset
    logger.info(f"Step [1/4]: Ingesting {dataset_name} sequences …")
    seqs = load_dataset(dataset_name)
    logger.info(f"Loaded {len(seqs)} taxa successfully.")

    # 2. Distance Matrix Calculation
    logger.info(f"Step [2/4]: Computing pairwise {args.model.upper()} distances …")
    labels, dist_mat, aligned_pairs = build_distance_matrix(seqs, model=args.model)

    # 3. Neighbor-Joining Tree
    logger.info("Step [3/4]: Reconstructing Neighbor-Joining tree …")
    joins = neighbor_joining(labels, dist_mat)

    # 4. Bootstrap Resampling
    split_support = None
    if args.bootstrap > 0:
        logger.info(f"Executing {args.bootstrap} Felsenstein bootstrap replicates …")
        split_support = bootstrap_neighbor_joining(labels, aligned_pairs, model=args.model, n_boot=args.bootstrap)
        for sp, pct in split_support.items():
            logger.info(f"  Clade support: {list(sp)[:3]}... = {pct:.1f}%")

    # Format Newick Tree
    newick_tree = format_newick(joins, all_taxa=labels, split_support=split_support)
    logger.info(f"Reconstructed Newick Tree:\n{newick_tree}")

    # Export Newick File
    nwk_file = os.path.join(output_dir, f"{dataset_name}_nj_tree_{args.model}.nwk")
    with open(nwk_file, "w") as f:
        f.write(newick_tree + "\n")
    logger.info(f"Exported Newick tree to: {nwk_file}")

    # Export Figure
    png_file = os.path.join(output_dir, f"{dataset_name}_nj_tree_{args.model}.png")
    plot_nj_tree(
        joins=joins,
        labels_orig=labels,
        output_path=png_file,
        dist_matrix=dist_mat,
        model=args.model,
        dataset_title=dataset_name.upper(),
        split_support=split_support
    )

    logger.info(f"✅ Pipeline completed for {dataset_name} ({args.model.upper()}).")


if __name__ == "__main__":
    main()
