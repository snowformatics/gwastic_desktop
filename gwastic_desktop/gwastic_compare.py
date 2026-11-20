#!/usr/bin/env python3
"""
gwastic_compare.py  –  Compare GWAS results across methods and phenotypes
=========================================================================
Reads all gwas_results.csv files produced by gwastic_batch.py and builds
a cross-method / cross-phenotype comparison showing which significant SNPs
(peaks) appear in which algorithm.

Usage
-----
  python gwastic_compare.py --results-dir Z:\free_wheat\From_scratch\gbs [OPTIONS]

Output files (written to --results-dir or --output)
----------------------------------------------------
  comparison_all_significant_snps.csv   – every sig. SNP with one column per method
  comparison_summary.csv                – per-phenotype hit-count per method
  comparison_consensus_snps.csv         – SNPs significant in >= --min-methods methods
  comparison_method_overlap.csv         – pairwise method overlap counts
  plots/venn_<phenotype>.png            – Venn / UpSet overlap diagram per phenotype
  plots/heatmap_<phenotype>.png         – SNP × method significance heatmap
  plots/manhattan_comparison_<pheno>.png– stacked Manhattan per phenotype
"""

import argparse
import glob
import logging
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("gwastic_compare")

# Expected column names from gwastic_desktop output
# We try several common variants for robustness
SNP_COL_CANDIDATES   = ["SNP", "snp", "rs", "MarkerName", "ID", "sid"]
CHR_COL_CANDIDATES   = ["Chr", "CHR", "chr", "chromosome", "Chromosome"]
POS_COL_CANDIDATES   = ["Pos", "POS", "pos", "BP", "bp", "position", "Position"]
PVAL_COL_CANDIDATES  = ["PValue", "pvalue", "p_value", "P", "p", "Pvalue", "pval"]
SCORE_COL_CANDIDATES = ["Score", "score", "effect", "Effect", "importance", "Importance",
                         "coefficient", "Coefficient", "NegLogP", "neg_log_p"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _pick_col(df: pd.DataFrame, candidates: list, label: str):
    """Return the first matching column name or None."""
    for c in candidates:
        if c in df.columns:
            return c
    log.warning("  Could not find %s column. Tried: %s. Available: %s",
                label, candidates, list(df.columns))
    return None


def _safe_read_csv(path: str):
    try:
        df = pd.read_csv(path, sep=None, engine="python")
        return df
    except Exception as e:
        log.error("  Failed to read %s: %s", path, e)
        return None


def _parse_run_dir_name(dirname: str):
    """
    Expected format: <timestamp>_<phenotype>_<Algorithm_Name>
    e.g.  20261119_143201_Browning_WGS_SAMEA_GWAS_FaST-LMM
    Returns (phenotype, algorithm) strings.
    """
    parts = dirname.split("_")
    # first two parts are date + time
    if len(parts) < 4:
        return dirname, "Unknown"
    # algorithm is the last token(s) after the phenotype name
    # We recognise known algorithm tags
    known_algs = {
        "FaST-LMM": "FaST-LMM",
        "Linear": "Linear regression",
        "Ridge": "Ridge Regression",
        "Random": "Random Forest (AI)",
        "XGBoost": "XGBoost (AI)",
    }
    alg_name = "Unknown"
    cut = len(parts)
    for i in range(len(parts) - 1, 1, -1):
        for key, full in known_algs.items():
            if parts[i].startswith(key):
                alg_name = full
                cut = i
                break
        if alg_name != "Unknown":
            break
    pheno = "_".join(parts[2:cut])
    return pheno, alg_name


def _significance_threshold(df: pd.DataFrame, pval_col, score_col,
                              pval_threshold: float, top_n: int):
    """
    Return a boolean mask for significant rows.
    Uses p-value threshold if a p-value column exists,
    otherwise falls back to top-N by score.
    """
    if pval_col and pval_col in df.columns:
        try:
            pvals = pd.to_numeric(df[pval_col], errors="coerce")
            return pvals <= pval_threshold
        except Exception:
            pass
    if score_col and score_col in df.columns:
        log.info("    No p-value column – using top-%d by score.", top_n)
        scores = pd.to_numeric(df[score_col], errors="coerce").fillna(0)
        threshold = scores.nlargest(top_n).min()
        return scores >= threshold
    return pd.Series([False] * len(df))


# ---------------------------------------------------------------------------
# Main discovery + loading
# ---------------------------------------------------------------------------

def discover_runs(results_dir: str):
    """
    Walk results_dir and collect all gwas_results.csv files.
    Returns list of dicts with keys: path, phenotype, algorithm, run_dir.
    """
    runs = []
    for run_dir in sorted(glob.glob(os.path.join(results_dir, "*"))):
        if not os.path.isdir(run_dir):
            continue
        csv_path = os.path.join(run_dir, "gwas_results.csv")
        if not os.path.isfile(csv_path):
            continue
        dirname = os.path.basename(run_dir)
        pheno, alg = _parse_run_dir_name(dirname)

        # Also try to read phenotype / algorithm from log.txt for accuracy
        log_path = os.path.join(run_dir, "log.txt")
        if os.path.isfile(log_path):
            with open(log_path) as f:
                for line in f:
                    if line.startswith("Algorithm"):
                        alg = line.split(":", 1)[1].strip()
                    if line.startswith("Phenotype name"):
                        pheno = line.split(":", 1)[1].strip()

        runs.append({"path": csv_path, "phenotype": pheno,
                     "algorithm": alg, "run_dir": run_dir})
    return runs


def load_all_runs(runs: list, pval_threshold: float, top_n: int):
    """
    For each run, load CSV, detect columns, flag significant SNPs.
    Returns dict: {phenotype: {algorithm: DataFrame_of_sig_snps}}
    """
    data = defaultdict(dict)

    for run in runs:
        pheno = run["phenotype"]
        alg   = run["algorithm"]
        log.info("  Loading  %-40s  [%s]", pheno, alg)

        df = _safe_read_csv(run["path"])
        if df is None or df.empty:
            continue

        snp_col   = _pick_col(df, SNP_COL_CANDIDATES,   "SNP")
        chr_col   = _pick_col(df, CHR_COL_CANDIDATES,   "Chr")
        pos_col   = _pick_col(df, POS_COL_CANDIDATES,   "Pos")
        pval_col  = _pick_col(df, PVAL_COL_CANDIDATES,  "PValue")
        score_col = _pick_col(df, SCORE_COL_CANDIDATES, "Score")

        if snp_col is None:
            log.warning("    No SNP column found – skipping run.")
            continue

        mask = _significance_threshold(df, pval_col, score_col,
                                        pval_threshold, top_n)
        sig = df[mask].copy()

        # Standardise column names
        rename = {}
        if snp_col:   rename[snp_col]   = "SNP"
        if chr_col:   rename[chr_col]   = "Chr"
        if pos_col:   rename[pos_col]   = "Pos"
        if pval_col:  rename[pval_col]  = "PValue"
        if score_col: rename[score_col] = "Score"
        sig = sig.rename(columns=rename)

        # Keep only useful columns
        keep = [c for c in ["SNP", "Chr", "Pos", "PValue", "Score"] if c in sig.columns]
        sig = sig[keep].drop_duplicates(subset=["SNP"])

        log.info("    Significant SNPs: %d", len(sig))
        data[pheno][alg] = sig

    return data


# ---------------------------------------------------------------------------
# Analysis functions
# ---------------------------------------------------------------------------

def build_comparison_table(data: dict):
    """
    For each phenotype: build a wide table SNP × algorithm (1/0 presence).
    Returns dict {phenotype: DataFrame}
    """
    tables = {}
    for pheno, alg_dict in data.items():
        all_snps = pd.concat(alg_dict.values(), ignore_index=True)
        if all_snps.empty:
            continue

        # Base SNP info
        base_cols = [c for c in ["SNP", "Chr", "Pos"] if c in all_snps.columns]
        base = all_snps[base_cols].drop_duplicates(subset=["SNP"]).set_index("SNP")

        for alg, df in alg_dict.items():
            base[alg] = base.index.isin(df["SNP"].values).astype(int)

        base = base.reset_index()
        base["methods_count"] = base[list(alg_dict.keys())].sum(axis=1)
        base = base.sort_values("methods_count", ascending=False)
        tables[pheno] = base

    return tables


def build_summary(data: dict):
    """Per-phenotype, per-algorithm significant SNP counts."""
    rows = []
    for pheno, alg_dict in data.items():
        row = {"Phenotype": pheno}
        for alg, df in alg_dict.items():
            row[alg] = len(df)
        rows.append(row)
    return pd.DataFrame(rows)


def build_consensus(tables: dict, min_methods: int):
    """SNPs found in >= min_methods methods."""
    frames = []
    for pheno, df in tables.items():
        sub = df[df["methods_count"] >= min_methods].copy()
        sub.insert(0, "Phenotype", pheno)
        frames.append(sub)
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame()


def build_pairwise_overlap(data: dict):
    """Pairwise overlap count between algorithms (summed across phenotypes)."""
    all_algs = sorted({alg for alg_dict in data.values() for alg in alg_dict})
    matrix = pd.DataFrame(0, index=all_algs, columns=all_algs)

    for pheno, alg_dict in data.items():
        algs = list(alg_dict.keys())
        for i, a1 in enumerate(algs):
            for a2 in algs[i:]:
                snps1 = set(alg_dict[a1]["SNP"].values)
                snps2 = set(alg_dict[a2]["SNP"].values)
                overlap = len(snps1 & snps2)
                matrix.loc[a1, a2] += overlap
                if a1 != a2:
                    matrix.loc[a2, a1] += overlap

    return matrix


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

COLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
          "#a65628", "#f781bf", "#999999"]


def plot_venn_or_upset(pheno: str, alg_dict: dict, out_dir: str):
    """UpSet-style bar chart of method intersections."""
    algs = list(alg_dict.keys())
    snp_sets = {alg: set(df["SNP"].values) for alg, df in alg_dict.items()}
    all_snps = sorted(set.union(*snp_sets.values()) if snp_sets else set())

    if not all_snps:
        return

    # Build presence matrix
    mat = np.array([[1 if snp in snp_sets[alg] else 0 for alg in algs]
                    for snp in all_snps])

    # Unique intersection patterns
    patterns, counts = np.unique(mat, axis=0, return_counts=True)
    order = np.argsort(-counts)
    patterns = patterns[order]
    counts   = counts[order]

    fig, axes = plt.subplots(2, 1, figsize=(max(8, len(patterns) * 0.6 + 2), 7),
                              gridspec_kw={"height_ratios": [3, 1]})
    fig.suptitle(f"Method overlap – {pheno}", fontsize=13, fontweight="bold")

    # Top: bar chart of intersection sizes
    ax_bar = axes[0]
    bar_colors = [COLORS[i % len(COLORS)] for i in range(len(patterns))]
    ax_bar.bar(range(len(counts)), counts, color=bar_colors, edgecolor="white")
    ax_bar.set_ylabel("# significant SNPs")
    ax_bar.set_xticks([])
    ax_bar.set_title("Intersection sizes")

    # Bottom: dot matrix showing which methods are in each intersection
    ax_dot = axes[1]
    ax_dot.set_xlim(-0.5, len(patterns) - 0.5)
    ax_dot.set_ylim(-0.5, len(algs) - 0.5)
    ax_dot.set_yticks(range(len(algs)))
    ax_dot.set_yticklabels(algs, fontsize=8)
    ax_dot.set_xticks([])
    ax_dot.set_xlabel("Intersection (each column = one combination)")

    for col_i, pat in enumerate(patterns):
        active = [j for j, v in enumerate(pat) if v]
        inactive = [j for j, v in enumerate(pat) if not v]
        for j in active:
            ax_dot.plot(col_i, j, "o", color="#222222", markersize=8)
        for j in inactive:
            ax_dot.plot(col_i, j, "o", color="#dddddd", markersize=8)
        if len(active) > 1:
            ax_dot.plot([col_i, col_i], [min(active), max(active)],
                        color="#222222", linewidth=2, zorder=0)

    plt.tight_layout()
    out_path = os.path.join(out_dir, f"venn_{_safe_fname(pheno)}.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("  Saved: %s", out_path)


def plot_heatmap(pheno: str, table: pd.DataFrame, out_dir: str):
    """SNP × method heatmap for top SNPs found in multiple methods."""
    alg_cols = [c for c in table.columns
                if c not in ("SNP", "Chr", "Pos", "methods_count")]
    if not alg_cols:
        return

    # Show top 50 SNPs by methods_count
    top = table.nlargest(50, "methods_count")
    mat = top[alg_cols].values.astype(float)

    snp_labels = top["SNP"].values
    if "Chr" in top.columns and "Pos" in top.columns:
        snp_labels = [f"{row.SNP} ({row.Chr}:{row.Pos})"
                      for row in top.itertuples()]

    fig, ax = plt.subplots(figsize=(max(5, len(alg_cols) * 1.4),
                                    max(4, len(top) * 0.3 + 2)))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_xticks(range(len(alg_cols)))
    ax.set_xticklabels(alg_cols, rotation=35, ha="right", fontsize=9)
    ax.set_yticks(range(len(snp_labels)))
    ax.set_yticklabels(snp_labels, fontsize=7)
    ax.set_title(f"SNP × Method significance – {pheno}\n(top {len(top)} SNPs by method count)",
                 fontsize=11)
    plt.colorbar(im, ax=ax, label="Significant (1=yes)")
    plt.tight_layout()
    out_path = os.path.join(out_dir, f"heatmap_{_safe_fname(pheno)}.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("  Saved: %s", out_path)


def plot_manhattan_comparison(pheno: str, alg_dict: dict, out_dir: str):
    """Stacked mini-Manhattan plots, one row per algorithm."""
    algs = list(alg_dict.keys())
    if not algs:
        return

    fig, axes = plt.subplots(len(algs), 1,
                              figsize=(14, 2.5 * len(algs)),
                              sharex=False)
    if len(algs) == 1:
        axes = [axes]

    fig.suptitle(f"Manhattan comparison – {pheno}", fontsize=13, fontweight="bold")

    for ax, alg in zip(axes, algs):
        df = alg_dict[alg]
        if df.empty or "Chr" not in df.columns or "Pos" not in df.columns:
            ax.text(0.5, 0.5, f"{alg}\n(no position data)",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_ylabel(alg, fontsize=8)
            continue

        df = df.copy()
        df["Chr"] = pd.to_numeric(df["Chr"], errors="coerce")
        df["Pos"] = pd.to_numeric(df["Pos"], errors="coerce")
        df = df.dropna(subset=["Chr", "Pos"])
        df = df.sort_values(["Chr", "Pos"])

        # Cumulative position for x-axis
        chrom_offsets = {}
        offset = 0
        for chrom, grp in df.groupby("Chr", sort=True):
            chrom_offsets[chrom] = offset
            offset += grp["Pos"].max() + 1_000_000

        df["cum_pos"] = df.apply(
            lambda r: chrom_offsets.get(r["Chr"], 0) + r["Pos"], axis=1)

        # Score or -log10(p)
        if "PValue" in df.columns:
            pvals = pd.to_numeric(df["PValue"], errors="coerce").clip(lower=1e-300)
            y = -np.log10(pvals)
            ylabel = "-log10(p)"
        elif "Score" in df.columns:
            y = pd.to_numeric(df["Score"], errors="coerce").fillna(0)
            ylabel = "Score"
        else:
            ax.text(0.5, 0.5, f"{alg}\n(no score column)",
                    ha="center", va="center", transform=ax.transAxes)
            continue

        # Color by chromosome
        chrom_nums = df["Chr"].astype(int)
        point_colors = [COLORS[c % len(COLORS)] for c in chrom_nums]

        ax.scatter(df["cum_pos"], y, c=point_colors, s=6, alpha=0.7, linewidths=0)
        ax.set_ylabel(f"{alg}\n{ylabel}", fontsize=7)

        # Chromosome tick labels
        chrom_mids = {}
        for chrom, grp in df.groupby("Chr", sort=True):
            cum = grp["cum_pos"]
            chrom_mids[chrom] = (cum.min() + cum.max()) / 2
        ax.set_xticks(list(chrom_mids.values()))
        ax.set_xticklabels([str(int(c)) for c in chrom_mids.keys()], fontsize=6)
        ax.set_xlabel("Chromosome", fontsize=7)

    plt.tight_layout()
    out_path = os.path.join(out_dir, f"manhattan_comparison_{_safe_fname(pheno)}.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("  Saved: %s", out_path)


def plot_pairwise_heatmap(overlap_df: pd.DataFrame, out_dir: str):
    """Heatmap of pairwise method overlap counts."""
    if overlap_df.empty:
        return
    fig, ax = plt.subplots(figsize=(max(5, len(overlap_df) * 1.2),
                                    max(4, len(overlap_df) * 1.0)))
    vals = overlap_df.values.astype(float)
    im = ax.imshow(vals, cmap="Blues")
    ax.set_xticks(range(len(overlap_df.columns)))
    ax.set_xticklabels(overlap_df.columns, rotation=35, ha="right", fontsize=9)
    ax.set_yticks(range(len(overlap_df.index)))
    ax.set_yticklabels(overlap_df.index, fontsize=9)
    ax.set_title("Pairwise method overlap\n(# shared significant SNPs, all phenotypes)",
                 fontsize=11)
    plt.colorbar(im, ax=ax)
    for i in range(len(overlap_df.index)):
        for j in range(len(overlap_df.columns)):
            ax.text(j, i, str(int(vals[i, j])), ha="center", va="center", fontsize=8)
    plt.tight_layout()
    out_path = os.path.join(out_dir, "pairwise_method_overlap_heatmap.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("  Saved: %s", out_path)


def _safe_fname(s: str):
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description=(
            "Compare GWAS results across methods and phenotypes.\n"
            "Reads all gwas_results.csv files produced by gwastic_batch.py."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--results-dir", "-r", required=True,
        help="Root results folder (contains the timestamped sub-folders).",
    )
    p.add_argument(
        "--output", "-o", default=None,
        help="Where to save comparison files. Defaults to --results-dir/comparison/.",
    )
    p.add_argument(
        "--pval-threshold", type=float, default=5e-8,
        help="P-value threshold for significance (default: 5e-8 genome-wide).",
    )
    p.add_argument(
        "--top-n", type=int, default=100,
        help="If no p-value column found, use top-N SNPs by score (default: 100).",
    )
    p.add_argument(
        "--min-methods", type=int, default=2,
        help="Minimum number of methods a SNP must appear in for consensus (default: 2).",
    )
    p.add_argument(
        "--skip-plots", action="store_true",
        help="Skip all plot generation.",
    )
    return p


def main():
    args = build_parser().parse_args()

    results_dir = os.path.abspath(args.results_dir)
    out_dir     = os.path.abspath(args.output or os.path.join(results_dir, "comparison"))
    plot_dir    = os.path.join(out_dir, "plots")
    os.makedirs(out_dir,  exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    log.info("GWAStic Compare")
    log.info("  Results dir     : %s", results_dir)
    log.info("  Output dir      : %s", out_dir)
    log.info("  P-val threshold : %s", args.pval_threshold)
    log.info("  Top-N fallback  : %d", args.top_n)
    log.info("  Min methods     : %d", args.min_methods)

    # ------------------------------------------------------------------
    # Discover & load
    # ------------------------------------------------------------------
    log.info("\nDiscovering runs …")
    runs = discover_runs(results_dir)
    if not runs:
        log.error("No gwas_results.csv files found under %s", results_dir)
        sys.exit(1)
    log.info("  Found %d run(s) across %d phenotype(s)",
             len(runs),
             len({r["phenotype"] for r in runs}))

    log.info("\nLoading significant SNPs …")
    data = load_all_runs(runs, args.pval_threshold, args.top_n)

    if not data:
        log.error("No significant SNPs loaded. Check --pval-threshold or --top-n.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Build tables
    # ------------------------------------------------------------------
    log.info("\nBuilding comparison tables …")
    comp_tables = build_comparison_table(data)
    summary_df  = build_summary(data)
    consensus   = build_consensus(comp_tables, args.min_methods)
    overlap_df  = build_pairwise_overlap(data)

    # ------------------------------------------------------------------
    # Save CSVs
    # ------------------------------------------------------------------
    all_frames = []
    for pheno, df in comp_tables.items():
        df_out = df.copy()
        df_out.insert(0, "Phenotype", pheno)
        all_frames.append(df_out)

    if all_frames:
        all_sig = pd.concat(all_frames, ignore_index=True)
        p = os.path.join(out_dir, "comparison_all_significant_snps.csv")
        all_sig.to_csv(p, index=False)
        log.info("Saved: %s  (%d rows)", p, len(all_sig))

    p = os.path.join(out_dir, "comparison_summary.csv")
    summary_df.to_csv(p, index=False)
    log.info("Saved: %s", p)

    if not consensus.empty:
        p = os.path.join(out_dir, "comparison_consensus_snps.csv")
        consensus.to_csv(p, index=False)
        log.info("Saved: %s  (%d SNPs in >= %d methods)",
                 p, len(consensus), args.min_methods)
    else:
        log.info("No consensus SNPs found (>= %d methods).", args.min_methods)

    p = os.path.join(out_dir, "comparison_method_overlap.csv")
    overlap_df.to_csv(p)
    log.info("Saved: %s", p)

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    if not args.skip_plots:
        log.info("\nGenerating plots …")

        log.info("  Pairwise overlap heatmap …")
        plot_pairwise_heatmap(overlap_df, plot_dir)

        for pheno, alg_dict in data.items():
            log.info("  Phenotype: %s", pheno)
            plot_venn_or_upset(pheno, alg_dict, plot_dir)
            if pheno in comp_tables:
                plot_heatmap(pheno, comp_tables[pheno], plot_dir)
            plot_manhattan_comparison(pheno, alg_dict, plot_dir)

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------
    log.info("\n" + "=" * 70)
    log.info("SUMMARY")
    log.info("=" * 70)
    print(summary_df.to_string(index=False))
    log.info("\nConsensus SNPs (>= %d methods): %d", args.min_methods, len(consensus))
    log.info("All output written to: %s", out_dir)


if __name__ == "__main__":
    main()
