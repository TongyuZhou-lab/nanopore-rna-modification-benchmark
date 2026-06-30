#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SWARM/m6Anet deep analysis with fully command-line configurable paths.

Main analyses
-------------
1. SWARM probability and stoichiometry threshold scans.
2. Complementarity analysis between SWARM and m6Anet.
3. Optional f5c-SWARM vs nanopolish-SWARM comparison.

All input and output paths are provided through command-line arguments.
No benchmark-specific file paths are hard-coded in this script.
"""

import argparse
import re
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def log(message: str) -> None:
    """Print a flushed progress message."""
    print(f"[SWARM-DEEP] {message}", flush=True)


def check(path: Optional[str], name: str, required: bool = True) -> Optional[Path]:
    """Validate an input path and return it as a Path object."""
    if path is None or str(path).strip() == "":
        if required:
            raise FileNotFoundError(f"Missing required path: {name}")
        return None

    p = Path(path)
    if not p.exists():
        if required:
            raise FileNotFoundError(f"{name} not found: {p}")
        log(f"Optional file not found; skipping {name}: {p}")
        return None
    return p


def read_bed(path: str) -> pd.DataFrame:
    """Read a BED-like file and build a chrom:start site key from the first three columns."""
    path_obj = Path(path)
    df = pd.read_csv(path_obj, sep="\t", header=None, comment="#")
    if df.shape[1] < 3:
        raise ValueError(f"BED-like file needs at least 3 columns: {path_obj}")

    df = df.iloc[:, :3].copy()
    df.columns = ["chrom", "start", "end"]
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["site_key"] = df["chrom"] + ":" + df["start"].astype(str)
    return df.drop_duplicates("site_key")


def read_scorebed(path: str, score_col: str) -> pd.DataFrame:
    """Read a score BED file and keep the maximum score for duplicated sites."""
    path_obj = Path(path)
    df = pd.read_csv(path_obj, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        raise ValueError(f"Score BED file needs at least 4 columns: {path_obj}")

    df = df.iloc[:, :4].copy()
    df.columns = ["chrom", "start", "end", score_col]
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df["site_key"] = df["chrom"] + ":" + df["start"].astype(str)

    df = df.groupby(["chrom", "start", "end", "site_key"], as_index=False)[score_col].max()
    return df


def infer_col(cols: Iterable[str], patterns: Iterable[str], avoid: Iterable[str] = ()) -> Optional[str]:
    """Infer a column name by regex patterns while avoiding coordinate or label columns."""
    avoid_lower = {str(x).lower() for x in avoid}
    for pattern in patterns:
        hits = []
        for col in cols:
            col_lower = str(col).lower()
            if col_lower in avoid_lower:
                continue
            if re.search(pattern, col_lower):
                hits.append(col)
        if hits:
            return hits[0]
    return None


def read_site_table(path: str, prob_name: str = "prob", stoich_name: str = "stoich") -> pd.DataFrame:
    """
    Read a site-level/detail table with a header.

    The function tries to infer coordinate columns, probability/score columns,
    and stoichiometry/modification-ratio columns automatically. It supports
    either explicit coordinate columns or a site_id-like column.
    """
    path_obj = Path(path)
    df = pd.read_csv(path_obj, sep="\t")
    if df.shape[1] == 1:
        df = pd.read_csv(path_obj)

    cols = list(df.columns)
    lower = {str(c).lower(): c for c in cols}

    chrom_col = None
    for candidate in ["chrom", "chr", "contig", "reference"]:
        if candidate in lower:
            chrom_col = lower[candidate]
            break

    start_col = None
    for candidate in ["start", "pos", "position", "genome_pos", "genomic_position"]:
        if candidate in lower:
            start_col = lower[candidate]
            break

    end_col = None
    for candidate in ["end", "stop"]:
        if candidate in lower:
            end_col = lower[candidate]
            break

    if chrom_col is None or start_col is None:
        sid_col = None
        for candidate in ["site_id", "site", "site_key"]:
            if candidate in lower:
                sid_col = lower[candidate]
                break

        if sid_col is not None:
            tmp = df[sid_col].astype(str).str.split(":", expand=True)
            if tmp.shape[1] >= 2:
                df["chrom"] = tmp[0]
                df["start"] = tmp[1].astype(int)
                df["end"] = tmp[2].astype(int) if tmp.shape[1] >= 3 else df["start"] + 1
            else:
                raise ValueError(f"Cannot parse site_id column from {path_obj}")
        else:
            if df.shape[1] >= 3:
                df = df.rename(columns={df.columns[0]: "chrom", df.columns[1]: "start", df.columns[2]: "end"})
            else:
                raise ValueError(f"Cannot infer coordinate columns in {path_obj}. Columns={cols}")
    else:
        df = df.rename(columns={chrom_col: "chrom", start_col: "start"})
        if end_col is not None:
            df = df.rename(columns={end_col: "end"})
        else:
            df["end"] = pd.to_numeric(df["start"], errors="coerce").astype(int) + 1

    prob_col = infer_col(
        df.columns,
        [r"probability", r"prob", r"score", r"prediction"],
        avoid=("chrom", "chr", "start", "end", "pos", "position", "label"),
    )
    stoich_col = infer_col(
        df.columns,
        [r"stoich", r"mod.*ratio", r"ratio"],
        avoid=("chrom", "chr", "start", "end", "pos", "position", "label"),
    )

    df["chrom"] = df["chrom"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype(int)
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype(int)
    df["site_key"] = df["chrom"] + ":" + df["start"].astype(str)

    keep = ["chrom", "start", "end", "site_key"]
    if prob_col is not None:
        df[prob_name] = pd.to_numeric(df[prob_col], errors="coerce")
        keep.append(prob_name)
    if stoich_col is not None and stoich_col != prob_col:
        df[stoich_name] = pd.to_numeric(df[stoich_col], errors="coerce")
        keep.append(stoich_name)

    out = df[keep].copy()
    agg = {c: "max" for c in out.columns if c not in ["chrom", "start", "end", "site_key"]}
    if agg:
        out = out.groupby(["chrom", "start", "end", "site_key"], as_index=False).agg(agg)
    else:
        out = out.drop_duplicates("site_key")
    return out


def label_common_sites(
    common: pd.DataFrame,
    positive: pd.DataFrame,
    exclude: pd.DataFrame,
    exclude_buffer: int = 2,
) -> pd.DataFrame:
    """
    Label common sites using the positive reference and remove ambiguous negative sites.

    A site is labeled as positive if its site_key appears in the positive BED.
    A non-positive site is removed as ambiguous if it is close to an excluded site
    within the configured coordinate buffer.
    """
    pos_set = set(positive["site_key"])
    common = common.copy()
    common["label"] = common["site_key"].isin(pos_set).astype(int)

    exc_by_chr = {}
    for chrom, sub in exclude.groupby("chrom"):
        exc_by_chr[chrom] = np.sort(sub["start"].astype(int).values)

    ambiguous = []
    for _, row in common.iterrows():
        if row["label"] == 1:
            ambiguous.append(False)
            continue

        arr = exc_by_chr.get(row["chrom"], None)
        if arr is None or len(arr) == 0:
            ambiguous.append(False)
            continue

        x = int(row["start"])
        idx = np.searchsorted(arr, x)
        near = False
        for j in [idx - 1, idx, idx + 1]:
            if 0 <= j < len(arr) and abs(int(arr[j]) - x) <= exclude_buffer:
                near = True
                break
        ambiguous.append(near)

    common["ambiguous"] = ambiguous
    common = common[~common["ambiguous"]].copy()
    common["label"] = common["label"].astype(int)
    return common


def metrics(y, pred) -> dict:
    """Calculate binary classification metrics from labels and predictions."""
    y = np.asarray(y).astype(int)
    pred = np.asarray(pred).astype(bool)

    tp = int(((y == 1) & pred).sum())
    fp = int(((y == 0) & pred).sum())
    fn = int(((y == 1) & (~pred)).sum())
    tn = int(((y == 0) & (~pred)).sum())

    precision = tp / (tp + fp) if tp + fp else np.nan
    recall = tp / (tp + fn) if tp + fn else np.nan
    specificity = tn / (tn + fp) if tn + fp else np.nan
    f1 = 2 * precision * recall / (precision + recall) if precision == precision and recall == recall and precision + recall > 0 else np.nan
    acc = (tp + tn) / (tp + fp + fn + tn) if tp + fp + fn + tn else np.nan

    return {
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "Precision": precision,
        "Recall": recall,
        "Specificity": specificity,
        "F1": f1,
        "Accuracy": acc,
        "Predicted_positive": int(pred.sum()),
        "N": int(len(y)),
    }


def roc_pr_auc(y, score) -> dict:
    """Compute AUROC and AUPRC without requiring scikit-learn."""
    y = np.asarray(y).astype(int)
    score = np.asarray(score).astype(float)
    mask = np.isfinite(score)
    y = y[mask]
    score = score[mask]

    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return {"AUROC": np.nan, "AUPRC": np.nan, "n": len(y), "n_pos": n_pos, "n_neg": n_neg}

    order = np.argsort(-score, kind="mergesort")
    y = y[order]
    score = score[order]

    distinct = np.r_[np.where(np.diff(score) != 0)[0], len(score) - 1]
    tp = np.cumsum(y == 1)[distinct]
    fp = np.cumsum(y == 0)[distinct]

    tpr = tp / n_pos
    fpr = fp / n_neg
    precision = tp / np.maximum(tp + fp, 1)
    recall = tpr

    auroc = float(np.trapz(np.r_[0, tpr, 1], np.r_[0, fpr, 1]))
    auprc = float(np.sum((np.r_[0, recall][1:] - np.r_[0, recall][:-1]) * precision))

    return {"AUROC": auroc, "AUPRC": auprc, "n": len(y), "n_pos": n_pos, "n_neg": n_neg}


def scan_threshold(y, val, thresholds) -> pd.DataFrame:
    """Scan one-dimensional thresholds and return metric rows."""
    rows = []
    for threshold in thresholds:
        m = metrics(y, val >= threshold)
        m["threshold"] = float(threshold)
        rows.append(m)
    return pd.DataFrame(rows)


def plot_scan(df: pd.DataFrame, out: Path, title: str) -> None:
    """Plot threshold-dependent precision, recall, F1, and predicted-positive counts."""
    fig, ax1 = plt.subplots(figsize=(6.5, 4.2))
    for metric_name in ["Precision", "Recall", "F1"]:
        ax1.plot(df["threshold"], df[metric_name], label=metric_name, linewidth=1.8)
    ax1.set_ylim(0, 1.02)
    ax1.set_xlabel("threshold")
    ax1.set_ylabel("metric")
    ax1.grid(alpha=0.25)

    ax2 = ax1.twinx()
    ax2.plot(df["threshold"], df["Predicted_positive"], linestyle="--", linewidth=1.2, label="Predicted positive")
    ax2.set_ylabel("predicted positive sites")

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=8, loc="best")
    ax1.set_title(title)
    fig.tight_layout()
    fig.savefig(str(out) + ".png", dpi=300)
    fig.savefig(str(out) + ".svg")
    plt.close(fig)


def plot_heatmap(df: pd.DataFrame, value: str, out: Path, title: str) -> None:
    """Plot a heatmap for a joint probability-stoichiometry threshold scan."""
    mat = df.pivot_table(index="stoich_threshold", columns="prob_threshold", values=value, aggfunc="max")
    fig, ax = plt.subplots(figsize=(7.5, 5.3))
    im = ax.imshow(mat.values, aspect="auto", origin="lower")
    ax.set_title(title)
    ax.set_xlabel("probability threshold")
    ax.set_ylabel("stoichiometry threshold")

    xs = np.linspace(0, mat.shape[1] - 1, min(6, mat.shape[1])).astype(int)
    ys = np.linspace(0, mat.shape[0] - 1, min(6, mat.shape[0])).astype(int)
    ax.set_xticks(xs)
    ax.set_yticks(ys)
    ax.set_xticklabels([f"{mat.columns[i]:.4g}" for i in xs], rotation=45, ha="right")
    ax.set_yticklabels([f"{mat.index[i]:.4g}" for i in ys])

    cb = fig.colorbar(im, ax=ax)
    cb.set_label(value)
    fig.tight_layout()
    fig.savefig(str(out) + ".png", dpi=300)
    fig.savefig(str(out) + ".svg")
    plt.close(fig)


def plot_bar_metrics(df: pd.DataFrame, out: Path, title: str) -> None:
    """Plot precision, recall, and F1 for multiple call-set strategies."""
    fig, ax = plt.subplots(figsize=(7, 4.2))
    x = np.arange(len(df))
    width = 0.25
    for i, metric_name in enumerate(["Precision", "Recall", "F1"]):
        ax.bar(x + (i - 1) * width, df[metric_name], width=width, label=metric_name)
    ax.set_xticks(x)
    ax.set_xticklabels(df["strategy"], rotation=25, ha="right")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("metric")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(str(out) + ".png", dpi=300)
    fig.savefig(str(out) + ".svg")
    plt.close(fig)


def plot_composition(df: pd.DataFrame, out: Path, title: str) -> None:
    """Plot reference-positive and pseudo-negative composition for each call-set category."""
    fig, ax = plt.subplots(figsize=(7, 4.2))
    x = np.arange(len(df))
    ax.bar(x, df["TP"], label="reference positive")
    ax.bar(x, df["FP"], bottom=df["TP"], label="pseudo-negative")
    ax.set_xticks(x)
    ax.set_xticklabels(df["set"], rotation=25, ha="right")
    ax.set_ylabel("sites")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(str(out) + ".png", dpi=300)
    fig.savefig(str(out) + ".svg")
    plt.close(fig)


def plot_scatter(df: pd.DataFrame, x: str, y: str, out: Path, title: str) -> None:
    """Plot a scatter comparison and report Pearson/Spearman correlations in the title."""
    sub = df[[x, y]].dropna()
    if len(sub) == 0:
        log(f"No valid paired values for scatter plot: {x} vs {y}")
        return

    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    ax.scatter(sub[x], sub[y], s=8, alpha=0.35)
    lo = min(sub[x].min(), sub[y].min())
    hi = max(sub[x].max(), sub[y].max())
    ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1)
    pearson = sub[x].corr(sub[y], method="pearson") if len(sub) > 1 else np.nan
    spearman = sub[x].corr(sub[y], method="spearman") if len(sub) > 1 else np.nan
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(f"{title}\nPearson={pearson:.3f}, Spearman={spearman:.3f}, n={len(sub)}")
    fig.tight_layout()
    fig.savefig(str(out) + ".png", dpi=300)
    fig.savefig(str(out) + ".svg")
    plt.close(fig)


def build_table(args: argparse.Namespace) -> pd.DataFrame:
    """Build the labeled evaluation table with SWARM and m6Anet scores."""
    common = read_bed(args.common_bed)
    pos = read_bed(args.positive_bed)
    exc = read_bed(args.exclude_bed)

    labeled = label_common_sites(common, pos, exc, exclude_buffer=args.exclude_buffer)

    swarm = read_scorebed(args.swarm_scorebed, "swarm_prob")
    m6a = read_scorebed(args.m6anet_scorebed, "m6anet_prob")

    df = labeled.merge(swarm[["site_key", "swarm_prob"]], on="site_key", how="left")
    df = df.merge(m6a[["site_key", "m6anet_prob"]], on="site_key", how="left")

    if args.swarm_stoich_table and Path(args.swarm_stoich_table).exists():
        sto = read_site_table(args.swarm_stoich_table, prob_name="swarm_prob_detail", stoich_name="swarm_stoich")
        keep = ["site_key"]
        if "swarm_stoich" in sto.columns:
            keep.append("swarm_stoich")
        if "swarm_prob_detail" in sto.columns:
            keep.append("swarm_prob_detail")

        sto = sto[keep].drop_duplicates("site_key")
        df = df.merge(sto, on="site_key", how="left")
        if "swarm_prob_detail" in df.columns:
            df["swarm_prob"] = df["swarm_prob"].fillna(df["swarm_prob_detail"])

    return df


def run_threshold(df: pd.DataFrame, outdir: Path, args: argparse.Namespace) -> None:
    """Run SWARM probability, stoichiometry, and joint threshold scans."""
    sub = df[df["swarm_prob"].notna()].copy()
    if len(sub) == 0:
        log("No SWARM probability values detected. Skipping SWARM threshold scans.")
        return

    y = sub["label"].values
    prob_grid = np.unique(np.r_[np.linspace(0, 1, 101), np.linspace(0.90, 0.9999, 101), [args.swarm_prob_cutoff]])
    prob_grid.sort()

    prob_scan = scan_threshold(y, sub["swarm_prob"].values, prob_grid)
    prob_scan.to_csv(outdir / "SWARM_probability_threshold_scan.tsv", sep="\t", index=False)
    plot_scan(prob_scan, outdir / "SWARM_probability_threshold_scan", "SWARM probability threshold scan")

    if "swarm_stoich" not in sub.columns or sub["swarm_stoich"].notna().sum() == 0:
        log("No SWARM stoichiometry detected. Skipping stoichiometry and joint scans.")
        return

    sub2 = sub[sub["swarm_stoich"].notna()].copy()
    y2 = sub2["label"].values

    sto_grid = np.unique(np.r_[np.linspace(0, 1, 101), [args.swarm_stoich_cutoff]])
    sto_grid.sort()
    sto_scan = scan_threshold(y2, sub2["swarm_stoich"].values, sto_grid)
    sto_scan.to_csv(outdir / "SWARM_stoichiometry_threshold_scan.tsv", sep="\t", index=False)
    plot_scan(sto_scan, outdir / "SWARM_stoichiometry_threshold_scan", "SWARM stoichiometry threshold scan")

    joint_prob_grid = np.unique(np.r_[np.linspace(0.90, 0.9999, 81), [args.swarm_prob_cutoff]])
    joint_sto_grid = np.unique(np.r_[np.linspace(0, 0.6, 61), [args.swarm_stoich_cutoff]])
    joint_prob_grid.sort()
    joint_sto_grid.sort()

    rows = []
    p = sub2["swarm_prob"].values
    s = sub2["swarm_stoich"].values

    for prob_threshold in joint_prob_grid:
        for stoich_threshold in joint_sto_grid:
            pred = (p >= prob_threshold) & (s > stoich_threshold)
            m = metrics(y2, pred)
            m["prob_threshold"] = float(prob_threshold)
            m["stoich_threshold"] = float(stoich_threshold)
            rows.append(m)

    joint = pd.DataFrame(rows)
    joint.to_csv(outdir / "SWARM_joint_probability_stoichiometry_threshold_scan.tsv", sep="\t", index=False)
    plot_heatmap(joint, "F1", outdir / "SWARM_joint_threshold_F1_heatmap", "SWARM joint threshold scan: F1")
    plot_heatmap(joint, "Precision", outdir / "SWARM_joint_threshold_Precision_heatmap", "SWARM joint threshold scan: Precision")


def run_complementarity(df: pd.DataFrame, outdir: Path, args: argparse.Namespace) -> None:
    """Compare SWARM, m6Anet, their intersection, and their union."""
    valid = df["swarm_prob"].notna() & df["m6anet_prob"].notna()
    if "swarm_stoich" in df.columns:
        valid = valid & df["swarm_stoich"].notna()

    sub = df[valid].copy()
    if len(sub) == 0:
        log("No common callable sites for SWARM and m6Anet. Skipping complementarity analysis.")
        return

    y = sub["label"].values

    if "swarm_stoich" in sub.columns:
        swarm_pos = (sub["swarm_prob"] >= args.swarm_prob_cutoff) & (sub["swarm_stoich"] > args.swarm_stoich_cutoff)
    else:
        swarm_pos = sub["swarm_prob"] >= args.swarm_prob_cutoff

    m6anet_pos = sub["m6anet_prob"] >= args.m6anet_cutoff

    strategies = {
        "SWARM": swarm_pos.values,
        "m6Anet": m6anet_pos.values,
        "intersection": (swarm_pos & m6anet_pos).values,
        "union": (swarm_pos | m6anet_pos).values,
    }

    rows = []
    for name, pred in strategies.items():
        m = metrics(y, pred)
        m["strategy"] = name
        rows.append(m)

    met = pd.DataFrame(rows)
    met = met[["strategy", "N", "Predicted_positive", "TP", "FP", "FN", "TN", "Precision", "Recall", "Specificity", "F1", "Accuracy"]]
    met.to_csv(outdir / "SWARM_m6Anet_strategy_metrics.tsv", sep="\t", index=False)
    plot_bar_metrics(met, outdir / "SWARM_m6Anet_strategy_metrics", "SWARM and m6Anet strategy metrics")

    cat = np.where(
        swarm_pos & m6anet_pos,
        "intersection",
        np.where(swarm_pos & (~m6anet_pos), "SWARM_only", np.where((~swarm_pos) & m6anet_pos, "m6Anet_only", "neither")),
    )

    sub["category"] = cat

    comp_rows = []
    for category in ["intersection", "SWARM_only", "m6Anet_only", "neither"]:
        ss = sub[sub["category"] == category]
        yy = ss["label"].values
        tp = int((yy == 1).sum())
        fp = int((yy == 0).sum())
        comp_rows.append(
            {
                "set": category,
                "n_sites": int(len(ss)),
                "TP": tp,
                "FP": fp,
                "reference_positive_fraction": tp / len(ss) if len(ss) else np.nan,
            }
        )

    comp = pd.DataFrame(comp_rows)
    comp.to_csv(outdir / "SWARM_m6Anet_set_composition.tsv", sep="\t", index=False)
    plot_composition(comp, outdir / "SWARM_m6Anet_set_composition", "SWARM/m6Anet positive-call composition")

    keep = ["chrom", "start", "end", "site_key", "label", "swarm_prob", "m6anet_prob", "category"]
    if "swarm_stoich" in sub.columns:
        keep.append("swarm_stoich")
    sub[keep].to_csv(outdir / "SWARM_m6Anet_labeled_site_categories.tsv", sep="\t", index=False)


def run_nanopolish(df: pd.DataFrame, outdir: Path, args: argparse.Namespace) -> None:
    """Optionally compare f5c-based SWARM results against nanopolish-based SWARM results."""
    if not args.swarm_nanopolish_sitelevel and not args.swarm_nanopolish_scorebed:
        log("No nanopolish SWARM file provided. Skipping f5c vs nanopolish analysis.")
        return

    if args.swarm_nanopolish_sitelevel:
        nano = read_site_table(args.swarm_nanopolish_sitelevel, prob_name="nano_prob", stoich_name="nano_stoich")
    else:
        nano = read_scorebed(args.swarm_nanopolish_scorebed, "nano_prob")

    sub = df.merge(nano, on="site_key", how="inner", suffixes=("", "_nano"))
    if len(sub) == 0:
        log("No overlap between f5c SWARM and nanopolish SWARM. Check coordinates or shifts.")
        return

    counts = pd.DataFrame(
        [
            {"source": "f5c_SWARM_eval_available", "n_sites": int(df["swarm_prob"].notna().sum())},
            {"source": "nanopolish_SWARM_available", "n_sites": int(nano["site_key"].nunique())},
            {"source": "common_f5c_nanopolish", "n_sites": int(len(sub))},
        ]
    )
    counts.to_csv(outdir / "SWARM_f5c_nanopolish_callable_counts.tsv", sep="\t", index=False)

    if "nano_prob" in sub.columns:
        plot_scatter(sub, "swarm_prob", "nano_prob", outdir / "SWARM_f5c_vs_nanopolish_probability_scatter", "SWARM probability: f5c vs nanopolish")

    if "swarm_stoich" in sub.columns and "nano_stoich" in sub.columns:
        plot_scatter(sub, "swarm_stoich", "nano_stoich", outdir / "SWARM_f5c_vs_nanopolish_stoichiometry_scatter", "SWARM stoichiometry: f5c vs nanopolish")

    y = sub["label"].values

    auc_rows = []
    if "swarm_prob" in sub.columns:
        result = roc_pr_auc(y, sub["swarm_prob"].values)
        result["source"] = "f5c_SWARM"
        auc_rows.append(result)
    if "nano_prob" in sub.columns:
        result = roc_pr_auc(y, sub["nano_prob"].values)
        result["source"] = "nanopolish_SWARM"
        auc_rows.append(result)

    pd.DataFrame(auc_rows).to_csv(outdir / "SWARM_f5c_nanopolish_AUROC_AUPRC_common_sites.tsv", sep="\t", index=False)

    fixed_rows = []
    if "swarm_stoich" in sub.columns:
        pred_f5c = (sub["swarm_prob"] >= args.swarm_prob_cutoff) & (sub["swarm_stoich"] > args.swarm_stoich_cutoff)
    else:
        pred_f5c = sub["swarm_prob"] >= args.swarm_prob_cutoff

    m = metrics(y, pred_f5c.values)
    m["source"] = "f5c_SWARM"
    fixed_rows.append(m)

    if "nano_prob" in sub.columns:
        if "nano_stoich" in sub.columns:
            pred_nano = (sub["nano_prob"] >= args.swarm_prob_cutoff) & (sub["nano_stoich"] > args.swarm_stoich_cutoff)
        else:
            pred_nano = sub["nano_prob"] >= args.swarm_prob_cutoff
        m = metrics(y, pred_nano.values)
        m["source"] = "nanopolish_SWARM"
        fixed_rows.append(m)

    pd.DataFrame(fixed_rows).to_csv(outdir / "SWARM_f5c_nanopolish_fixed_cutoff_metrics_common_sites.tsv", sep="\t", index=False)

    keep = [c for c in ["chrom", "start", "end", "site_key", "label", "swarm_prob", "swarm_stoich", "nano_prob", "nano_stoich"] if c in sub.columns]
    sub[keep].to_csv(outdir / "SWARM_f5c_nanopolish_common_site_scores.tsv", sep="\t", index=False)


def write_run_summary(df: pd.DataFrame, outdir: Path, args: argparse.Namespace) -> None:
    """Write a run-specific summary file into the output directory."""
    lines = []
    lines.append("SWARM/m6Anet deep analysis run summary")
    lines.append("")
    lines.append(f"evaluation sites: {len(df)}")
    lines.append(f"reference positives: {int((df['label'] == 1).sum())}")
    lines.append(f"pseudo-negatives: {int((df['label'] == 0).sum())}")
    lines.append(f"SWARM probability available: {int(df['swarm_prob'].notna().sum())}")
    lines.append(f"m6Anet score available: {int(df['m6anet_prob'].notna().sum())}")
    if "swarm_stoich" in df.columns:
        lines.append(f"SWARM stoichiometry available: {int(df['swarm_stoich'].notna().sum())}")
    else:
        lines.append("SWARM stoichiometry available: 0")

    lines.append("")
    lines.append("cutoffs:")
    lines.append(f"SWARM probability >= {args.swarm_prob_cutoff}")
    lines.append(f"SWARM stoichiometry > {args.swarm_stoich_cutoff}")
    lines.append(f"m6Anet probability_modified >= {args.m6anet_cutoff}")

    lines.append("")
    lines.append("input paths:")
    for key in [
        "common_bed",
        "positive_bed",
        "exclude_bed",
        "swarm_scorebed",
        "swarm_stoich_table",
        "m6anet_scorebed",
        "swarm_nanopolish_sitelevel",
        "swarm_nanopolish_scorebed",
    ]:
        lines.append(f"{key}: {getattr(args, key)}")

    (outdir / "README.txt").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Run SWARM threshold scans, SWARM/m6Anet complementarity analysis, "
            "and optional f5c-vs-nanopolish SWARM comparison."
        )
    )

    required = parser.add_argument_group("required input/output paths")
    required.add_argument("--common-bed", required=True, help="BED file defining the common evaluation universe.")
    required.add_argument("--positive-bed", required=True, help="BED file defining reference-positive sites.")
    required.add_argument("--exclude-bed", required=True, help="BED file defining sites to exclude from pseudo-negative construction.")
    required.add_argument("--swarm-scorebed", required=True, help="SWARM score BED file. The 4th column is interpreted as SWARM probability/score.")
    required.add_argument("--m6anet-scorebed", required=True, help="m6Anet score BED file. The 4th column is interpreted as m6Anet probability_modified/score.")
    required.add_argument("--outdir", required=True, help="Output directory for tables and figures.")

    optional_paths = parser.add_argument_group("optional input paths")
    optional_paths.add_argument("--swarm-stoich-table", default=None, help="Optional SWARM site-level/detail table containing stoichiometry or modification ratio.")
    optional_paths.add_argument("--swarm-nanopolish-sitelevel", default=None, help="Optional nanopolish-SWARM site-level table.")
    optional_paths.add_argument("--swarm-nanopolish-scorebed", default=None, help="Optional nanopolish-SWARM score BED file. Used only if --swarm-nanopolish-sitelevel is not provided.")

    settings = parser.add_argument_group("analysis settings")
    settings.add_argument("--exclude-buffer", type=int, default=2, help="Coordinate buffer used to remove ambiguous pseudo-negative sites near excluded sites. Default: 2.")
    settings.add_argument("--swarm-prob-cutoff", type=float, default=0.9972, help="Fixed SWARM probability cutoff. Default: 0.9972.")
    settings.add_argument("--swarm-stoich-cutoff", type=float, default=0.1, help="Fixed SWARM stoichiometry cutoff. Default: 0.1.")
    settings.add_argument("--m6anet-cutoff", type=float, default=0.9, help="Fixed m6Anet probability cutoff. Default: 0.9.")

    return parser.parse_args()


def main() -> None:
    """Run the complete analysis workflow."""
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for name in ["common_bed", "positive_bed", "exclude_bed", "swarm_scorebed", "m6anet_scorebed"]:
        check(getattr(args, name), name, required=True)

    for name in ["swarm_stoich_table", "swarm_nanopolish_sitelevel", "swarm_nanopolish_scorebed"]:
        check(getattr(args, name), name, required=False)

    log("Building labeled common-site table...")
    df = build_table(args)
    df.to_csv(outdir / "base_eval_table.with_SWARM_m6Anet.tsv", sep="\t", index=False)

    log(f"Evaluation sites: {len(df)}")
    log(f"Reference positives: {int((df['label'] == 1).sum())}")
    log(f"Pseudo-negatives: {int((df['label'] == 0).sum())}")

    write_run_summary(df, outdir, args)

    log("Running SWARM threshold sensitivity analysis...")
    run_threshold(df, outdir, args)

    log("Running SWARM-m6Anet complementarity analysis...")
    run_complementarity(df, outdir, args)

    log("Running optional f5c-vs-nanopolish SWARM analysis...")
    run_nanopolish(df, outdir, args)

    log("Done.")
    log(f"Output directory: {outdir}")


if __name__ == "__main__":
    main()
