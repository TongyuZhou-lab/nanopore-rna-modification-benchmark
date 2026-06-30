#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot cross-species m6A ROC/PR curves with a unified style.

This command-line version removes hard-coded absolute paths from the original script.
It supports:
  1) mouse: read existing curve-point and metric tables
  2) human: recompute curves from common-universe BED/table + reference BEDs + scorebeds
  3) Arabidopsis: read existing curve-point and metric tables, or recompute from a labeled input table

Example:
python plot_m6a_cross_species_rocp.py \
  --outroot /path/to/output \
  --mouse-curve /path/to/ROC_PR_curve_points_common.tsv \
  --mouse-metrics /path/to/ROC_PR_metrics_common.tsv \
  --arab-curve /path/to/ROC_PR_curve_points.primary_nanopore.TandemMod_p95rate.tsv \
  --arab-metrics /path/to/ROC_PR_metrics.primary_nanopore.TandemMod_p95rate.tsv \
  --human-common-universe /path/to/human_four_tool_common_universe.site_ids.tsv \
  --human-truth-bed /path/to/strict_consensus_2of3.bed \
  --human-loose-bed /path/to/loose_union_3way.bed \
  --human-tandemmod-scorebed /path/to/TandemMod.score.bed \
  --human-m6anet-scorebed /path/to/m6Anet.score.bed \
  --human-cheui-scorebed /path/to/CHEUI.score.bed \
  --human-swarm-scorebed /path/to/SWARM.score.bed
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.metrics import (
    roc_curve,
    precision_recall_curve,
    roc_auc_score,
    average_precision_score,
)


TOOL_ORDER = ["TandemMod", "m6Anet", "CHEUI", "SWARM"]
TOOL_COLORS = {
    "TandemMod": "#1f77b4",
    "m6Anet": "#ff7f0e",
    "CHEUI": "#2ca02c",
    "SWARM": "#d62728",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot cross-species m6A ROC/PR curves without hard-coded paths."
    )

    # Output
    parser.add_argument("--outroot", required=True, type=Path, help="Output root directory.")

    # Mouse: existing curve and metrics
    parser.add_argument("--mouse-curve", required=True, type=Path, help="Mouse ROC/PR curve-point TSV.")
    parser.add_argument("--mouse-metrics", required=True, type=Path, help="Mouse ROC/PR metrics TSV.")

    # Arabidopsis: existing curves, or recompute from labeled table
    parser.add_argument("--arab-curve", type=Path, default=None, help="Arabidopsis ROC/PR curve-point TSV.")
    parser.add_argument("--arab-metrics", type=Path, default=None, help="Arabidopsis ROC/PR metrics TSV.")
    parser.add_argument(
        "--arab-input",
        type=Path,
        default=None,
        help=(
            "Optional Arabidopsis labeled input table for recomputing curves. "
            "Required only when --arab-curve/--arab-metrics are not provided."
        ),
    )
    parser.add_argument("--arab-truth-name", default="primary_nanopore", help="Arabidopsis truth-set name.")
    parser.add_argument("--arab-label-col", default="label", help="Arabidopsis label column name.")
    parser.add_argument("--arab-tandemmod-col", default="mod_rate_p95", help="Arabidopsis TandemMod score column.")
    parser.add_argument("--arab-m6anet-col", default="m6Anet_txShift0", help="Arabidopsis m6Anet score column.")
    parser.add_argument("--arab-cheui-col", default="CHEUI_origin1_shiftMinus5", help="Arabidopsis CHEUI score column.")
    parser.add_argument("--arab-swarm-col", default="SWARM_txCenterPlus5", help="Arabidopsis SWARM score column.")

    # Human: recompute from common universe + reference BEDs + scorebeds
    parser.add_argument("--human-common-universe", required=True, type=Path, help="Human common-universe table/BED.")
    parser.add_argument("--human-truth-bed", required=True, type=Path, help="Human positive reference BED.")
    parser.add_argument("--human-loose-bed", required=True, type=Path, help="Human loose-union BED for ambiguous-site exclusion.")
    parser.add_argument("--human-truth-name", default="strict_consensus_2of3", help="Human truth-set name.")
    parser.add_argument("--human-buffer-bp", type=int, default=2, help="Buffer size around loose-union sites for exclusion.")
    parser.add_argument("--human-tandemmod-scorebed", required=True, type=Path, help="Human TandemMod score BED.")
    parser.add_argument("--human-m6anet-scorebed", required=True, type=Path, help="Human m6Anet score BED.")
    parser.add_argument("--human-cheui-scorebed", required=True, type=Path, help="Human CHEUI score BED.")
    parser.add_argument("--human-swarm-scorebed", required=True, type=Path, help="Human SWARM score BED.")
    parser.add_argument(
        "--human-cache-dir",
        type=Path,
        default=None,
        help="Directory to save recomputed human curve/metric/eval tables. Default: OUTROOT/tables/human_cache",
    )

    # Plot style
    parser.add_argument("--fig-width", type=float, default=4.0, help="Figure width in inches.")
    parser.add_argument("--fig-height", type=float, default=4.3, help="Figure height in inches.")
    parser.add_argument("--line-width", type=float, default=1.8, help="Main curve line width.")
    parser.add_argument("--ref-line-width", type=float, default=0.9, help="ROC diagonal / PR baseline line width.")
    parser.add_argument("--axis-line-width", type=float, default=0.8, help="Axis spine line width.")
    parser.add_argument("--tick-line-width", type=float, default=0.8, help="Tick line width.")
    parser.add_argument("--tick-length", type=float, default=3.0, help="Tick length.")
    parser.add_argument("--tick-labelsize", type=float, default=12, help="Tick-number font size.")
    parser.add_argument("--axis-labelsize", type=float, default=9, help="Axis-label font size for full version.")
    parser.add_argument("--title-size", type=float, default=10, help="Title font size for full version.")
    parser.add_argument("--legend-size", type=float, default=7.5, help="Legend font size for full version.")
    parser.add_argument("--dpi", type=int, default=300, help="PNG output DPI.")
    parser.add_argument(
        "--tick-values",
        default="0,0.25,0.5,0.75,1.0",
        help="Comma-separated tick values, default: 0,0.25,0.5,0.75,1.0",
    )
    parser.add_argument(
        "--keep-full-box",
        action="store_true",
        help="Keep all four spines. By default, top/right spines are removed.",
    )
    parser.add_argument(
        "--clean-keep-axis-labels",
        action="store_true",
        help="Keep x/y axis labels in clean figures. Default: clean figures hide axis labels.",
    )

    return parser.parse_args()


def check_file(path: Path | None, label: str) -> None:
    if path is None:
        raise ValueError(f"{label} is required.")
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def canon_tool_name(x):
    s = str(x).strip()
    low = s.lower()
    if "tandemmod" in low:
        return "TandemMod"
    if "m6anet" in low:
        return "m6Anet"
    if "cheui" in low:
        return "CHEUI"
    if "swarm" in low:
        return "SWARM"
    return s


def safe_div(a, b):
    return np.nan if b == 0 else a / b


def read_table_auto(path):
    path = Path(path)
    try:
        df = pd.read_csv(path, sep="\t")
        if df.shape[1] > 1:
            return df
    except Exception:
        pass
    return pd.read_csv(path)


def normalize_curve_df(df):
    out = df.copy()
    out.columns = [str(c).strip() for c in out.columns]
    lower = {c.lower(): c for c in out.columns}

    rename = {}
    if "tool" not in out.columns and "tool" in lower:
        rename[lower["tool"]] = "tool"
    if "curve" not in out.columns and "curve" in lower:
        rename[lower["curve"]] = "curve"

    if "x" not in out.columns:
        if "fpr" in lower:
            rename[lower["fpr"]] = "x"
        elif "recall" in lower:
            rename[lower["recall"]] = "x"
    if "y" not in out.columns:
        if "tpr" in lower:
            rename[lower["tpr"]] = "y"
        elif "precision" in lower:
            rename[lower["precision"]] = "y"

    if rename:
        out = out.rename(columns=rename)

    required = ["tool", "curve", "x", "y"]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise ValueError(f"Curve table missing columns: {missing}")

    out["tool"] = out["tool"].map(canon_tool_name)
    out["curve"] = out["curve"].astype(str).str.upper().replace({"PRC": "PR", "PRECISION_RECALL": "PR"})
    out["x"] = pd.to_numeric(out["x"], errors="coerce")
    out["y"] = pd.to_numeric(out["y"], errors="coerce")
    out = out.dropna(subset=["tool", "curve", "x", "y"]).copy()
    return out


def normalize_metrics_df(df):
    out = df.copy()
    out.columns = [str(c).strip() for c in out.columns]
    lower = {c.lower(): c for c in out.columns}

    rename = {}
    if "tool" not in out.columns and "tool" in lower:
        rename[lower["tool"]] = "tool"
    for lc, raw in lower.items():
        if lc == "auroc":
            rename[raw] = "AUROC"
        elif lc in {"auprc", "aupr", "auprc_average_precision", "average_precision"}:
            rename[raw] = "AUPRC"
        elif lc == "positive_prevalence":
            rename[raw] = "positive_prevalence"
    if rename:
        out = out.rename(columns=rename)

    if "tool" not in out.columns:
        raise ValueError("Metrics table missing tool column")
    out["tool"] = out["tool"].map(canon_tool_name)

    out["AUROC"] = pd.to_numeric(out["AUROC"], errors="coerce") if "AUROC" in out.columns else np.nan
    out["AUPRC"] = pd.to_numeric(out["AUPRC"], errors="coerce") if "AUPRC" in out.columns else np.nan
    out["positive_prevalence"] = (
        pd.to_numeric(out["positive_prevalence"], errors="coerce")
        if "positive_prevalence" in out.columns
        else np.nan
    )
    return out


def read_bed_sites(path):
    sites = set()
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fs = line.rstrip("\n").split("\t")
            if len(fs) < 3:
                continue
            try:
                sites.add((fs[0], int(fs[1]), int(fs[2])))
            except Exception:
                continue
    return sites


def read_common_universe(path):
    path = Path(path)
    try:
        df = pd.read_csv(path, sep="\t")
        if {"chrom", "start", "end"}.issubset(df.columns):
            return set(zip(df["chrom"], df["start"].astype(int), df["end"].astype(int)))
    except Exception:
        pass

    sites = set()
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fs = line.rstrip("\n").split("\t")
            if len(fs) >= 3:
                try:
                    sites.add((fs[0], int(fs[1]), int(fs[2])))
                except Exception:
                    continue
    return sites


def read_scorebed(path):
    score = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fs = line.rstrip("\n").split("\t")
            if len(fs) < 4:
                continue
            try:
                site = (fs[0], int(fs[1]), int(fs[2]))
                val = float(fs[3])
            except Exception:
                continue
            if site not in score or val > score[site]:
                score[site] = val
    return score


def make_buffered_sites(sites, bp):
    out = set()
    for chrom, start, end in sites:
        for d in range(-bp, bp + 1):
            ns = start + d
            ne = end + d
            if ns < 0:
                continue
            out.add((chrom, ns, ne))
    return out


def calc_metrics_and_curves_from_scores(y, score_by_tool, score_name_by_tool=None, fillna_zero=False):
    y = np.asarray(y).astype(int)
    n = len(y)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    prevalence = safe_div(n_pos, n)

    metric_rows = []
    curve_rows = []

    for raw_tool, raw_scores in score_by_tool.items():
        tool = canon_tool_name(raw_tool)
        scores = pd.to_numeric(pd.Series(raw_scores), errors="coerce")
        if fillna_zero:
            scores = scores.fillna(0)

        mask = ~scores.isna().to_numpy()
        yy = y[mask]
        ss = scores.to_numpy(dtype=float)[mask]

        if len(np.unique(yy)) < 2:
            auroc = np.nan
            auprc = np.nan
        else:
            fpr, tpr, _ = roc_curve(yy, ss)
            precision, recall, _ = precision_recall_curve(yy, ss)
            auroc = roc_auc_score(yy, ss)
            auprc = average_precision_score(yy, ss)

            for x, val in zip(fpr, tpr):
                curve_rows.append({"tool": tool, "curve": "ROC", "x": x, "y": val})

            order = np.argsort(recall)
            for x, val in zip(recall[order], precision[order]):
                curve_rows.append({"tool": tool, "curve": "PR", "x": x, "y": val})

        score_name = "" if score_name_by_tool is None else score_name_by_tool.get(raw_tool, "")
        metric_rows.append({
            "tool": tool,
            "raw_tool": raw_tool,
            "score_name": score_name,
            "n_kept_sites": int(mask.sum()),
            "n_positive": n_pos,
            "n_pseudo_negative": n_neg,
            "positive_prevalence": prevalence,
            "AUROC": auroc,
            "AUPRC": auprc,
        })

    return pd.DataFrame(metric_rows), pd.DataFrame(curve_rows)


def load_existing_curves_and_metrics(curve_path, metric_path):
    check_file(curve_path, "curve table")
    check_file(metric_path, "metrics table")
    curves = normalize_curve_df(read_table_auto(curve_path))
    metrics = normalize_metrics_df(read_table_auto(metric_path))
    return curves, metrics, str(curve_path), str(metric_path)


def compute_arabidopsis(args, table_dir: Path):
    if args.arab_curve is not None or args.arab_metrics is not None:
        if args.arab_curve is None or args.arab_metrics is None:
            raise ValueError("Please provide both --arab-curve and --arab-metrics, or provide neither and use --arab-input.")
        return load_existing_curves_and_metrics(args.arab_curve, args.arab_metrics)

    check_file(args.arab_input, "Arabidopsis labeled input table")
    df = pd.read_csv(args.arab_input, sep="\t")

    tool_to_col = {
        "TandemMod_p95rate": args.arab_tandemmod_col,
        "m6Anet_txShift0": args.arab_m6anet_col,
        "CHEUI_origin1_shiftMinus5": args.arab_cheui_col,
        "SWARM_txCenterPlus5": args.arab_swarm_col,
    }

    required = [args.arab_label_col] + list(tool_to_col.values())
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Arabidopsis input missing columns: {missing}")

    y = df[args.arab_label_col].astype(int).to_numpy()
    score_by_tool = {
        tool: pd.to_numeric(df[col], errors="coerce").fillna(0).to_numpy()
        for tool, col in tool_to_col.items()
    }

    metrics, curves = calc_metrics_and_curves_from_scores(y, score_by_tool, tool_to_col, fillna_zero=True)
    metrics.insert(0, "truth_set", args.arab_truth_name)
    curves.insert(0, "truth_set", args.arab_truth_name)

    curve_path = table_dir / f"arabidopsis_ROC_PR_curve_points.{args.arab_truth_name}.tsv"
    metric_path = table_dir / f"arabidopsis_ROC_PR_metrics.{args.arab_truth_name}.tsv"
    curves.to_csv(curve_path, sep="\t", index=False)
    metrics.to_csv(metric_path, sep="\t", index=False)

    return normalize_curve_df(curves), normalize_metrics_df(metrics), str(curve_path), str(metric_path)


def compute_human(args, table_dir: Path):
    required_paths = [
        args.human_common_universe,
        args.human_truth_bed,
        args.human_loose_bed,
        args.human_tandemmod_scorebed,
        args.human_m6anet_scorebed,
        args.human_cheui_scorebed,
        args.human_swarm_scorebed,
    ]
    for path in required_paths:
        check_file(path, "human input file")

    cache_dir = args.human_cache_dir if args.human_cache_dir is not None else table_dir / "human_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    human_tools = {
        "TandemMod": {
            "score_bed": args.human_tandemmod_scorebed,
            "score_name": "max_ratio = max(p_0.95/total)",
        },
        "m6Anet": {
            "score_bed": args.human_m6anet_scorebed,
            "score_name": "probability_modified",
        },
        "CHEUI": {
            "score_bed": args.human_cheui_scorebed,
            "score_name": "probability",
        },
        "SWARM": {
            "score_bed": args.human_swarm_scorebed,
            "score_name": "probability",
        },
    }

    common_universe = read_common_universe(args.human_common_universe)
    truth_sites = read_bed_sites(args.human_truth_bed)
    loose_sites = read_bed_sites(args.human_loose_bed)
    loose_buffered = make_buffered_sites(loose_sites, args.human_buffer_bp)

    positives = common_universe & truth_sites
    pseudo_negatives = common_universe - loose_buffered
    kept = sorted(positives | pseudo_negatives, key=lambda x: (x[0], x[1], x[2]))
    y = np.array([1 if site in positives else 0 for site in kept], dtype=int)

    eval_df = pd.DataFrame({
        "chrom": [site[0] for site in kept],
        "start": [site[1] for site in kept],
        "end": [site[2] for site in kept],
        "site_id": [f"{site[0]}:{site[1]}-{site[2]}" for site in kept],
        "label": y,
    })

    score_by_tool = {}
    score_name_by_tool = {}
    for tool, info in human_tools.items():
        smap = read_scorebed(info["score_bed"])
        vals = [smap.get(site, np.nan) for site in kept]
        eval_df[tool] = vals
        score_by_tool[tool] = vals
        score_name_by_tool[tool] = info["score_name"]

    metrics, curves = calc_metrics_and_curves_from_scores(
        y,
        score_by_tool,
        score_name_by_tool,
        fillna_zero=False,
    )

    metrics.insert(0, "truth_set", args.human_truth_name)
    curves.insert(0, "truth_set", args.human_truth_name)

    eval_path = cache_dir / f"human_ROCPR_eval_table.{args.human_truth_name}.tsv"
    metric_path = cache_dir / f"human_ROCPR_metrics.{args.human_truth_name}.tsv"
    curve_path = cache_dir / f"human_ROC_PR_curve_points.{args.human_truth_name}.tsv"

    eval_df.to_csv(eval_path, sep="\t", index=False)
    metrics.to_csv(metric_path, sep="\t", index=False)
    curves.to_csv(curve_path, sep="\t", index=False)

    return normalize_curve_df(curves), normalize_metrics_df(metrics), str(curve_path), str(metric_path)


def get_prevalence(metrics):
    if "positive_prevalence" not in metrics.columns:
        return np.nan
    vals = pd.to_numeric(metrics["positive_prevalence"], errors="coerce").dropna()
    if len(vals) == 0:
        return np.nan
    return float(vals.iloc[0])


def metric_dict(metrics):
    out = {}
    for _, row in metrics.iterrows():
        tool = canon_tool_name(row["tool"])
        out[tool] = {
            "AUROC": row.get("AUROC", np.nan),
            "AUPRC": row.get("AUPRC", np.nan),
        }
    return out


def configure_matplotlib(args):
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": args.title_size,
        "axes.labelsize": args.axis_labelsize,
        "xtick.labelsize": args.tick_labelsize,
        "ytick.labelsize": args.tick_labelsize,
        "legend.fontsize": args.legend_size,
        "axes.linewidth": args.axis_line_width,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })


def style_ax(ax, args, xlabel, ylabel, title=None, full=True, roc=False):
    tick_values = [float(x) for x in str(args.tick_values).split(",") if str(x).strip() != ""]

    if full and title:
        ax.set_title(title, pad=5)
    else:
        ax.set_title("")

    if full or args.clean_keep_axis_labels:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    else:
        ax.set_xlabel("")
        ax.set_ylabel("")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks(tick_values)
    ax.set_yticks(tick_values)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    if args.keep_full_box:
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(args.axis_line_width)
        top = True
        right = True
    else:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(True)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_linewidth(args.axis_line_width)
        ax.spines["bottom"].set_linewidth(args.axis_line_width)
        top = False
        right = False

    ax.tick_params(
        axis="both",
        which="both",
        width=args.tick_line_width,
        length=args.tick_length,
        top=top,
        right=right,
        bottom=True,
        left=True,
        labelbottom=True,
        labelleft=True,
    )

    if roc:
        ax.plot([0, 1], [0, 1], linestyle="--", linewidth=args.ref_line_width, color="0.6", zorder=0)


def save_all(fig, outbase, dpi):
    outbase = Path(outbase)
    fig.savefig(str(outbase) + ".svg", bbox_inches="tight", transparent=True)
    fig.savefig(str(outbase) + ".png", dpi=dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)


def plot_one(args, dirs, species_key, species_label, curves, metrics, curve_type, full=True):
    full_dir, clean_dir, _ = dirs
    md = metric_dict(metrics)
    prevalence = get_prevalence(metrics)

    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))

    for tool in TOOL_ORDER:
        sub = curves[(curves["tool"] == tool) & (curves["curve"] == curve_type)].copy()
        if sub.empty:
            continue
        sub = sub.sort_values("x")

        if full:
            if curve_type == "ROC":
                value = md.get(tool, {}).get("AUROC", np.nan)
                label = f"{tool} (AUROC={value:.3f})" if pd.notna(value) else tool
            else:
                value = md.get(tool, {}).get("AUPRC", np.nan)
                label = f"{tool} (AUPRC={value:.3f})" if pd.notna(value) else tool
        else:
            label = None

        ax.plot(
            sub["x"].to_numpy(),
            sub["y"].to_numpy(),
            linewidth=args.line_width,
            color=TOOL_COLORS[tool],
            label=label,
            zorder=3,
        )

    if curve_type == "ROC":
        style_ax(args=args, ax=ax, xlabel="False Positive Rate", ylabel="True Positive Rate",
                 title=f"{species_label} m6A ROC", full=full, roc=True)
        if full:
            ax.legend(frameon=False, loc="lower right")
    else:
        if pd.notna(prevalence):
            ax.axhline(prevalence, linestyle="--", linewidth=args.ref_line_width, color="0.6", zorder=0)
        style_ax(args=args, ax=ax, xlabel="Recall", ylabel="Precision",
                 title=f"{species_label} m6A PR", full=full, roc=False)
        if full:
            ax.legend(frameon=False, loc="lower left")

    outdir = full_dir if full else clean_dir
    suffix = "full" if full else "clean"
    save_all(fig, outdir / f"{species_key}_{curve_type}.{suffix}", dpi=args.dpi)


def main():
    args = parse_args()
    configure_matplotlib(args)

    full_dir = args.outroot / "full_version"
    clean_dir = args.outroot / "clean_version"
    table_dir = args.outroot / "tables"
    for directory in [full_dir, clean_dir, table_dir]:
        directory.mkdir(parents=True, exist_ok=True)
    dirs = (full_dir, clean_dir, table_dir)

    all_metrics = []
    source_rows = []
    datasets = []

    print("\n=== Loading Mouse ===")
    curves, metrics, curve_path, metric_path = load_existing_curves_and_metrics(args.mouse_curve, args.mouse_metrics)
    datasets.append(("mouse", "Mouse", curves, metrics))
    source_rows.append({"species": "Mouse", "curve_file": curve_path, "metrics_file": metric_path})
    print("Curve:", curve_path)
    print("Metrics:", metric_path)

    print("\n=== Computing Human ===")
    curves, metrics, curve_path, metric_path = compute_human(args, table_dir)
    datasets.append(("human", "Human", curves, metrics))
    source_rows.append({"species": "Human", "curve_file": curve_path, "metrics_file": metric_path})
    print("Curve:", curve_path)
    print("Metrics:", metric_path)

    print("\n=== Loading / Computing Arabidopsis ===")
    curves, metrics, curve_path, metric_path = compute_arabidopsis(args, table_dir)
    datasets.append(("arabidopsis", "Arabidopsis", curves, metrics))
    source_rows.append({"species": "Arabidopsis", "curve_file": curve_path, "metrics_file": metric_path})
    print("Curve:", curve_path)
    print("Metrics:", metric_path)

    for species_key, species_label, curves, metrics in datasets:
        curves.to_csv(table_dir / f"{species_key}.curve_points.used_for_plot.tsv", sep="\t", index=False)
        metrics.to_csv(table_dir / f"{species_key}.metrics.used_for_plot.tsv", sep="\t", index=False)

        metrics2 = metrics.copy()
        metrics2.insert(0, "species", species_label)
        all_metrics.append(metrics2)

        for curve_type in ["ROC", "PR"]:
            plot_one(args, dirs, species_key, species_label, curves, metrics, curve_type, full=True)
            plot_one(args, dirs, species_key, species_label, curves, metrics, curve_type, full=False)

    pd.concat(all_metrics, ignore_index=True).to_csv(
        table_dir / "cross_species_m6A_ROCPR_metrics.used_for_plot.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(source_rows).to_csv(
        table_dir / "cross_species_m6A_ROCPR_source_files.tsv",
        sep="\t",
        index=False,
    )

    print("\n===== Done =====")
    print("Output root:", args.outroot)
    print("Full version:", full_dir)
    print("Clean version:", clean_dir)
    print("Tables:", table_dir)


if __name__ == "__main__":
    main()
