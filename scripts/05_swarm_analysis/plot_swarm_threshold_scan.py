#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot SWARM threshold-scan figures with full and clean versions.

This command-line version removes hard-coded paths from the original script.
It reads four input tables:
  1. SWARM probability-threshold scan
  2. SWARM stoichiometry-threshold scan
  3. SWARM joint probability/stoichiometry threshold scan
  4. SWARM/m6Anet strategy metrics

Main outputs:
  OUTDIR/full/
  OUTDIR/clean/
  OUTDIR/tables/

Example:
python plot_swarm_threshold_scan.py \
  --indir /path/to/swarm_m6anet_deep_analysis_mouse \
  --outdir /path/to/output_figures

Or specify files separately:
python plot_swarm_threshold_scan.py \
  --prob-tsv /path/to/SWARM_probability_threshold_scan.tsv \
  --stoich-tsv /path/to/SWARM_stoichiometry_threshold_scan.tsv \
  --joint-tsv /path/to/SWARM_joint_probability_stoichiometry_threshold_scan.tsv \
  --strategy-tsv /path/to/SWARM_m6Anet_strategy_metrics.tsv \
  --outdir /path/to/output_figures
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot SWARM threshold-scan figures without hard-coded paths."
    )

    parser.add_argument(
        "--indir",
        type=Path,
        default=None,
        help=(
            "Input directory containing the four default TSV files. "
            "If provided, missing individual TSV paths will be inferred from this directory."
        ),
    )
    parser.add_argument("--prob-tsv", type=Path, default=None, help="SWARM probability threshold scan TSV.")
    parser.add_argument("--stoich-tsv", type=Path, default=None, help="SWARM stoichiometry threshold scan TSV.")
    parser.add_argument("--joint-tsv", type=Path, default=None, help="SWARM joint probability/stoichiometry threshold scan TSV.")
    parser.add_argument("--strategy-tsv", type=Path, default=None, help="SWARM/m6Anet strategy metrics TSV.")

    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument("--prob-official", type=float, default=0.9972, help="Official SWARM probability cutoff.")
    parser.add_argument("--stoich-official", type=float, default=0.1, help="Official SWARM stoichiometry cutoff.")
    parser.add_argument("--heatmap-value", default="F1", help="Metric to show in the joint-threshold heatmap.")

    # Unified style
    parser.add_argument("--fig-width", type=float, default=3.0, help="Single-figure width.")
    parser.add_argument("--fig-height", type=float, default=3.0, help="Single-figure height.")
    parser.add_argument("--heatmap-width", type=float, default=3.25, help="Heatmap figure width.")
    parser.add_argument("--heatmap-height", type=float, default=3.0, help="Heatmap figure height.")
    parser.add_argument("--line-width", type=float, default=1.8, help="Main line width.")
    parser.add_argument("--ref-line-width", type=float, default=0.9, help="Reference line width.")
    parser.add_argument("--axis-line-width", type=float, default=0.8, help="Axis spine width.")
    parser.add_argument("--tick-line-width", type=float, default=0.8, help="Tick line width.")
    parser.add_argument("--tick-length", type=float, default=3.0, help="Tick length.")
    parser.add_argument("--tick-labelsize", type=float, default=12, help="Tick-number font size.")
    parser.add_argument("--axis-labelsize", type=float, default=9, help="Axis-label font size for full figures.")
    parser.add_argument("--title-size", type=float, default=10, help="Title font size for full figures.")
    parser.add_argument("--legend-size", type=float, default=7.0, help="Legend font size for full figures.")
    parser.add_argument("--dpi", type=int, default=300, help="PNG output DPI.")
    parser.add_argument(
        "--tick-values",
        default="0,0.25,0.5,0.75,1.0",
        help="Comma-separated tick values for 0-1 metric axes.",
    )
    parser.add_argument(
        "--keep-full-box",
        action="store_true",
        help="Keep all four spines. By default, top and right spines are removed.",
    )
    parser.add_argument(
        "--clean-keep-axis-labels",
        action="store_true",
        help="Keep x/y axis labels in clean figures. Default: clean figures hide axis labels.",
    )

    return parser.parse_args()


def configure_paths(args: argparse.Namespace) -> None:
    if args.indir is not None:
        if args.prob_tsv is None:
            args.prob_tsv = args.indir / "SWARM_probability_threshold_scan.tsv"
        if args.stoich_tsv is None:
            args.stoich_tsv = args.indir / "SWARM_stoichiometry_threshold_scan.tsv"
        if args.joint_tsv is None:
            args.joint_tsv = args.indir / "SWARM_joint_probability_stoichiometry_threshold_scan.tsv"
        if args.strategy_tsv is None:
            args.strategy_tsv = args.indir / "SWARM_m6Anet_strategy_metrics.tsv"

    required = {
        "--prob-tsv": args.prob_tsv,
        "--stoich-tsv": args.stoich_tsv,
        "--joint-tsv": args.joint_tsv,
        "--strategy-tsv": args.strategy_tsv,
    }
    missing_args = [name for name, path in required.items() if path is None]
    if missing_args:
        raise SystemExit(
            "Missing input files. Provide --indir or specify these arguments: "
            + ", ".join(missing_args)
        )

    for name, path in required.items():
        if not path.exists():
            raise FileNotFoundError(f"{name} not found: {path}")


def configure_matplotlib(args: argparse.Namespace) -> None:
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


def make_output_dirs(args: argparse.Namespace):
    full_dir = args.outdir / "full"
    clean_dir = args.outdir / "clean"
    table_dir = args.outdir / "tables"
    for directory in [full_dir, clean_dir, table_dir]:
        directory.mkdir(parents=True, exist_ok=True)
    return full_dir, clean_dir, table_dir


def tick_values(args: argparse.Namespace):
    return [float(x) for x in str(args.tick_values).split(",") if str(x).strip() != ""]


def save_all(fig, args: argparse.Namespace, outdirs, base, full=True):
    full_dir, clean_dir, _ = outdirs
    outdir = full_dir if full else clean_dir
    suffix = "full" if full else "clean"
    fig.savefig(outdir / f"{base}.{suffix}.svg", bbox_inches="tight", transparent=True)
    fig.savefig(outdir / f"{base}.{suffix}.pdf", bbox_inches="tight", transparent=True)
    fig.savefig(outdir / f"{base}.{suffix}.png", dpi=args.dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)


def norm_col(c):
    return str(c).strip().lower().replace(" ", "_").replace("-", "_").replace(".", "_")


def find_col(df, candidates, contains_any=None):
    cols = list(df.columns)
    norm_map = {norm_col(c): c for c in cols}

    for cand in candidates:
        key = norm_col(cand)
        if key in norm_map:
            return norm_map[key]

    if contains_any:
        for c in cols:
            nc = norm_col(c)
            if all(x.lower() in nc for x in contains_any):
                return c

    raise KeyError(f"Cannot find column. candidates={candidates}, columns={cols}")


def maybe_col(df, candidates, contains_any=None):
    try:
        return find_col(df, candidates, contains_any)
    except Exception:
        return None


def get_metric_cols(df):
    precision_col = find_col(df, ["Precision", "precision"])
    recall_col = find_col(df, ["Recall", "recall"])
    f1_col = find_col(df, ["F1", "f1", "F1_score", "f1_score"])

    pred_col = maybe_col(
        df,
        [
            "predicted_positive",
            "predicted_positive_sites",
            "pred_positive",
            "pred_positive_sites",
            "n_predicted_positive",
            "positive_calls",
            "n_positive_calls",
            "predicted_sites",
            "n_predicted_sites",
            "n_predicted",
            "call_count",
            "n_calls",
        ],
        contains_any=["pred"],
    )

    return precision_col, recall_col, f1_col, pred_col


def apply_spine_style(ax, args: argparse.Namespace):
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
        labelsize=args.tick_labelsize,
        labelbottom=True,
        labelleft=True,
        bottom=True,
        left=True,
        top=top,
        right=right,
        direction="out",
    )


def style_left_metric_ax(ax, args: argparse.Namespace, xlabel, ylabel="Metric", title=None, full=True):
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

    tv = tick_values(args)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks(tv)
    ax.set_yticks(tv)
    ax.set_xticklabels([f"{x:g}" if x in {0, 1} else str(x) for x in tv])
    ax.set_yticklabels([f"{x:g}" if x in {0, 1} else str(x) for x in tv])

    ax.grid(False)
    apply_spine_style(ax, args)


def style_general_ax(ax, args: argparse.Namespace):
    ax.grid(False)
    apply_spine_style(ax, args)


def nearest_row(df, x_col, x):
    tmp = df.dropna(subset=[x_col]).copy()
    tmp["_dist"] = (tmp[x_col] - x).abs()
    return tmp.sort_values("_dist").iloc[0]


def plot_threshold_scan_with_predicted(args, outdirs, tsv, x_candidates, title, xlabel, outbase, official_x):
    df = pd.read_csv(tsv, sep="\t")

    x_col = None
    for cand in x_candidates:
        if cand in df.columns:
            x_col = cand
            break

    if x_col is None:
        threshold_cols = [c for c in df.columns if "threshold" in norm_col(c)]
        if not threshold_cols:
            raise KeyError(f"Cannot find threshold column in {tsv}: {df.columns.tolist()}")
        x_col = threshold_cols[0]

    precision_col, recall_col, f1_col, pred_col = get_metric_cols(df)

    for c in [x_col, precision_col, recall_col, f1_col]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    if pred_col is not None:
        df[pred_col] = pd.to_numeric(df[pred_col], errors="coerce")

    df = df.dropna(subset=[x_col, precision_col, recall_col, f1_col]).copy()
    df = df.sort_values(x_col)

    best = df.dropna(subset=[f1_col]).sort_values(
        [f1_col, precision_col, recall_col],
        ascending=[False, False, False],
    ).iloc[0]

    best_x = float(best[x_col])
    best_f1 = float(best[f1_col])
    best_p = float(best[precision_col])
    best_r = float(best[recall_col])

    official = nearest_row(df, x_col, official_x)
    off_f1 = float(official[f1_col])
    off_p = float(official[precision_col])
    off_r = float(official[recall_col])

    used_cols = [x_col, precision_col, recall_col, f1_col]
    if pred_col is not None:
        used_cols.append(pred_col)
    df[used_cols].to_csv(outdirs[2] / f"{outbase}.used_table.tsv", sep="\t", index=False)

    for full in [True, False]:
        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))

        l1, = ax.plot(df[x_col], df[precision_col], lw=args.line_width, label="Precision")
        l2, = ax.plot(df[x_col], df[recall_col], lw=args.line_width, label="Recall")
        l3, = ax.plot(df[x_col], df[f1_col], lw=args.line_width, label="F1")

        ax.axvline(official_x, ls="--", lw=args.ref_line_width, color="0.55")
        ax.axvline(best_x, ls=":", lw=args.ref_line_width, color="0.35")

        style_left_metric_ax(ax, args, xlabel=xlabel, title=title, full=full)

        ax2 = None
        l4 = None
        if pred_col is not None:
            ax2 = ax.twinx()
            l4, = ax2.plot(
                df[x_col],
                df[pred_col],
                lw=args.ref_line_width,
                ls="--",
                color=l1.get_color(),
                alpha=0.95,
                label="Predicted positive",
            )

            if full:
                ax2.set_ylabel("Predicted positive sites")
            else:
                ax2.set_ylabel("")

            if args.keep_full_box:
                for spine in ax2.spines.values():
                    spine.set_visible(True)
                    spine.set_linewidth(args.axis_line_width)
                ax2.tick_params(
                    width=args.tick_line_width,
                    length=args.tick_length,
                    labelsize=args.tick_labelsize,
                )
            else:
                ax2.spines["top"].set_visible(False)
                ax2.spines["left"].set_visible(False)
                ax2.spines["bottom"].set_visible(False)
                ax2.spines["right"].set_visible(False)
                ax2.tick_params(
                    width=args.tick_line_width,
                    length=args.tick_length,
                    labelsize=args.tick_labelsize,
                    top=False,
                    right=False,
                    left=False,
                    bottom=False,
                )

        if full:
            handles = [l1, l2, l3]
            if l4 is not None:
                handles.append(l4)
            labels = [h.get_label() for h in handles]
            ax.legend(handles, labels, frameon=False, loc="center right")

            text = (
                f"Official: {official_x:g}\n"
                f"P={off_p:.3f}, R={off_r:.3f}, F1={off_f1:.3f}\n\n"
                f"Best F1: {best_x:g}\n"
                f"P={best_p:.3f}, R={best_r:.3f}, F1={best_f1:.3f}"
            )

            ax.text(
                0.98,
                0.04,
                text,
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=6.5,
                bbox=dict(facecolor="white", edgecolor="0.45", alpha=0.88, pad=2.5),
            )

        save_all(fig, args, outdirs, outbase, full=full)


def pivot_heatmap(df, x_col, y_col, value_col):
    tmp = df[[x_col, y_col, value_col]].copy()
    tmp[x_col] = pd.to_numeric(tmp[x_col], errors="coerce")
    tmp[y_col] = pd.to_numeric(tmp[y_col], errors="coerce")
    tmp[value_col] = pd.to_numeric(tmp[value_col], errors="coerce")
    tmp = tmp.dropna(subset=[x_col, y_col, value_col])

    mat = tmp.pivot_table(index=y_col, columns=x_col, values=value_col, aggfunc="mean")
    mat = mat.sort_index(ascending=True)
    mat = mat.reindex(sorted(mat.columns), axis=1)
    return mat


def plot_joint_heatmap(args, outdirs, value_name="F1", outbase="C_SWARM_joint_F1_heatmap"):
    df = pd.read_csv(args.joint_tsv, sep="\t")

    prob_col = maybe_col(
        df,
        ["probability_threshold", "prob_threshold", "prob_th", "p_threshold"],
        contains_any=["prob", "threshold"],
    )
    stoich_col = maybe_col(
        df,
        ["stoichiometry_threshold", "stoich_threshold", "stoich_th"],
        contains_any=["stoich", "threshold"],
    )
    val_col = find_col(df, [value_name, value_name.lower(), value_name.upper()])
    precision_col = maybe_col(df, ["Precision", "precision"])
    recall_col = maybe_col(df, ["Recall", "recall"])

    if prob_col is None or stoich_col is None:
        threshold_cols = [c for c in df.columns if "threshold" in norm_col(c)]
        if len(threshold_cols) < 2:
            raise KeyError(f"Cannot identify threshold columns in {args.joint_tsv}: {df.columns.tolist()}")
        prob_col = threshold_cols[0]
        stoich_col = threshold_cols[1]

    for c in [prob_col, stoich_col, val_col]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    if precision_col:
        df[precision_col] = pd.to_numeric(df[precision_col], errors="coerce")
    if recall_col:
        df[recall_col] = pd.to_numeric(df[recall_col], errors="coerce")

    df = df.dropna(subset=[prob_col, stoich_col, val_col]).copy()

    best = df.sort_values([val_col], ascending=False).iloc[0]
    best_prob = float(best[prob_col])
    best_stoich = float(best[stoich_col])
    best_val = float(best[val_col])
    best_p = float(best[precision_col]) if precision_col else np.nan
    best_r = float(best[recall_col]) if recall_col else np.nan

    mat = pivot_heatmap(df, prob_col, stoich_col, val_col)
    mat.to_csv(outdirs[2] / f"{outbase}.matrix.tsv", sep="\t")

    xvals = np.array(mat.columns, dtype=float)
    yvals = np.array(mat.index, dtype=float)

    vmin = np.floor(np.nanmin(mat.values) * 100) / 100
    vmax = np.ceil(np.nanmax(mat.values) * 100) / 100
    if vmax <= vmin:
        vmax = vmin + 0.01

    for full in [True, False]:
        fig, ax = plt.subplots(figsize=(args.heatmap_width, args.heatmap_height))

        im = ax.imshow(
            mat.values,
            origin="lower",
            aspect="auto",
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            extent=[xvals.min(), xvals.max(), yvals.min(), yvals.max()],
        )

        ax.scatter(
            [args.prob_official],
            [args.stoich_official],
            s=38,
            marker="o",
            color="#1f77b4",
            edgecolor="white",
            linewidth=0.4,
            zorder=5,
        )
        ax.scatter(
            [best_prob],
            [best_stoich],
            s=78,
            marker="*",
            color="#ff7f0e",
            edgecolor="white",
            linewidth=0.35,
            zorder=6,
        )

        xticks = np.linspace(xvals.min(), xvals.max(), 5)
        yticks = np.linspace(yvals.min(), yvals.max(), 5)

        ax.set_xlim(xvals.min(), xvals.max())
        ax.set_ylim(yvals.min(), yvals.max())
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.set_xticklabels([f"{x:.4g}" for x in xticks], rotation=45, ha="right")
        ax.set_yticklabels([f"{y:.2g}" for y in yticks])

        if full:
            ax.set_title(f"SWARM joint threshold scan ({value_name} heatmap)", pad=5)
            ax.set_xlabel("Probability threshold")
            ax.set_ylabel("Stoichiometry threshold")

            text = (
                f"Best F1\n"
                f"prob={best_prob:.4g}, stoich={best_stoich:.2g}\n"
                f"F1={best_val:.3f}, P={best_p:.3f}, R={best_r:.3f}\n\n"
                f"Circle = official cutoff\n"
                f"prob={args.prob_official}, stoich={args.stoich_official}"
            )

            ax.text(
                0.04,
                0.06,
                text,
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=6.5,
                bbox=dict(facecolor="white", edgecolor="0.45", alpha=0.88, pad=2.5),
            )
        else:
            if args.clean_keep_axis_labels:
                ax.set_xlabel("Probability threshold")
                ax.set_ylabel("Stoichiometry threshold")
            else:
                ax.set_xlabel("")
                ax.set_ylabel("")
            fig.subplots_adjust(bottom=0.22)

        style_general_ax(ax, args)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=args.tick_labelsize)
        if full:
            cbar.set_label(value_name)
        else:
            cbar.set_label("")
            cbar.ax.set_yticklabels([])

        save_all(fig, args, outdirs, outbase, full=full)


def plot_strategy_metrics(args, outdirs):
    df = pd.read_csv(args.strategy_tsv, sep="\t")

    strategy_col = maybe_col(df, ["strategy", "method", "Strategy", "Method", "call_set"])
    if strategy_col is None:
        strategy_col = df.columns[0]

    precision_col, recall_col, f1_col, _ = get_metric_cols(df)

    for c in [precision_col, recall_col, f1_col]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df[[strategy_col, precision_col, recall_col, f1_col]].to_csv(
        outdirs[2] / "D_SWARM_strategy_metrics.used_table.tsv",
        sep="\t",
        index=False,
    )

    x = np.arange(len(df))
    width = 0.24

    for full in [True, False]:
        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))

        ax.bar(x - width, df[precision_col], width, label="Precision")
        ax.bar(x, df[recall_col], width, label="Recall")
        ax.bar(x + width, df[f1_col], width, label="F1")

        ax.set_ylim(0, 1)
        ax.set_yticks(tick_values(args))
        ax.grid(False)

        if full:
            ax.set_title("SWARM and m6Anet strategy comparison", pad=5)
            ax.set_ylabel("Metric")
            ax.set_xticks(x)
            ax.set_xticklabels(df[strategy_col].astype(str), rotation=35, ha="right")
            ax.legend(frameon=False, loc="best")
        else:
            ax.set_ylabel("")
            ax.set_xticks(x)
            ax.set_xticklabels([])

        style_general_ax(ax, args)
        save_all(fig, args, outdirs, "D_SWARM_strategy_metrics", full=full)


def main():
    args = parse_args()
    configure_paths(args)
    configure_matplotlib(args)
    outdirs = make_output_dirs(args)

    plot_threshold_scan_with_predicted(
        args,
        outdirs,
        args.prob_tsv,
        x_candidates=["probability_threshold", "prob_threshold", "threshold"],
        title="SWARM probability threshold scan",
        xlabel="Probability threshold",
        outbase="A_SWARM_probability_scan",
        official_x=args.prob_official,
    )

    plot_threshold_scan_with_predicted(
        args,
        outdirs,
        args.stoich_tsv,
        x_candidates=["stoichiometry_threshold", "stoich_threshold", "threshold"],
        title="SWARM stoichiometry threshold scan",
        xlabel="Stoichiometry threshold",
        outbase="B_SWARM_stoichiometry_scan",
        official_x=args.stoich_official,
    )

    plot_joint_heatmap(
        args,
        outdirs,
        value_name=args.heatmap_value,
        outbase=f"C_SWARM_joint_{args.heatmap_value}_heatmap",
    )

    plot_strategy_metrics(args, outdirs)

    print("Done.")
    print("Output directory:", args.outdir)
    print("Full figures:", outdirs[0])
    print("Clean figures:", outdirs[1])
    print("Tables:", outdirs[2])


if __name__ == "__main__":
    main()
