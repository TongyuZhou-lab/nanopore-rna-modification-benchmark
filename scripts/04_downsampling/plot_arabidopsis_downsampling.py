#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot Arabidopsis m6A downsampling / read-depth robustness figures.

This command-line version removes hard-coded paths from the original script.
It reads two prepared downsampling result tables and generates unified-style
full and clean figures.

Main panels:
  A. Reference-hit recall
  B. Full-depth TP-site retention
  C. Recovered reference-hit count
  D. Positive-call composition

Expected outputs:
  OUTDIR/full/
  OUTDIR/clean/
  OUTDIR/tables/

Example:
python plot_arabidopsis_downsampling.py \
  --acd-table path/to/downsample_four_figure_metrics.v3.tsv \
  --b-table path/to/B_fullDepth_TPsite_retention.v3.tsv \
  --outdir path/to/downsample_single_unified
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


TOOL_ORDER_DEFAULT = ["TandemMod", "m6Anet", "CHEUI", "SWARM"]

# Keep the original line colors.
TOOL_COLORS = {
    "TandemMod": "#7E57C2",
    "m6Anet": "#2E86DE",
    "CHEUI": "#009E73",
    "SWARM": "#E67E22",
}

# Keep the original bar colors.
BAR_REF_COLOR = "#4C78A8"
BAR_NONREF_COLOR = "#F58518"

FRACTION_ORDER_DEFAULT = [
    "full_depth",
    "frac0.8_seed11",
    "frac0.6_seed11",
    "frac0.4_seed11",
    "frac0.2_seed11",
]

FRACTION_LABELS_DEFAULT = {
    "full_depth": "Full",
    "frac0.8_seed11": "80%",
    "frac0.6_seed11": "60%",
    "frac0.4_seed11": "40%",
    "frac0.2_seed11": "20%",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Arabidopsis m6A downsampling figures without hard-coded paths."
    )

    parser.add_argument(
        "--acd-table",
        required=True,
        type=Path,
        help=(
            "Input table containing metrics for panels A/C/D, "
            "for example downsample_four_figure_metrics.v3.tsv."
        ),
    )
    parser.add_argument(
        "--b-table",
        required=True,
        type=Path,
        help=(
            "Input table containing full-depth TP-site retention for panel B, "
            "for example B_fullDepth_TPsite_retention.v3.tsv."
        ),
    )
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")

    # Figure style
    parser.add_argument("--fig-width", type=float, default=4.0, help="Single-figure width.")
    parser.add_argument("--fig-height", type=float, default=4.3, help="Single-figure height.")
    parser.add_argument("--combined-width", type=float, default=8.0, help="Combined 2x2 figure width.")
    parser.add_argument("--combined-height", type=float, default=8.0, help="Combined 2x2 figure height.")
    parser.add_argument("--line-width", type=float, default=1.8, help="Line width.")
    parser.add_argument("--marker-size", type=float, default=4.0, help="Line marker size.")
    parser.add_argument("--axis-line-width", type=float, default=0.8, help="Axis spine width.")
    parser.add_argument("--tick-line-width", type=float, default=0.8, help="Tick line width.")
    parser.add_argument("--tick-length", type=float, default=3.0, help="Tick length.")
    parser.add_argument("--tick-labelsize", type=float, default=12, help="Tick-label font size.")
    parser.add_argument("--axis-labelsize", type=float, default=12, help="Axis-label font size.")
    parser.add_argument("--title-size", type=float, default=12, help="Title font size for full figures.")
    parser.add_argument("--legend-size", type=float, default=8, help="Legend font size for full figures.")
    parser.add_argument("--dpi", type=int, default=300, help="PNG output DPI.")

    parser.add_argument(
        "--clean-hide-axis-labels",
        action="store_true",
        help=(
            "Hide x/y axis labels in clean figures. "
            "By default, clean figures remove titles and legends but keep axis labels and tick labels."
        ),
    )
    parser.add_argument(
        "--keep-full-box",
        action="store_true",
        help="Keep all four spines. By default, only left and bottom spines are kept.",
    )
    parser.add_argument(
        "--tool-order",
        default=",".join(TOOL_ORDER_DEFAULT),
        help="Comma-separated tool order. Default: TandemMod,m6Anet,CHEUI,SWARM",
    )
    parser.add_argument(
        "--fraction-order",
        default=",".join(FRACTION_ORDER_DEFAULT),
        help=(
            "Comma-separated fraction order. Default: "
            "full_depth,frac0.8_seed11,frac0.6_seed11,frac0.4_seed11,frac0.2_seed11"
        ),
    )
    parser.add_argument(
        "--fraction-labels",
        default="",
        help=(
            "Optional comma-separated mapping for fraction labels, format: key=label,key=label. "
            "If omitted, default labels are used."
        ),
    )

    # Panel y-axis options
    parser.add_argument("--a-ymin", type=float, default=0.0, help="Panel A y-axis minimum.")
    parser.add_argument("--a-ymax", type=float, default=0.8, help="Panel A y-axis maximum.")
    parser.add_argument("--b-ymin", type=float, default=0.0, help="Panel B y-axis minimum.")
    parser.add_argument("--b-ymax", type=float, default=1.05, help="Panel B y-axis maximum.")
    parser.add_argument(
        "--c-auto-y",
        action="store_true",
        help="Use automatic y-axis range for panel C. This is the default behavior.",
    )

    return parser.parse_args()


def parse_csv_list(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def parse_fraction_labels(value: str) -> dict[str, str]:
    labels = dict(FRACTION_LABELS_DEFAULT)
    if not value:
        return labels

    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"Invalid --fraction-labels item: {item}. Expected key=label.")
        key, label = item.split("=", 1)
        labels[key.strip()] = label.strip()
    return labels


def check_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def configure_matplotlib(args: argparse.Namespace) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": args.axis_labelsize,
        "axes.titlesize": args.title_size,
        "axes.labelsize": args.axis_labelsize,
        "xtick.labelsize": args.tick_labelsize,
        "ytick.labelsize": args.tick_labelsize,
        "legend.fontsize": args.legend_size,
        "axes.linewidth": args.axis_line_width,
        "axes.grid": False,
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


def save_all(fig, args: argparse.Namespace, outdir: Path, base: str):
    """Save transparent SVG/PDF/PNG while preserving the figure canvas ratio."""
    fig.patch.set_alpha(0)
    fig.savefig(outdir / f"{base}.svg", transparent=True)
    fig.savefig(outdir / f"{base}.pdf", transparent=True)
    fig.savefig(outdir / f"{base}.png", dpi=args.dpi, transparent=True)
    plt.close(fig)


def norm_col(c):
    return str(c).strip().lower().replace(" ", "_").replace("-", "_")


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


def standardize_fraction(df, fraction_order, fraction_labels):
    df = df.copy()

    frac_col = maybe_col(df, ["fraction", "fraction_tag", "frac", "sample"])
    if frac_col is None:
        raise KeyError(f"Cannot find fraction column in {df.columns.tolist()}")

    df = df.rename(columns={frac_col: "fraction"})

    label_col = maybe_col(df, ["fraction_label", "label"])
    if label_col:
        df["fraction_label"] = df[label_col].astype(str)
    else:
        df["fraction_label"] = df["fraction"].map(fraction_labels).fillna(df["fraction"].astype(str))

    df["fraction"] = df["fraction"].astype(str)
    df["fraction_rank"] = df["fraction"].map({v: i for i, v in enumerate(fraction_order)})

    # Fallback if fraction names differ.
    if df["fraction_rank"].isna().any():
        def infer_rank(x):
            s = str(x).lower()
            if "full" in s:
                return 0
            if "0.8" in s or "80" in s:
                return 1
            if "0.6" in s or "60" in s:
                return 2
            if "0.4" in s or "40" in s:
                return 3
            if "0.2" in s or "20" in s:
                return 4
            return np.nan

        df["fraction_rank"] = df["fraction_rank"].fillna(df["fraction"].map(infer_rank))

    df = df.dropna(subset=["fraction_rank"]).copy()
    df["fraction_rank"] = df["fraction_rank"].astype(int)
    df["fraction_label"] = df["fraction"].map(fraction_labels).fillna(df["fraction_label"])
    return df


def standardize_tool(df, tool_order):
    df = df.copy()
    tool_col = maybe_col(df, ["tool", "method"])
    if tool_col is None:
        raise KeyError(f"Cannot find tool column in {df.columns.tolist()}")
    df = df.rename(columns={tool_col: "tool"})
    df["tool"] = df["tool"].astype(str)
    df = df[df["tool"].isin(tool_order)].copy()
    df["tool_rank"] = df["tool"].map({tool: i for i, tool in enumerate(tool_order)})
    return df


def style_ax(args: argparse.Namespace, ax, xlabel, ylabel, title=None, show_title=False, y_lim=None, clean=False):
    """Shared axis style for both full and clean versions."""
    if show_title and title:
        ax.set_title(title, pad=5)
    else:
        ax.set_title("")

    if clean and args.clean_hide_axis_labels:
        ax.set_xlabel("")
        ax.set_ylabel("")
    else:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    if y_lim is not None:
        ax.set_ylim(*y_lim)

    ax.grid(False)
    ax.set_facecolor("none")

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
        top=top,
        right=right,
        direction="out",
    )


def line_plot_by_tool(args, outdirs, df, y_col, title, ylabel, outbase, y_lim=None, legend_loc="best"):
    full_dir, clean_dir, table_dir = outdirs
    tool_order = parse_csv_list(args.tool_order)
    fraction_order = parse_csv_list(args.fraction_order)
    fraction_labels = parse_fraction_labels(args.fraction_labels)

    df = standardize_fraction(standardize_tool(df, tool_order), fraction_order, fraction_labels)
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")

    used_cols = ["fraction", "fraction_label", "fraction_rank", "tool", y_col]
    if "tool_rank" in df.columns:
        used_cols.append("tool_rank")
    used = df[used_cols].copy()
    used = used.sort_values(["tool_rank" if "tool_rank" in used.columns else "tool", "fraction_rank"])
    used.to_csv(table_dir / f"{outbase}.used_table.tsv", sep="\t", index=False)

    labels = (
        df[["fraction_rank", "fraction_label"]]
        .drop_duplicates()
        .sort_values("fraction_rank")
    )

    for full in [True, False]:
        fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
        fig.patch.set_alpha(0)

        for tool in tool_order:
            sub = df[df["tool"] == tool].sort_values("fraction_rank")
            if sub.empty:
                continue
            ax.plot(
                sub["fraction_rank"],
                sub[y_col],
                marker="o",
                markersize=args.marker_size,
                linewidth=args.line_width,
                label=tool,
                color=TOOL_COLORS.get(tool, None),
            )

        ax.set_xticks(labels["fraction_rank"])
        ax.set_xticklabels(labels["fraction_label"])

        style_ax(
            args,
            ax,
            title=title,
            xlabel="Reads fraction",
            ylabel=ylabel,
            show_title=full,
            y_lim=y_lim,
            clean=not full,
        )

        # Full version keeps legend; clean version removes title and legend.
        if full:
            ax.legend(frameon=False, loc=legend_loc)

        fig.tight_layout()
        outdir = full_dir if full else clean_dir
        suffix = "full" if full else "clean"
        save_all(fig, args, outdir, f"{outbase}.{suffix}")


def positive_call_composition_per_tool(args, outdirs, df):
    full_dir, clean_dir, table_dir = outdirs
    tool_order = parse_csv_list(args.tool_order)
    fraction_order = parse_csv_list(args.fraction_order)
    fraction_labels = parse_fraction_labels(args.fraction_labels)

    df = standardize_fraction(standardize_tool(df, tool_order), fraction_order, fraction_labels)

    ref_col = find_col(df, ["TP_primary", "reference_supported_calls", "reference_calls", "TP"])
    nonref_col = find_col(df, ["non_reference_calls_primary", "non_reference_calls", "nonref_calls"])
    pos_col = maybe_col(df, ["n_primary_eval_positive_calls", "positive_calls", "n_positive_calls"])

    df[ref_col] = pd.to_numeric(df[ref_col], errors="coerce").fillna(0)
    df[nonref_col] = pd.to_numeric(df[nonref_col], errors="coerce").fillna(0)

    if pos_col:
        df[pos_col] = pd.to_numeric(df[pos_col], errors="coerce")
    else:
        df["total_positive_calls"] = df[ref_col] + df[nonref_col]
        pos_col = "total_positive_calls"

    used_cols = ["fraction", "fraction_label", "fraction_rank", "tool", ref_col, nonref_col, pos_col]
    if "tool_rank" in df.columns:
        used_cols.append("tool_rank")
    used = df[used_cols].copy()
    used = used.sort_values(["tool_rank" if "tool_rank" in used.columns else "tool", "fraction_rank"])
    used.to_csv(table_dir / "D_positive_call_composition.used_table.tsv", sep="\t", index=False)

    # Per-tool single figures.
    for tool in tool_order:
        sub = df[df["tool"] == tool].sort_values("fraction_rank")
        if sub.empty:
            continue

        for full in [True, False]:
            fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
            fig.patch.set_alpha(0)

            x = sub["fraction_rank"].to_numpy()
            ref = sub[ref_col].to_numpy()
            nonref = sub[nonref_col].to_numpy()

            ax.bar(x, ref, width=0.65, label="Reference-supported calls", color=BAR_REF_COLOR)
            ax.bar(x, nonref, bottom=ref, width=0.65, label="Non-reference calls", color=BAR_NONREF_COLOR)

            ax.set_xticks(sub["fraction_rank"])
            ax.set_xticklabels(sub["fraction_label"])

            style_ax(
                args,
                ax,
                title=f"{tool}: positive-call composition",
                xlabel="Reads fraction",
                ylabel="Positive calls",
                show_title=full,
                y_lim=None,
                clean=not full,
            )

            if full:
                ax.legend(frameon=False, loc="best")

            fig.tight_layout()
            outdir = full_dir if full else clean_dir
            suffix = "full" if full else "clean"
            safe_tool = tool.replace(" ", "_")
            save_all(fig, args, outdir, f"D_positive_call_composition_{safe_tool}.{suffix}")

    # Combined 2x2 figure as optional.
    for full in [True, False]:
        fig, axes = plt.subplots(2, 2, figsize=(args.combined_width, args.combined_height))
        fig.patch.set_alpha(0)
        axes = axes.ravel()

        for ax, tool in zip(axes, tool_order):
            sub = df[df["tool"] == tool].sort_values("fraction_rank")
            if sub.empty:
                ax.axis("off")
                continue

            x = sub["fraction_rank"].to_numpy()
            ref = sub[ref_col].to_numpy()
            nonref = sub[nonref_col].to_numpy()

            ax.bar(x, ref, width=0.65, label="Reference-supported calls", color=BAR_REF_COLOR)
            ax.bar(x, nonref, bottom=ref, width=0.65, label="Non-reference calls", color=BAR_NONREF_COLOR)
            ax.set_xticks(sub["fraction_rank"])
            ax.set_xticklabels(sub["fraction_label"], rotation=0)

            style_ax(
                args,
                ax,
                title=tool,
                xlabel="Reads fraction",
                ylabel="Positive calls",
                show_title=full,
                y_lim=None,
                clean=not full,
            )

        if full and len(axes) > 0:
            handles, labels = axes[0].get_legend_handles_labels()
            fig.legend(
                handles,
                labels,
                frameon=False,
                loc="upper center",
                ncol=2,
                bbox_to_anchor=(0.5, 0.985),
            )
            fig.tight_layout(rect=(0, 0, 1, 0.94))
        else:
            fig.tight_layout()

        outdir = full_dir if full else clean_dir
        suffix = "full" if full else "clean"
        save_all(fig, args, outdir, f"D_positive_call_composition_combined.{suffix}")


def main():
    args = parse_args()
    check_file(args.acd_table, "--acd-table")
    check_file(args.b_table, "--b-table")

    configure_matplotlib(args)
    outdirs = make_output_dirs(args)

    acd = pd.read_csv(args.acd_table, sep="\t")
    b = pd.read_csv(args.b_table, sep="\t")

    # A. Reference-hit recall
    recall_col = find_col(acd, ["Reference_hit_recall", "reference_hit_recall"])
    line_plot_by_tool(
        args,
        outdirs,
        acd,
        y_col=recall_col,
        title="Reference-hit recall across read fractions",
        ylabel="Reference-hit recall",
        outbase="A_reference_hit_recall",
        y_lim=(args.a_ymin, args.a_ymax),
        legend_loc="best",
    )

    # B. Full-depth TP-site retention
    b_col = find_col(b, ["full_depth_TP_site_retention", "Full_depth_TP_site_retention"])
    line_plot_by_tool(
        args,
        outdirs,
        b,
        y_col=b_col,
        title="Full-depth TP-site retention",
        ylabel="Full-depth TP-site retention",
        outbase="B_fullDepth_TPsite_retention",
        y_lim=(args.b_ymin, args.b_ymax),
        legend_loc="best",
    )

    # C. Recovered reference-hit count
    tp_col = find_col(acd, ["TP_primary", "recovered_reference_hits", "TP"])
    line_plot_by_tool(
        args,
        outdirs,
        acd,
        y_col=tp_col,
        title="Recovered reference-hit count",
        ylabel="Recovered reference hits",
        outbase="C_recovered_reference_hits_TP_count",
        y_lim=None,
        legend_loc="best",
    )

    # D. Positive-call composition
    positive_call_composition_per_tool(args, outdirs, acd)

    print("Done.")
    print("Output directory:", args.outdir)
    print("Full figures:", outdirs[0])
    print("Clean figures:", outdirs[1])
    print("Tables:", outdirs[2])


if __name__ == "__main__":
    main()
