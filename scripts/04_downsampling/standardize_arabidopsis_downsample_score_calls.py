#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Standardize Arabidopsis full-depth and downsampled m6A tool outputs.

This command-line version removes hard-coded benchmark paths from the original
standardization script. It converts raw outputs from TandemMod, m6Anet, CHEUI,
and SWARM into unified genome-level score/call tables and positive-call site
lists for downstream downsampling analysis.

Required inputs
---------------
1. A label table with site_id and label columns.
2. A full-depth common-universe site list.
3. A GTF file for transcript-to-genome coordinate conversion.
4. A raw-path table listing input files for each fraction/tool pair.

Raw-path table format
---------------------
TSV with at least these columns:

    fraction    tool    path

Example:

    full_depth        TandemMod    /path/to/full_depth/TandemMod/m6A.site_level.tsv
    full_depth        m6Anet       /path/to/full_depth/m6Anet/data.site_proba.csv
    full_depth        CHEUI        /path/to/full_depth/CHEUI/site_level_m6A_predictions.txt
    full_depth        SWARM        /path/to/full_depth/SWARM/pred.tsv
    frac0.8_seed11    TandemMod    /path/to/frac0.8_seed11/TandemMod/m6A.site_level.tsv

Main outputs
------------
OUTDIR/
  score_call_tables/
  positive_call_site_lists/
  qc/

Example
-------
python standardize_arabidopsis_downsample_score_calls.py \
  --label-path path/to/labels.primary_nanopore.full_depth_common_universe.tsv \
  --universe-path path/to/full_depth_four_tool_common_universe.site_ids.tsv \
  --raw-paths path/to/arabidopsis_downsample_raw_paths.tsv \
  --gtf path/to/Arabidopsis_thaliana.TAIR10.60.gtf \
  --outdir path/to/standardized_full_and_downsample_score_call_tables
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


DRACH_RE = re.compile(r"[AGT][AG]AC[ACT]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Standardize Arabidopsis full-depth and downsampled outputs from "
            "TandemMod, m6Anet, CHEUI, and SWARM into genome-level score/call tables."
        )
    )

    parser.add_argument("--label-path", required=True, type=Path, help="Label TSV with site_id and label columns.")
    parser.add_argument("--universe-path", required=True, type=Path, help="Full common-universe site list.")
    parser.add_argument(
        "--raw-paths",
        required=True,
        type=Path,
        help="TSV with columns: fraction, tool, path.",
    )
    parser.add_argument("--gtf", required=True, type=Path, help="GTF/GTF.GZ used for transcript-to-genome mapping.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")

    parser.add_argument(
        "--tools",
        default="TandemMod,m6Anet,CHEUI,SWARM",
        help="Comma-separated tools to process. Default: TandemMod,m6Anet,CHEUI,SWARM",
    )
    parser.add_argument(
        "--fractions",
        default="full_depth,frac0.8_seed11,frac0.6_seed11,frac0.4_seed11,frac0.2_seed11",
        help=(
            "Comma-separated fractions to process. Default: "
            "full_depth,frac0.8_seed11,frac0.6_seed11,frac0.4_seed11,frac0.2_seed11"
        ),
    )

    parser.add_argument("--m6anet-threshold", type=float, default=0.9, help="m6Anet call threshold. Default: 0.9")
    parser.add_argument("--cheui-threshold", type=float, default=0.9999, help="CHEUI call threshold. Default: 0.9999")
    parser.add_argument("--swarm-prob-threshold", type=float, default=0.9972, help="SWARM probability threshold. Default: 0.9972")
    parser.add_argument("--swarm-stoich-threshold", type=float, default=0.1, help="SWARM stoichiometry threshold. Default: 0.1")
    parser.add_argument("--tandemmod-p95-min", type=float, default=10, help="TandemMod p_0.95 minimum. Default: >10")
    parser.add_argument("--tandemmod-ratio-min", type=float, default=0.2, help="TandemMod p_0.95/total minimum. Default: >=0.2")

    parser.add_argument(
        "--expected-full-depth",
        type=Path,
        default=None,
        help=(
            "Optional sanity-check TSV with columns: tool, expected_n_primary_eval_positive_calls, expected_TP_primary. "
            "If omitted, the expected-count sanity check is skipped."
        ),
    )
    parser.add_argument(
        "--no-drach-filter",
        action="store_true",
        help="Disable DRACH motif filtering even when motif/kmer columns are available.",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Skip missing fraction/tool input files instead of raising an error.",
    )

    return parser.parse_args()


def parse_csv_list(x: str) -> list[str]:
    return [v.strip() for v in str(x).split(",") if v.strip()]


def check_file(path: Path, name: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{name} not found: {path}")


def sniff_sep(path: Path) -> str:
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if line.strip():
                return "," if line.count(",") > line.count("\t") else "\t"
    return "\t"


def read_table(path: Path) -> pd.DataFrame:
    sep = sniff_sep(path)
    return pd.read_csv(path, sep=sep, low_memory=False)


def normalize_chrom(x):
    s = str(x)
    s = s.replace("Chr", "").replace("chr", "")
    return s


def to_bool(s):
    if s.dtype == object:
        return s.astype(str).str.lower().isin(["true", "1", "yes"])
    return s.fillna(False).astype(bool)


def safe_numeric(s):
    return pd.to_numeric(s, errors="coerce")


def is_drach(x):
    if pd.isna(x):
        return False
    s = str(x).upper().replace("U", "T")
    if len(s) < 5:
        return False
    return any(DRACH_RE.fullmatch(s[i:i + 5]) for i in range(len(s) - 4))


def maybe_filter_drach(df: pd.DataFrame, candidates: list[str], tool: str, no_drach_filter: bool) -> pd.DataFrame:
    if no_drach_filter:
        print(f"  {tool}: DRACH filter disabled; keeping all rows")
        return df.copy()

    for c in candidates:
        if c in df.columns:
            before = len(df)
            keep = df[c].map(is_drach)
            after = int(keep.sum())
            if after > 0:
                print(f"  {tool}: DRACH filter by column '{c}': {before} -> {after}")
                return df.loc[keep].copy()

    print(f"  {tool}: no usable DRACH motif column found; keeping all rows")
    return df.copy()


def pick_col(df: pd.DataFrame, candidates: list[str], tool: str, desc: str) -> str:
    cols = list(df.columns)
    lower = {str(c).lower(): c for c in cols}
    for c in candidates:
        if c in cols:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    raise ValueError(f"{tool}: cannot find {desc}. candidates={candidates}; existing columns={cols}")


def load_site_set(path: Path) -> set[str]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "site_id" in df.columns:
        return set(df["site_id"].astype(str))
    if {"chrom", "start"}.issubset(df.columns):
        return set(df["chrom"].astype(str) + ":" + df["start"].astype(int).astype(str))
    if df.shape[1] == 1:
        return set(df.iloc[:, 0].astype(str))
    raise ValueError(f"Cannot parse site set: {path}")


def parse_attrs(attr):
    d = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', str(attr)):
        d[m.group(1)] = m.group(2)
    return d


def load_gtf_models(gtf_path: Path):
    print(f"[GTF] Loading {gtf_path}")
    opener = gzip.open if str(gtf_path).endswith(".gz") else open

    exons = {}
    with opener(gtf_path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attr = parts
            if feature != "exon":
                continue

            attrs = parse_attrs(attr)
            tx = attrs.get("transcript_id")
            if not tx:
                continue

            chrom = normalize_chrom(chrom)
            start0 = int(start) - 1
            end0 = int(end)
            exons.setdefault(tx, []).append((chrom, start0, end0, strand))

    models = {}
    alias = {}

    for tx, items in exons.items():
        strand = items[0][3]
        items = sorted(items, key=lambda x: x[1], reverse=(strand == "-"))

        segs = []
        cur = 0
        for chrom, start0, end0, strand in items:
            length = end0 - start0
            segs.append({
                "tx_start": cur,
                "tx_end": cur + length,
                "chrom": chrom,
                "start0": start0,
                "end0": end0,
                "strand": strand,
            })
            cur += length

        models[tx] = segs
        no_version = tx.split(".")[0]
        alias.setdefault(no_version, set()).add(tx)

    alias_unique = {}
    for k, v in alias.items():
        if len(v) == 1:
            alias_unique[k] = list(v)[0]

    print(f"[GTF] transcript models: {len(models)}")
    print(f"[GTF] unique no-version aliases: {len(alias_unique)}")
    return models, alias_unique


def tx_lookup(tx, models, alias_unique):
    tx = str(tx)
    if tx in models:
        return tx
    if tx in alias_unique:
        return alias_unique[tx]
    no_version = tx.split(".")[0]
    if no_version in alias_unique:
        return alias_unique[no_version]
    return None


def tx_to_genome(tx, pos, models, alias_unique, origin=0):
    tx2 = tx_lookup(tx, models, alias_unique)
    if tx2 is None:
        return None

    try:
        p = int(float(pos))
    except Exception:
        return None

    q = p if origin == 0 else p - 1
    if q < 0:
        return None

    for seg in models[tx2]:
        if seg["tx_start"] <= q < seg["tx_end"]:
            offset = q - seg["tx_start"]
            if seg["strand"] == "-":
                g = seg["end0"] - 1 - offset
            else:
                g = seg["start0"] + offset
            return seg["chrom"], int(g)

    return None


def aggregate_by_max_score_any_call(df: pd.DataFrame) -> pd.DataFrame:
    df = df.dropna(subset=["site_id"]).copy()
    df["site_id"] = df["site_id"].astype(str)
    df["score"] = safe_numeric(df["score"]).fillna(0)
    df["call"] = to_bool(df["call"])

    idx = df.groupby("site_id")["score"].idxmax()
    best = df.loc[idx].copy()

    any_call = df.groupby("site_id", as_index=False)["call"].max()
    best = best.drop(columns=["call"]).merge(any_call, on="site_id", how="left")

    return best


def standardize_tandemmod(path: Path, fraction: str, args: argparse.Namespace):
    tool = "TandemMod"
    df = read_table(path)

    required = ["chr", "site.1", "p_0.95", "total"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"{tool} {fraction}: missing required column {c}; columns={list(df.columns)}")

    df = maybe_filter_drach(df, ["motif", "kmer", "site"], tool, args.no_drach_filter)

    df = df.copy()
    df["chrom"] = df["chr"].map(normalize_chrom)
    df["start"] = safe_numeric(df["site.1"]) - 1
    df["p95"] = safe_numeric(df["p_0.95"]).fillna(0)
    df["total_num"] = safe_numeric(df["total"]).fillna(0)

    df = df.dropna(subset=["chrom", "start"]).copy()
    df = df[df["total_num"] > 0].copy()
    df["start"] = df["start"].astype(int)
    df["site_id"] = df["chrom"].astype(str) + ":" + df["start"].astype(str)
    df["ratio"] = (df["p95"] / df["total_num"]).replace([np.inf, -np.inf], np.nan).fillna(0)

    grouped = df.groupby("site_id").agg(
        chrom=("chrom", "first"),
        start=("start", "first"),
        score=("ratio", "max"),
        max_p95=("p95", "max"),
        max_ratio=("ratio", "max"),
        max_total=("total_num", "max"),
        n_source_rows=("site_id", "size"),
    ).reset_index()

    grouped["call"] = (grouped["max_p95"] > args.tandemmod_p95_min) & (grouped["max_ratio"] >= args.tandemmod_ratio_min)
    return grouped


def standardize_m6anet(path: Path, models, alias_unique, fraction: str, args: argparse.Namespace):
    tool = "m6Anet"
    df = read_table(path)

    tx_col = pick_col(df, ["transcript_id", "transcript", "contig", "tx"], tool, "transcript column")
    pos_col = pick_col(df, ["transcript_position", "position", "pos"], tool, "position column")
    score_col = pick_col(df, ["probability_modified", "prob_modified", "probability", "score"], tool, "score column")

    df = maybe_filter_drach(df, ["kmer", "motif", "site"], tool, args.no_drach_filter)

    df = df.copy()
    df["score"] = safe_numeric(df[score_col]).fillna(0)
    df["call"] = df["score"] >= args.m6anet_threshold

    mapped = [tx_to_genome(tx, pos, models, alias_unique, origin=0) for tx, pos in zip(df[tx_col], df[pos_col])]

    df["chrom"] = [x[0] if x else None for x in mapped]
    df["start"] = [x[1] if x else None for x in mapped]

    df = df.dropna(subset=["chrom", "start"]).copy()
    df["chrom"] = df["chrom"].map(normalize_chrom)
    df["start"] = df["start"].astype(int)
    df["site_id"] = df["chrom"].astype(str) + ":" + df["start"].astype(str)

    keep = ["site_id", "chrom", "start", "score", "call"]
    for extra in [tx_col, pos_col, "kmer", "motif"]:
        if extra in df.columns and extra not in keep:
            keep.append(extra)

    return aggregate_by_max_score_any_call(df[keep])


def standardize_cheui(path: Path, models, alias_unique, fraction: str, args: argparse.Namespace):
    tool = "CHEUI"
    df = read_table(path)

    tx_col = pick_col(df, ["contig", "transcript_id", "transcript", "tx"], tool, "transcript column")
    pos_col = pick_col(df, ["position", "pos", "transcript_position"], tool, "position column")
    score_col = pick_col(df, ["probability", "prob", "score"], tool, "score column")

    df = maybe_filter_drach(df, ["site", "kmer", "motif"], tool, args.no_drach_filter)

    df = df.copy()
    df["score"] = safe_numeric(df[score_col]).fillna(0)
    df["call"] = df["score"] >= args.cheui_threshold

    mapped = [tx_to_genome(tx, pos, models, alias_unique, origin=1) for tx, pos in zip(df[tx_col], df[pos_col])]

    df["chrom"] = [x[0] if x else None for x in mapped]
    df["start"] = [x[1] - 5 if x else None for x in mapped]

    df = df.dropna(subset=["chrom", "start"]).copy()
    df["chrom"] = df["chrom"].map(normalize_chrom)
    df["start"] = df["start"].astype(int)
    df["site_id"] = df["chrom"].astype(str) + ":" + df["start"].astype(str)

    keep = ["site_id", "chrom", "start", "score", "call"]
    for extra in [tx_col, pos_col, "site", "kmer", "motif"]:
        if extra in df.columns and extra not in keep:
            keep.append(extra)

    return aggregate_by_max_score_any_call(df[keep])


def standardize_swarm(path: Path, models, alias_unique, fraction: str, args: argparse.Namespace):
    tool = "SWARM"
    df = read_table(path)

    tx_col = pick_col(df, ["contig", "transcript_id", "transcript", "tx"], tool, "transcript column")
    pos_col = pick_col(df, ["position", "pos", "transcript_position"], tool, "position column")
    score_col = pick_col(df, ["probability", "prob", "score"], tool, "score column")
    stoich_col = pick_col(df, ["stoichiometry", "stochiometry", "mod_ratio", "ratio"], tool, "stoichiometry column")

    df = maybe_filter_drach(df, ["site", "kmer", "motif"], tool, args.no_drach_filter)

    df = df.copy()
    df["score"] = safe_numeric(df[score_col]).fillna(0)
    df["stoichiometry"] = safe_numeric(df[stoich_col]).fillna(0)
    df["call"] = (df["score"] >= args.swarm_prob_threshold) & (df["stoichiometry"] > args.swarm_stoich_threshold)

    # Verified SWARM coordinate convention:
    # transcript position + 4, then origin=0 mapping to genome.
    pos_shifted = safe_numeric(df[pos_col]) + 4
    mapped = [tx_to_genome(tx, pos, models, alias_unique, origin=0) for tx, pos in zip(df[tx_col], pos_shifted)]

    df["chrom"] = [x[0] if x else None for x in mapped]
    df["start"] = [x[1] if x else None for x in mapped]

    df = df.dropna(subset=["chrom", "start"]).copy()
    df["chrom"] = df["chrom"].map(normalize_chrom)
    df["start"] = df["start"].astype(int)
    df["site_id"] = df["chrom"].astype(str) + ":" + df["start"].astype(str)

    keep = ["site_id", "chrom", "start", "score", "stoichiometry", "call"]
    for extra in [tx_col, pos_col, "site", "kmer", "motif"]:
        if extra in df.columns and extra not in keep:
            keep.append(extra)

    return aggregate_by_max_score_any_call(df[keep])


def standardize_one(tool: str, fraction: str, path: Path, models, alias_unique, args: argparse.Namespace):
    if tool == "TandemMod":
        return standardize_tandemmod(path, fraction, args)
    if tool == "m6Anet":
        return standardize_m6anet(path, models, alias_unique, fraction, args)
    if tool == "CHEUI":
        return standardize_cheui(path, models, alias_unique, fraction, args)
    if tool == "SWARM":
        return standardize_swarm(path, models, alias_unique, fraction, args)
    raise ValueError(f"Unsupported tool: {tool}")


def read_raw_paths_table(path: Path) -> dict[tuple[str, str], Path]:
    df = pd.read_csv(path, sep="\t")
    required = ["fraction", "tool", "path"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"--raw-paths missing columns: {missing}. Required columns: {required}")

    raw = {}
    for _, row in df.iterrows():
        fraction = str(row["fraction"]).strip()
        tool = str(row["tool"]).strip()
        p = Path(str(row["path"]).strip())
        raw[(fraction, tool)] = p
    return raw


def read_expected_full_depth(path: Path | None) -> dict[str, dict[str, int]]:
    if path is None:
        return {}

    df = pd.read_csv(path, sep="\t")
    lower = {c.lower(): c for c in df.columns}
    tool_col = lower.get("tool")
    n_col = lower.get("expected_n_primary_eval_positive_calls") or lower.get("n_called") or lower.get("expected_n_called")
    tp_col = lower.get("expected_tp_primary") or lower.get("tp") or lower.get("expected_tp")

    if tool_col is None or n_col is None or tp_col is None:
        raise ValueError(
            "--expected-full-depth must contain tool and expected count columns. "
            "Suggested columns: tool, expected_n_primary_eval_positive_calls, expected_TP_primary"
        )

    out = {}
    for _, row in df.iterrows():
        tool = str(row[tool_col]).strip()
        out[tool] = {
            "n_called": int(row[n_col]),
            "TP": int(row[tp_col]),
        }
    return out


def main():
    args = parse_args()

    check_file(args.label_path, "--label-path")
    check_file(args.universe_path, "--universe-path")
    check_file(args.raw_paths, "--raw-paths")
    check_file(args.gtf, "--gtf")
    if args.expected_full_depth is not None:
        check_file(args.expected_full_depth, "--expected-full-depth")

    tools = parse_csv_list(args.tools)
    fractions = parse_csv_list(args.fractions)

    score_call_dir = args.outdir / "score_call_tables"
    pos_call_dir = args.outdir / "positive_call_site_lists"
    qc_dir = args.outdir / "qc"
    for d in [score_call_dir, pos_call_dir, qc_dir]:
        d.mkdir(parents=True, exist_ok=True)

    raw_paths = read_raw_paths_table(args.raw_paths)
    expected_full_depth = read_expected_full_depth(args.expected_full_depth)

    for fraction in fractions:
        for tool in tools:
            p = raw_paths.get((fraction, tool))
            if p is None:
                msg = f"Missing raw path entry for fraction={fraction}, tool={tool}"
                if args.allow_missing:
                    print("WARNING:", msg)
                    continue
                raise FileNotFoundError(msg)
            if not p.exists():
                msg = f"{fraction} {tool} raw file not found: {p}"
                if args.allow_missing:
                    print("WARNING:", msg)
                    continue
                raise FileNotFoundError(msg)

    labels = pd.read_csv(args.label_path, sep="\t", low_memory=False)
    if "site_id" not in labels.columns or "label" not in labels.columns:
        raise ValueError("--label-path must contain columns: site_id, label")

    labels["site_id"] = labels["site_id"].astype(str)
    labels["label"] = labels["label"].astype(int)

    eval_space = set(labels["site_id"])
    positives = set(labels.loc[labels["label"] == 1, "site_id"])
    universe = load_site_set(args.universe_path)

    print("Primary eval-space sites:", len(eval_space))
    print("Primary reference positives:", len(positives))
    print("Full common universe sites:", len(universe))

    models, alias_unique = load_gtf_models(args.gtf)

    summary_rows = []
    sanity_rows = []

    for fraction in fractions:
        print("\n" + "#" * 100)
        print(f"FRACTION: {fraction}")

        for tool in tools:
            path = raw_paths.get((fraction, tool))
            if path is None or not path.exists():
                if args.allow_missing:
                    print(f"WARNING: skip missing {fraction} {tool}")
                    continue
                raise FileNotFoundError(f"Missing {fraction} {tool}: {path}")

            print("\n" + "=" * 90)
            print(f"Standardizing {fraction} {tool}")
            print(path)

            std = standardize_one(tool, fraction, path, models, alias_unique, args)

            std["in_full_common_universe"] = std["site_id"].isin(universe)
            std["in_primary_eval_space"] = std["site_id"].isin(eval_space)
            std["is_primary_reference_positive"] = std["site_id"].isin(positives)

            score_path = score_call_dir / f"{fraction}.{tool}.genome_DRACH.score_call.tsv"
            std.to_csv(score_path, sep="\t", index=False)

            calls_eval = set(std.loc[to_bool(std["call"]) & std["in_primary_eval_space"], "site_id"].astype(str))
            calls_universe = set(std.loc[to_bool(std["call"]) & std["in_full_common_universe"], "site_id"].astype(str))

            pos_path = pos_call_dir / f"{fraction}.{tool}.primary_eval_positive_call_sites.tsv"
            pd.DataFrame({"site_id": sorted(calls_eval)}).to_csv(pos_path, sep="\t", index=False)

            tp = len(calls_eval & positives)
            n_called = len(calls_eval)
            non_ref = len(calls_eval - positives)

            print(f"n standardized genome sites: {std['site_id'].nunique()}")
            print(f"n overlap full common universe: {std['in_full_common_universe'].sum()}")
            print(f"n overlap primary eval-space: {std['in_primary_eval_space'].sum()}")
            print(f"primary eval positive calls: {n_called}")
            print(f"primary TP: {tp}")
            print(f"primary non-reference calls: {non_ref}")
            print(f"score/call table: {score_path}")
            print(f"positive call list: {pos_path}")

            summary_rows.append({
                "fraction": fraction,
                "tool": tool,
                "raw_path": str(path),
                "score_call_path": str(score_path),
                "positive_call_list_path": str(pos_path),
                "n_standardized_genome_sites": std["site_id"].nunique(),
                "n_rows_overlap_full_common_universe": int(std["in_full_common_universe"].sum()),
                "n_rows_overlap_primary_eval_space": int(std["in_primary_eval_space"].sum()),
                "n_primary_eval_positive_calls": n_called,
                "TP_primary": tp,
                "non_reference_calls_primary": non_ref,
            })

            if fraction == "full_depth" and tool in expected_full_depth:
                exp = expected_full_depth[tool]
                status = "MATCH" if (n_called == exp["n_called"] and tp == exp["TP"]) else "MISMATCH"
                minor_delta_1_nonref = (
                    tool == "CHEUI"
                    and tp == exp["TP"]
                    and abs(n_called - exp["n_called"]) == 1
                )

                sanity_rows.append({
                    "tool": tool,
                    "n_primary_eval_positive_calls": n_called,
                    "expected_n_primary_eval_positive_calls": exp["n_called"],
                    "TP_primary": tp,
                    "expected_TP_primary": exp["TP"],
                    "delta_n_called": n_called - exp["n_called"],
                    "delta_TP": tp - exp["TP"],
                    "sanity": status,
                    "minor_delta_1_nonref": minor_delta_1_nonref,
                })

                print(f"expected primary eval positive calls: {exp['n_called']}")
                print(f"expected primary TP: {exp['TP']}")
                print(f"full-depth sanity: {status}")
                if minor_delta_1_nonref:
                    print("note: CHEUI differs by only 1 non-reference call; TP is identical.")

    summary = pd.DataFrame(summary_rows)
    summary_path = qc_dir / "standardization_summary.primary_eval_space.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    print("\n" + "#" * 100)
    print("Saved summary:")
    print(summary_path)

    if sanity_rows:
        sanity = pd.DataFrame(sanity_rows)
        sanity_path = qc_dir / "full_depth_sanity_check.primary_eval_space.tsv"
        sanity.to_csv(sanity_path, sep="\t", index=False)

        print("Saved full-depth sanity check:")
        print(sanity_path)

        print("\nFull-depth sanity summary:")
        print(sanity.to_string(index=False))

        hard_fail = sanity[(sanity["sanity"] != "MATCH") & (~sanity["minor_delta_1_nonref"])]
        if len(hard_fail) > 0:
            print("\nWARNING: Some full-depth tools have hard mismatch.")
            print("Do not draw final figures until these are resolved.")
            print(hard_fail.to_string(index=False))
        else:
            print("\nOK: Full-depth sanity check passed, except possible tolerated CHEUI ±1 non-reference call.")
    else:
        print("No expected full-depth sanity table was provided; skipped expected-count sanity check.")

    print("\nDone.")
    print("Output directory:", args.outdir)
    print("Score/call tables:", score_call_dir)
    print("Positive call site lists:", pos_call_dir)
    print("QC:", qc_dir)


if __name__ == "__main__":
    main()
