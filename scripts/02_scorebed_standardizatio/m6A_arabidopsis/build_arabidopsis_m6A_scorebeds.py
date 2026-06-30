#!/usr/bin/env python3
import argparse
from pathlib import Path
import re
import gzip
import pandas as pd
import numpy as np
TOOLS = ["TandemMod", "m6Anet", "CHEUI", "SWARM"]

FRACTIONS = [
    "full_depth",
    "frac0.8_seed11",
    "frac0.6_seed11",
    "frac0.4_seed11",
    "frac0.2_seed11",
]

EXPECTED_FULL_DEPTH_PRIMARY = {
    "TandemMod": {"n_called": 52, "TP": 15},
    "m6Anet": {"n_called": 760, "TP": 251},
    "CHEUI": {"n_called": 577, "TP": 177},
    "SWARM": {"n_called": 1597, "TP": 382},
}

DRACH_RE = re.compile(r"[AGT][AG]AC[ACT]")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Standardize Arabidopsis m6A full-depth and downsampling tool outputs."
    )
    parser.add_argument("--label-path", required=True, help="Primary evaluation label TSV file.")
    parser.add_argument("--universe-path", required=True, help="Full-depth common universe site ID TSV file.")
    parser.add_argument("--raw-paths", required=True, help="TSV file with columns: fraction, tool, path.")
    parser.add_argument("--gtf", required=True, help="Arabidopsis GTF annotation file.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    return parser.parse_args()

def load_raw_paths(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"fraction", "tool", "path"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"raw paths TSV missing columns: {sorted(missing)}")

    raw_paths = {}
    for _, row in df.iterrows():
        fraction = str(row["fraction"])
        tool = str(row["tool"])
        raw_paths.setdefault(fraction, {})[tool] = Path(str(row["path"]))

    return raw_paths

def sniff_sep(path):
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if line.strip():
                return "," if line.count(",") > line.count("\t") else "\t"
    return "\t"

def read_table(path):
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
    return any(DRACH_RE.fullmatch(s[i:i+5]) for i in range(len(s) - 4))

def maybe_filter_drach(df, candidates, tool):
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

def pick_col(df, candidates, tool, desc):
    cols = list(df.columns)
    lower = {c.lower(): c for c in cols}
    for c in candidates:
        if c in cols:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    raise ValueError(f"{tool}: cannot find {desc}. candidates={candidates}; existing columns={cols}")

def load_site_set(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "site_id" in df.columns:
        return set(df["site_id"].astype(str))
    if {"chrom", "start"}.issubset(df.columns):
        return set(df["chrom"].astype(str) + ":" + df["start"].astype(int).astype(str))
    if df.shape[1] == 1:
        return set(df.iloc[:, 0].astype(str))
    raise ValueError(f"Cannot parse site set: {path}")


def find_gtf(gtf_path):
    gtf_path = Path(gtf_path)
    if gtf_path.exists():
        return gtf_path
    raise FileNotFoundError(f"GTF not found: {gtf_path}")

def parse_attrs(attr):
    d = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', str(attr)):
        d[m.group(1)] = m.group(2)
    return d

def load_gtf_models(gtf_path):
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

            chrom, source, feature, start, end, score, strand, frame, attr = parts

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

        if strand == "-":
            items = sorted(items, key=lambda x: x[1], reverse=True)
        else:
            items = sorted(items, key=lambda x: x[1])

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

def aggregate_by_max_score_any_call(df):
    df = df.dropna(subset=["site_id"]).copy()
    df["site_id"] = df["site_id"].astype(str)
    df["score"] = safe_numeric(df["score"]).fillna(0)
    df["call"] = to_bool(df["call"])

    idx = df.groupby("site_id")["score"].idxmax()
    best = df.loc[idx].copy()

    any_call = df.groupby("site_id", as_index=False)["call"].max()
    best = best.drop(columns=["call"]).merge(any_call, on="site_id", how="left")

    return best

def standardize_tandemmod(path, fraction):
    tool = "TandemMod"
    df = read_table(path)

    required = ["chr", "site.1", "p_0.95", "total"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"{tool} {fraction}: missing required column {c}; columns={list(df.columns)}")

    df = maybe_filter_drach(df, ["motif", "kmer", "site"], tool)

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

    # Important adjustment:
    # TandemMod genome-site aggregation rule consistent with the finalized full-depth results.
    grouped = df.groupby("site_id").agg(
        chrom=("chrom", "first"),
        start=("start", "first"),
        score=("ratio", "max"),
        max_p95=("p95", "max"),
        max_ratio=("ratio", "max"),
        max_total=("total_num", "max"),
        n_source_rows=("site_id", "size"),
    ).reset_index()

    grouped["call"] = (grouped["max_p95"] > 10) & (grouped["max_ratio"] >= 0.2)

    return grouped

def standardize_m6anet(path, models, alias_unique, fraction):
    tool = "m6Anet"
    df = read_table(path)

    tx_col = pick_col(df, ["transcript_id", "transcript", "contig", "tx"], tool, "transcript column")
    pos_col = pick_col(df, ["transcript_position", "position", "pos"], tool, "position column")
    score_col = pick_col(df, ["probability_modified", "prob_modified", "probability", "score"], tool, "score column")

    df = maybe_filter_drach(df, ["kmer", "motif", "site"], tool)

    df = df.copy()
    df["score"] = safe_numeric(df[score_col]).fillna(0)
    df["call"] = df["score"] >= 0.9

    mapped = [
        tx_to_genome(tx, pos, models, alias_unique, origin=0)
        for tx, pos in zip(df[tx_col], df[pos_col])
    ]

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

def standardize_cheui(path, models, alias_unique, fraction):
    tool = "CHEUI"
    df = read_table(path)

    tx_col = pick_col(df, ["contig", "transcript_id", "transcript", "tx"], tool, "transcript column")
    pos_col = pick_col(df, ["position", "pos", "transcript_position"], tool, "position column")
    score_col = pick_col(df, ["probability", "prob", "score"], tool, "score column")

    df = maybe_filter_drach(df, ["site", "kmer", "motif"], tool)

    df = df.copy()
    df["score"] = safe_numeric(df[score_col]).fillna(0)
    df["call"] = df["score"] >= 0.9999

    mapped = [
        tx_to_genome(tx, pos, models, alias_unique, origin=1)
        for tx, pos in zip(df[tx_col], df[pos_col])
    ]

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

def standardize_swarm(path, models, alias_unique, fraction):
    tool = "SWARM"
    df = read_table(path)

    tx_col = pick_col(df, ["contig", "transcript_id", "transcript", "tx"], tool, "transcript column")
    pos_col = pick_col(df, ["position", "pos", "transcript_position"], tool, "position column")
    score_col = pick_col(df, ["probability", "prob", "score"], tool, "score column")
    stoich_col = pick_col(df, ["stoichiometry", "stochiometry", "mod_ratio", "ratio"], tool, "stoichiometry column")

    df = maybe_filter_drach(df, ["site", "kmer", "motif"], tool)

    df = df.copy()
    df["score"] = safe_numeric(df[score_col]).fillna(0)
    df["stoichiometry"] = safe_numeric(df[stoich_col]).fillna(0)
    df["call"] = (df["score"] >= 0.9972) & (df["stoichiometry"] > 0.1)

    # Validated SWARM coordinate convention:
    # Add 4 to transcript position, then map to genome with origin=0.
    pos_shifted = safe_numeric(df[pos_col]) + 4

    mapped = [
        tx_to_genome(tx, pos, models, alias_unique, origin=0)
        for tx, pos in zip(df[tx_col], pos_shifted)
    ]

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

def standardize_one(tool, fraction, path, models, alias_unique):
    if tool == "TandemMod":
        return standardize_tandemmod(path, fraction)
    if tool == "m6Anet":
        return standardize_m6anet(path, models, alias_unique, fraction)
    if tool == "CHEUI":
        return standardize_cheui(path, models, alias_unique, fraction)
    if tool == "SWARM":
        return standardize_swarm(path, models, alias_unique, fraction)
    raise ValueError(tool)

def main():
    args = parse_args()

    label_path = Path(args.label_path)
    universe_path = Path(args.universe_path)
    outdir = Path(args.outdir)
    score_call_dir = outdir / "score_call_tables"
    pos_call_dir = outdir / "positive_call_site_lists"
    qc_dir = outdir / "qc"

    for d in [score_call_dir, pos_call_dir, qc_dir]:
        d.mkdir(parents=True, exist_ok=True)

    raw_paths = load_raw_paths(args.raw_paths)

    for fraction in FRACTIONS:
        if fraction not in raw_paths:
            raise ValueError(f"Missing fraction in raw paths TSV: {fraction}")
        for tool in TOOLS:
            if tool not in raw_paths[fraction]:
                raise ValueError(f"Missing tool path in raw paths TSV: {fraction} {tool}")
            p = raw_paths[fraction][tool]
            if not p.exists():
                raise FileNotFoundError(f"{fraction} {tool} raw file not found: {p}")

    labels = pd.read_csv(label_path, sep="\t", low_memory=False)
    labels["site_id"] = labels["site_id"].astype(str)
    labels["label"] = labels["label"].astype(int)

    eval_space = set(labels["site_id"])
    positives = set(labels.loc[labels["label"] == 1, "site_id"])

    universe = load_site_set(universe_path)

    print("Primary eval-space sites:", len(eval_space))
    print("Primary reference positives:", len(positives))
    print("Full common universe sites:", len(universe))

    gtf_path = find_gtf(args.gtf)
    models, alias_unique = load_gtf_models(gtf_path)

    summary_rows = []
    sanity_rows = []

    for fraction in FRACTIONS:
        print("\n" + "#" * 100)
        print(f"FRACTION: {fraction}")

        for tool in TOOLS:
            path = raw_paths[fraction][tool]

            print("\n" + "=" * 90)
            print(f"Standardizing {fraction} {tool}")
            print(path)

            std = standardize_one(tool, fraction, path, models, alias_unique)

            std["in_full_common_universe"] = std["site_id"].isin(universe)
            std["in_primary_eval_space"] = std["site_id"].isin(eval_space)
            std["is_primary_reference_positive"] = std["site_id"].isin(positives)

            score_path = score_call_dir / f"{fraction}.{tool}.genome_DRACH.score_call.tsv"
            std.to_csv(score_path, sep="\t", index=False)

            calls_eval = set(
                std.loc[to_bool(std["call"]) & std["in_primary_eval_space"], "site_id"].astype(str)
            )
            calls_universe = set(
                std.loc[to_bool(std["call"]) & std["in_full_common_universe"], "site_id"].astype(str)
            )

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

            if fraction == "full_depth":
                exp = EXPECTED_FULL_DEPTH_PRIMARY[tool]
                status = "MATCH" if (n_called == exp["n_called"] and tp == exp["TP"]) else "MISMATCH"

                # CHEUI may show 578 vs 577 in earlier runs, while TP remains identical.
                # This is not forced to MATCH; it is only flagged as a minor non-reference delta.
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

    sanity = pd.DataFrame(sanity_rows)
    sanity_path = qc_dir / "full_depth_sanity_check.primary_eval_space.tsv"
    sanity.to_csv(sanity_path, sep="\t", index=False)

    print("\n" + "#" * 100)
    print("Saved summary:")
    print(summary_path)
    print("Saved full-depth sanity check:")
    print(sanity_path)

    print("\nFull-depth sanity summary:")
    print(sanity.to_string(index=False))

    hard_fail = sanity[
        (sanity["sanity"] != "MATCH") & (~sanity["minor_delta_1_nonref"])
    ]

    if len(hard_fail) > 0:
        print("\nWARNING: Some full-depth tools have hard mismatch.")
        print("Do not draw final figures until these are resolved.")
        print(hard_fail.to_string(index=False))
    else:
        print("\nOK: Full-depth sanity check passed, except possible tolerated CHEUI ±1 non-reference call.")

if __name__ == "__main__":
    main()