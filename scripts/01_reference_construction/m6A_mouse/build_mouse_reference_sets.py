#!/usr/bin/env python3

import argparse
import os
from collections import Counter


def read_bed3(path):
    sites = set()
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("track") or line.startswith("#") or line.startswith("browser"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            sites.add((chrom, start, end))

    return sites


def write_bed3(path, sites):
    with open(path, "w") as f:
        for chrom, start, end in sorted(sites, key=lambda x: (x[0], x[1], x[2])):
            f.write(f"{chrom}\t{start}\t{end}\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build mouse m6A composite reference sets from m6A-CLIP, "
            "m6Aboost, and GLORI-derived positive sites."
        )
    )

    parser.add_argument(
        "--truth-shixiong",
        required=True,
        help="m6A-CLIP / external truth BED file.",
    )

    parser.add_argument(
        "--truth-m6aboost",
        required=True,
        help="m6Aboost predicted m6A sites BED file.",
    )

    parser.add_argument(
        "--truth-glori02",
        required=True,
        help="GLORI loose positive BED file, e.g. meanRatio >= 0.2.",
    )

    parser.add_argument(
        "--truth-glori095",
        required=True,
        help="GLORI strict positive BED file, e.g. meanRatio >= 0.95.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for reference sets.",
    )

    args = parser.parse_args()

    truth_clip = args.truth_clip
    truth_m6aboost = args.truth_m6aboost
    truth_glori02 = args.truth_glori02
    truth_glori095 = args.truth_glori095
    outdir = args.outdir

    os.makedirs(outdir, exist_ok=True)

    S = read_bed3(truth_shixiong)
    M = read_bed3(truth_m6aboost)
    G02 = read_bed3(truth_glori02)
    G095 = read_bed3(truth_glori095)

    print("Counts:")
    print("Shixiong =", len(S))
    print("m6Aboost =", len(M))
    print("GLORI0.2 =", len(G02))
    print("GLORI0.95 =", len(G095))

    strict_intersection_3way = S & M & G095
    strict_union_3way = S | M | G095

    counter_strict = Counter()
    for x in S:
        counter_strict[x] += 1
    for x in M:
        counter_strict[x] += 1
    for x in G095:
        counter_strict[x] += 1

    strict_consensus_2of3 = {k for k, v in counter_strict.items() if v >= 2}

    loose_intersection_3way = S & M & G02
    loose_union_3way = S | M | G02

    counter_loose = Counter()
    for x in S:
        counter_loose[x] += 1
    for x in M:
        counter_loose[x] += 1
    for x in G02:
        counter_loose[x] += 1

    loose_consensus_2of3 = {k for k, v in counter_loose.items() if v >= 2}

    ref_sets = {
        "strict_intersection_3way": strict_intersection_3way,
        "strict_consensus_2of3": strict_consensus_2of3,
        "strict_union_3way": strict_union_3way,
        "loose_intersection_3way": loose_intersection_3way,
        "loose_consensus_2of3": loose_consensus_2of3,
        "loose_union_3way": loose_union_3way,
    }

    for name, sites in ref_sets.items():
        write_bed3(f"{outdir}/{name}.bed", sites)

    summary_path = f"{outdir}/reference_set_counts.tsv"
    with open(summary_path, "w") as f:
        f.write("set_name\tsite_count\n")
        for name, sites in ref_sets.items():
            f.write(f"{name}\t{len(sites)}\n")

    pairwise = [
        ("Shixiong", "m6Aboost", len(S & M)),
        ("Shixiong", "GLORI0.2", len(S & G02)),
        ("Shixiong", "GLORI0.95", len(S & G095)),
        ("m6Aboost", "GLORI0.2", len(M & G02)),
        ("m6Aboost", "GLORI0.95", len(M & G095)),
        ("GLORI0.2", "GLORI0.95", len(G02 & G095)),
    ]

    pairwise_path = f"{outdir}/pairwise_overlap_counts.tsv"
    with open(pairwise_path, "w") as f:
        f.write("set1\tset2\toverlap_count\n")
        for a, b, n in pairwise:
            f.write(f"{a}\t{b}\t{n}\n")

    print("\nPairwise overlaps:")
    for a, b, n in pairwise:
        print(f"{a} ∩ {b} = {n}")

    print("\nReference sets:")
    for name, sites in ref_sets.items():
        print(f"{name} = {len(sites)}")

    print("\nSaved to:", outdir)
    print("Saved summary:", summary_path)
    print("Saved pairwise:", pairwise_path)


if __name__ == "__main__":
    main()