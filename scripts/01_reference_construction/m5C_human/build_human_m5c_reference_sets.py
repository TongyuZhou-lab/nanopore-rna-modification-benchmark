#!/usr/bin/env python3
import argparse
import os
from collections import defaultdict


def read_bed(path):
    s = set()
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue

            fs = line.rstrip("\n").split("\t")
            if len(fs) < 3:
                continue

            chrom = fs[0]
            start = int(fs[1])
            end = int(fs[2])
            s.add((chrom, start, end))

    return s


def write_bed(keys, path):
    keys = sorted(keys, key=lambda x: (x[0], x[1], x[2]))

    with open(path, "w") as out:
        for chrom, start, end in keys:
            out.write(f"{chrom}\t{start}\t{end}\n")

    print("saved:", path, "n =", len(keys))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build human m5C reference sets from three HeLa m5C positive reference sources."
        )
    )

    parser.add_argument(
        "--gse122260",
        required=True,
        help="GSE122260 Ctrl GRCh38 BED file.",
    )

    parser.add_argument(
        "--gse140995",
        required=True,
        help="GSE140995 transcriptome-wide GRCh38 BED file.",
    )

    parser.add_argument(
        "--gse93749",
        required=True,
        help="GSE93749 HeLa siCTRL final GRCh38 BED file.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for human m5C reference sets.",
    )

    args = parser.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)


    sources = {
        "GSE122260_Ctrl": args.gse122260,
        "GSE140995": args.gse140995,
        "GSE93749": args.gse93749,
    }

    sets = {}

    for name, path in sources.items():
        if not os.path.exists(path):
            raise SystemExit(f"[ERROR] missing source: {name} {path}")

        sets[name] = read_bed(path)
        print(name, len(sets[name]))

    support = defaultdict(int)
    support_names = defaultdict(list)

    for name, s in sets.items():
        for k in s:
            support[k] += 1
            support_names[k].append(name)

    intersection_3way = {k for k, n in support.items() if n == 3}
    consensus_2of3 = {k for k, n in support.items() if n >= 2}
    union_3way = {k for k, n in support.items() if n >= 1}


    write_bed(intersection_3way, os.path.join(outdir, "strict_intersection_3way.bed"))
    write_bed(consensus_2of3, os.path.join(outdir, "strict_consensus_2of3.bed"))
    write_bed(union_3way, os.path.join(outdir, "loose_union_3way.bed"))

    write_bed(union_3way, os.path.join(outdir, "strict_union_3way.bed"))
    write_bed(intersection_3way, os.path.join(outdir, "loose_intersection_3way.bed"))
    write_bed(consensus_2of3, os.path.join(outdir, "loose_consensus_2of3.bed"))

    with open(os.path.join(outdir, "reference_set_counts.tsv"), "w") as out:
        out.write("set\tn\n")
        for name, s in sets.items():
            out.write(f"{name}\t{len(s)}\n")

        out.write(f"strict_intersection_3way\t{len(intersection_3way)}\n")
        out.write(f"strict_consensus_2of3\t{len(consensus_2of3)}\n")
        out.write(f"strict_union_3way\t{len(union_3way)}\n")
        out.write(f"loose_intersection_3way\t{len(intersection_3way)}\n")
        out.write(f"loose_consensus_2of3\t{len(consensus_2of3)}\n")
        out.write(f"loose_union_3way\t{len(union_3way)}\n")

    names = list(sets.keys())

    with open(os.path.join(outdir, "pairwise_overlap_counts.tsv"), "w") as out:
        out.write("source1\tsource2\tn1\tn2\toverlap\n")
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = names[i], names[j]
                out.write(
                    f"{a}\t{b}\t"
                    f"{len(sets[a])}\t{len(sets[b])}\t"
                    f"{len(sets[a] & sets[b])}\n"
                )

    print("summary:", os.path.join(outdir, "reference_set_counts.tsv"))
    print("pairwise:", os.path.join(outdir, "pairwise_overlap_counts.tsv"))


if __name__ == "__main__":
    main()