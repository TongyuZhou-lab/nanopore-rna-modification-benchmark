#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict


STD_CHRS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}


def open_maybe_gz(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_attr(attr):
    d = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr):
        d[m.group(1)] = m.group(2)
    return d


def load_gtf(gtf):
    exons = defaultdict(list)
    strand = {}
    chrom_of = {}

    print("Loading GTF:", gtf, file=sys.stderr)

    with open_maybe_gz(gtf) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t", 8)
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, st, frame, attr = parts

            if feature != "exon":
                continue
            if chrom not in STD_CHRS:
                continue

            a = parse_attr(attr)
            tid = a.get("transcript_id")
            if not tid:
                continue

            exons[tid].append((int(start), int(end)))
            strand[tid] = st
            chrom_of[tid] = chrom

    tx_info = {}

    for tid, xs in exons.items():
        st = strand[tid]
        chrom = chrom_of[tid]

        if st == "+":
            xs_sorted = sorted(xs, key=lambda x: x[0])
        else:
            xs_sorted = sorted(xs, key=lambda x: x[0], reverse=True)

        cum = []
        total = 0

        for s, e in xs_sorted:
            cum.append(total)
            total += e - s + 1

        tx_info[tid] = {
            "chrom": chrom,
            "strand": st,
            "exons": xs_sorted,
            "cum": cum,
            "length": total,
        }

    print("Loaded transcripts:", len(tx_info), file=sys.stderr)
    return tx_info


def normalize_tid(contig):
    return str(contig).split("|")[0]


def txpos_to_genome(tx_info, contig, tx_pos_1based):
    tid = normalize_tid(contig)

    if tid not in tx_info:
        return None

    info = tx_info[tid]
    pos = int(tx_pos_1based)

    if pos < 1 or pos > info["length"]:
        return None

    for (s, e), c in zip(info["exons"], info["cum"]):
        length = e - s + 1

        if c < pos <= c + length:
            offset = pos - c - 1

            if info["strand"] == "+":
                gpos = s + offset
            else:
                gpos = e - offset

            return info["chrom"], gpos

    return None


def is_center_c(seq):
    seq = str(seq).upper().replace("U", "T")

    if len(seq) == 9:
        return seq[4] == "C"
    if len(seq) == 5:
        return seq[2] == "C"

    return False


def bed_key(chrom, gpos_1based):
    start = int(gpos_1based) - 1
    end = int(gpos_1based)

    if start < 0:
        return None

    return chrom, start, end


def key_sorter(key):
    chrom, start, end = key

    def chrom_rank(c):
        if c.startswith("chr"):
            x = c[3:]
        else:
            x = c

        if x.isdigit():
            return (0, int(x))
        if x == "X":
            return (1, 23)
        if x == "Y":
            return (1, 24)
        if x in {"M", "MT"}:
            return (1, 25)

        return (2, x)

    return chrom_rank(chrom), start, end


def update_best(best, key, score):
    if key is None:
        return

    if key not in best or score > best[key]:
        best[key] = score


def save_scorebed(scores, outpath):
    with open(outpath, "w") as out:
        for chrom, start, end in sorted(scores.keys(), key=key_sorter):
            out.write(f"{chrom}\t{start}\t{end}\t{scores[(chrom, start, end)]:.10g}\n")

    print("saved:", outpath, "rows =", len(scores), file=sys.stderr)


def build_cheui_or_swarm(name, path, tx_info, out_prefix, outdir):
    print(f"\nBuilding {name}: {path}", file=sys.stderr)

    best_center_c = {}

    n = 0
    ok = 0
    missing = 0
    bad = 0
    center_c_rows = 0

    with open_maybe_gz(path) as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            n += 1

            try:
                contig = row["contig"]
                position = int(float(row["position"]))
                site = row.get("site", "")
                score = float(row["probability"])

                tx_center = position + 5

                conv = txpos_to_genome(tx_info, contig, tx_center)
                if conv is None:
                    missing += 1
                    continue

                chrom, gpos = conv
                key = bed_key(chrom, gpos)

                if is_center_c(site):
                    update_best(best_center_c, key, score)
                    center_c_rows += 1

                ok += 1

            except Exception:
                bad += 1

    print(
        f"{name}: total={n} ok={ok} missing={missing} bad={bad} centerC_rows={center_c_rows}",
        file=sys.stderr,
    )

    outpath = os.path.join(outdir, f"{out_prefix}.centerC.score.bed")
    save_scorebed(best_center_c, outpath)

    return {
        "tool": name,
        "input": path,
        "output": outpath,
        "raw_rows": n,
        "ok_rows": ok,
        "missing_rows": missing,
        "bad_rows": bad,
        "centerC_rows": center_c_rows,
        "unique_centerC_sites": len(best_center_c),
        "score_definition": "maximum site-level probability after transcript-to-genome mapping",
    }


def build_tandemmod_p95_ratio(path, outdir, assembly_label, p95_cutoff):
    print(f"\nBuilding TandemMod p95 ratio: {path}", file=sys.stderr)

    total = defaultdict(int)
    p95_count = defaultdict(int)

    n = 0
    bad = 0
    nonstd = 0
    noncenter = 0

    with open_maybe_gz(path) as f:
        for line in f:
            n += 1

            if not line.strip() or line.startswith("#"):
                continue

            fs = line.rstrip("\n").split("\t")

            if len(fs) < 7:
                bad += 1
                continue

            try:
                chrom = fs[2]
                gpos = int(float(fs[3]))
                motif = fs[4]
                prob = float(fs[6])
            except Exception:
                bad += 1
                continue

            if chrom not in STD_CHRS:
                nonstd += 1
                continue

            if not is_center_c(motif):
                noncenter += 1
                continue

            key = bed_key(chrom, gpos)
            if key is None:
                bad += 1
                continue

            total[key] += 1

            if prob > p95_cutoff:
                p95_count[key] += 1

            if n % 10000000 == 0:
                print("TandemMod processed rows:", n, "sites:", len(total), file=sys.stderr)

    ratio_scores = {}
    detail_rows = []

    for key in total:
        ratio = p95_count[key] / total[key] if total[key] > 0 else 0.0
        ratio_scores[key] = ratio
        chrom, start, end = key
        detail_rows.append((chrom, start, end, total[key], p95_count[key], ratio))

    score_out = os.path.join(outdir, f"TandemMod_human_m5C.{assembly_label}.centerC.p95_ratio.score.bed")
    detail_out = os.path.join(outdir, f"TandemMod_human_m5C.{assembly_label}.centerC.p95_ratio.detail.tsv")

    save_scorebed(ratio_scores, score_out)

    with open(detail_out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["chrom", "start", "end", "total", "p95", "p95_ratio"])

        for chrom, start, end, total_count, p95, ratio in sorted(detail_rows, key=lambda x: key_sorter((x[0], x[1], x[2]))):
            writer.writerow([chrom, start, end, total_count, p95, f"{ratio:.10g}"])

    print("saved:", detail_out, "rows =", len(detail_rows), file=sys.stderr)

    return {
        "tool": "TandemMod",
        "input": path,
        "output": score_out,
        "detail_output": detail_out,
        "raw_rows": n,
        "bad_rows": bad,
        "nonstandard_rows": nonstd,
        "noncenter_rows": noncenter,
        "unique_centerC_sites": len(ratio_scores),
        "score_definition": f"count(probability > {p95_cutoff}) / total",
    }


def write_summary(rows, outpath):
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with open(outpath, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for row in rows:
            writer.writerow(row)

    print("saved:", outpath, "rows =", len(rows), file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build human m5C genome-level continuous score BED files for TandemMod, "
            "CHEUI, and SWARM. TandemMod is summarized as p95_ratio "
            "to align with the m6A score definition."
        )
    )

    parser.add_argument("--gtf", required=True, help="GENCODE GRCh38 GTF annotation file.")
    parser.add_argument("--tandemmod", required=True, help="TandemMod genome-level read-level prediction TSV.")
    parser.add_argument("--cheui", required=True, help="CHEUI site-level m5C prediction TSV.")
    parser.add_argument("--swarm", required=True, help="SWARM m5C prediction TSV.")
    parser.add_argument("--outdir", required=True, help="Output directory for score BED files.")
    parser.add_argument("--assembly-label", default="GRCh38", help="Assembly label used in output filenames.")
    parser.add_argument("--tandemmod-p95-cutoff", type=float, default=0.95, help="Read-level probability cutoff for TandemMod p95_ratio.")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    tx_info = load_gtf(args.gtf)

    summaries = []

    summaries.append(
        build_tandemmod_p95_ratio(
            path=args.tandemmod,
            outdir=args.outdir,
            assembly_label=args.assembly_label,
            p95_cutoff=args.tandemmod_p95_cutoff,
        )
    )

    summaries.append(
        build_cheui_or_swarm(
            name="CHEUI",
            path=args.cheui,
            tx_info=tx_info,
            out_prefix=f"CHEUI_human_m5C.{args.assembly_label}.centerPlus5",
            outdir=args.outdir,
        )
    )

    summaries.append(
        build_cheui_or_swarm(
            name="SWARM",
            path=args.swarm,
            tx_info=tx_info,
            out_prefix=f"SWARM_human_m5C.{args.assembly_label}.centerPlus5",
            outdir=args.outdir,
        )
    )

    summary_out = os.path.join(args.outdir, f"human_m5C.{args.assembly_label}.scorebed_standardization.summary.tsv")
    write_summary(summaries, summary_out)

    print("\nAll done.", file=sys.stderr)


if __name__ == "__main__":
    main()
