#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import re
from collections import defaultdict


STD_CHRS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}


def open_maybe_gz(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def parse_transcript_id(attrs):
    m = re.search(r'transcript_id "([^"]+)"', attrs)
    return m.group(1) if m else None


def load_gtf_exons(gtf_path, keep_standard_chromosomes=True):
    exon_dict = {}
    chr_dict = {}
    strand_dict = {}

    with open_maybe_gz(gtf_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            if feature != "exon":
                continue
            if keep_standard_chromosomes and chrom not in STD_CHRS:
                continue

            tx_id = parse_transcript_id(attrs)
            if tx_id is None:
                continue

            exon_dict.setdefault(tx_id, []).append((int(start), int(end)))
            chr_dict[tx_id] = chrom
            strand_dict[tx_id] = strand

    for tx in exon_dict:
        exon_dict[tx] = sorted(exon_dict[tx], key=lambda x: x[0])

    # Add no-version aliases to match tool outputs such as ENSTxxx and ENSTxxx.10.
    for tx in list(exon_dict.keys()):
        tx_nover = tx.split(".")[0]
        if tx_nover not in exon_dict:
            exon_dict[tx_nover] = exon_dict[tx]
            chr_dict[tx_nover] = chr_dict[tx]
            strand_dict[tx_nover] = strand_dict[tx]

    print("Loaded transcript aliases:", len(exon_dict))
    return exon_dict, chr_dict, strand_dict


def txpos_to_genome(tx_id, tx_pos_1based, exon_dict, chr_dict, strand_dict):
    tx_key = tx_id if tx_id in exon_dict else tx_id.split(".")[0]
    if tx_key not in exon_dict:
        return None

    exons = exon_dict[tx_key]
    chrom = chr_dict[tx_key]
    strand = strand_dict[tx_key]

    if strand == "+":
        acc = 0
        for start, end in exons:
            exon_len = end - start + 1
            if tx_pos_1based > acc + exon_len:
                acc += exon_len
            else:
                genome_pos = start + (tx_pos_1based - acc) - 1
                return chrom, genome_pos
    else:
        acc = 0
        for start, end in exons[::-1]:
            exon_len = end - start + 1
            if tx_pos_1based > acc + exon_len:
                acc += exon_len
            else:
                genome_pos = end - (tx_pos_1based - acc) + 1
                return chrom, genome_pos

    return None


def key_sorter(key):
    chrom, start, end = key

    def chrom_rank(c):
        x = c[3:] if c.startswith("chr") else c
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


def write_score_bed(out_path, score_map):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as out:
        for chrom, start, end in sorted(score_map.keys(), key=key_sorter):
            out.write(f"{chrom}\t{start}\t{end}\t{score_map[(chrom, start, end)]:.10g}\n")
    print("Saved:", out_path)
    print("unique sites =", len(score_map))


def write_summary(out_path, rows):
    rows = [x for x in rows if x is not None]
    if not rows:
        return

    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with open(out_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print("Saved:", out_path)


def convert_transcript_site_tsv(input_path, out_path, tool_name, exon_dict, chr_dict, strand_dict, center_offset):
    """
    CHEUI/SWARM-style input:
    contig, position, site, coverage, stoichiometry, probability.

    Position is treated as the k-mer start, and the modified base is
    position + center_offset.
    """
    if not input_path:
        print(f"[SKIP] {tool_name}: input path not provided")
        return None
    if not os.path.exists(input_path):
        print(f"[SKIP] {tool_name}: input not found:", input_path)
        return None

    score_map = {}
    total = converted = missing_tx = bad = 0

    with open_maybe_gz(input_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            total += 1
            try:
                tx_id = row["contig"].split("|")[0]
                kmer_start = int(float(row["position"]))
                tx_center_pos = kmer_start + center_offset
                score = float(row["probability"])
            except Exception:
                bad += 1
                continue

            res = txpos_to_genome(tx_id, tx_center_pos, exon_dict, chr_dict, strand_dict)
            if res is None:
                missing_tx += 1
                continue

            chrom, genome_pos_1based = res
            start = genome_pos_1based - 1
            end = genome_pos_1based

            if start < 0:
                bad += 1
                continue

            key = (chrom, start, end)
            if key not in score_map or score > score_map[key]:
                score_map[key] = score
            converted += 1

    print(f"\n=== {tool_name} ===")
    print("input:", input_path)
    print("total rows:", total)
    print("converted rows:", converted)
    print("missing transcript:", missing_tx)
    print("bad rows:", bad)
    write_score_bed(out_path, score_map)

    return {
        "tool": tool_name,
        "input": input_path,
        "output": out_path,
        "raw_rows": total,
        "converted_rows": converted,
        "missing_transcript_rows": missing_tx,
        "bad_rows": bad,
        "unique_sites": len(score_map),
        "score_definition": "maximum probability after transcript-to-genome mapping",
    }


def convert_m6anet(input_path, out_path, exon_dict, chr_dict, strand_dict):
    if not input_path:
        print("[SKIP] m6Anet: input path not provided")
        return None
    if not os.path.exists(input_path):
        print("[SKIP] m6Anet input not found:", input_path)
        return None

    score_map = {}
    total = converted = missing_tx = bad = 0

    with open_maybe_gz(input_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            total += 1
            try:
                tx_id = row["transcript_id"].split("|")[0]
                tx_pos0 = int(float(row["transcript_position"]))
                tx_pos1 = tx_pos0 + 1
                score = float(row["probability_modified"])
            except Exception:
                bad += 1
                continue

            res = txpos_to_genome(tx_id, tx_pos1, exon_dict, chr_dict, strand_dict)
            if res is None:
                missing_tx += 1
                continue

            chrom, genome_pos_1based = res
            key = (chrom, genome_pos_1based - 1, genome_pos_1based)
            if key not in score_map or score > score_map[key]:
                score_map[key] = score
            converted += 1

    print("\n=== m6Anet_human_DRACH ===")
    print("input:", input_path)
    print("total rows:", total)
    print("converted rows:", converted)
    print("missing transcript:", missing_tx)
    print("bad rows:", bad)
    write_score_bed(out_path, score_map)

    return {
        "tool": "m6Anet_human_DRACH",
        "input": input_path,
        "output": out_path,
        "raw_rows": total,
        "converted_rows": converted,
        "missing_transcript_rows": missing_tx,
        "bad_rows": bad,
        "unique_sites": len(score_map),
        "score_definition": "maximum probability_modified using transcript_position + 1",
    }


def convert_tandemmod(input_path, out_path, exon_dict, chr_dict, strand_dict, center_offset):
    """
    Supports two TandemMod formats:
    1. Genome-level table with chr, genome_pos/site, p_0.95 and total.
       Score is p_0.95 / total.
    2. Transcript-level table with contig, position and probability.
       Score is maximum probability after transcript-to-genome mapping.
    """
    if not input_path:
        print("[SKIP] TandemMod: input path not provided")
        return None
    if not os.path.exists(input_path):
        print("[SKIP] TandemMod input not found:", input_path)
        return None

    with open_maybe_gz(input_path) as f:
        header = f.readline().rstrip("\n").split("\t")

    idx = {x: i for i, x in enumerate(header)}
    score_map = {}

    if "chr" in idx and ("p_0.95" in idx) and ("total" in idx):
        chr_i = idx["chr"]
        if "genome_pos" in idx:
            pos_i = idx["genome_pos"]
        else:
            site_idxs = [i for i, col in enumerate(header) if col == "site"]
            if site_idxs:
                pos_i = site_idxs[-1]
            else:
                raise ValueError(f"TandemMod genome-level table lacks genome_pos/site column: {header}")

        p_i = idx["p_0.95"]
        t_i = idx["total"]

        total_rows = bad = zero_total = 0
        with open_maybe_gz(input_path) as f:
            next(f)
            for line in f:
                if not line.strip():
                    continue
                total_rows += 1
                fs = line.rstrip("\n").split("\t")
                try:
                    chrom = fs[chr_i]
                    genome_pos = int(float(fs[pos_i]))
                    p95 = float(fs[p_i])
                    total = float(fs[t_i])
                    if total <= 0:
                        zero_total += 1
                        continue
                    score = p95 / total
                except Exception:
                    bad += 1
                    continue

                key = (chrom, genome_pos - 1, genome_pos)
                if key not in score_map or score > score_map[key]:
                    score_map[key] = score

        print("\n=== TandemMod_human_DRACH genome-level ===")
        print("input:", input_path)
        print("total rows:", total_rows)
        print("bad rows:", bad)
        print("zero-total rows:", zero_total)
        write_score_bed(out_path, score_map)

        return {
            "tool": "TandemMod_human_DRACH",
            "input": input_path,
            "output": out_path,
            "raw_rows": total_rows,
            "bad_rows": bad,
            "zero_total_rows": zero_total,
            "unique_sites": len(score_map),
            "score_definition": "p_0.95 / total",
        }

    if all(x in idx for x in ["contig", "position", "probability"]):
        return convert_transcript_site_tsv(
            input_path,
            out_path,
            "TandemMod_human_DRACH transcript-level",
            exon_dict,
            chr_dict,
            strand_dict,
            center_offset,
        )

    raise ValueError(f"Unrecognized TandemMod format. Header={header}")


def main():
    parser = argparse.ArgumentParser(
        description="Build human m6A genome-level score BED files from TandemMod, m6Anet, CHEUI, and SWARM outputs."
    )

    parser.add_argument("--gtf", required=True, help="GENCODE GRCh38 GTF annotation file.")
    parser.add_argument("--tandemmod", required=True, help="TandemMod site-level table.")
    parser.add_argument("--m6anet", required=True, help="m6Anet data.site_proba.csv.")
    parser.add_argument("--cheui-all", required=True, help="CHEUI all site-level prediction file.")
    parser.add_argument("--cheui-drach", required=True, help="CHEUI DRACH-only site-level prediction file.")
    parser.add_argument("--swarm-all", required=True, help="SWARM all prediction TSV.")
    parser.add_argument("--swarm-drach", required=True, help="SWARM DRACH-only prediction TSV.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--assembly-label", default="GRCh38", help="Assembly label in output filenames.")
    parser.add_argument("--center-offset", type=int, default=5, help="Offset from k-mer start to modified base for CHEUI/SWARM-style 9-mer outputs.")
    parser.add_argument("--keep-all-chromosomes", action="store_true", help="Do not restrict the GTF to chr1-22, chrX, chrY and chrM.")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    exon_dict, chr_dict, strand_dict = load_gtf_exons(
        args.gtf,
        keep_standard_chromosomes=not args.keep_all_chromosomes,
    )

    suffix = "score.bed" if args.keep_all_chromosomes else "stdchr.score.bed"

    summaries = [
        convert_tandemmod(
            args.tandemmod,
            os.path.join(args.outdir, f"TandemMod_human_DRACH.{args.assembly_label}.{suffix}"),
            exon_dict,
            chr_dict,
            strand_dict,
            args.center_offset,
        ),
        convert_m6anet(
            args.m6anet,
            os.path.join(args.outdir, f"m6Anet_human_DRACH.{args.assembly_label}.txPlus1.{suffix}"),
            exon_dict,
            chr_dict,
            strand_dict,
        ),
        convert_transcript_site_tsv(
            args.cheui_all,
            os.path.join(args.outdir, f"CHEUI_human_all.{args.assembly_label}.centerPlus5.{suffix}"),
            "CHEUI_human_all",
            exon_dict,
            chr_dict,
            strand_dict,
            args.center_offset,
        ),
        convert_transcript_site_tsv(
            args.cheui_drach,
            os.path.join(args.outdir, f"CHEUI_human_DRACH.{args.assembly_label}.centerPlus5.{suffix}"),
            "CHEUI_human_DRACH",
            exon_dict,
            chr_dict,
            strand_dict,
            args.center_offset,
        ),
        convert_transcript_site_tsv(
            args.swarm_all,
            os.path.join(args.outdir, f"SWARM_human_all.{args.assembly_label}.centerPlus5.{suffix}"),
            "SWARM_human_all",
            exon_dict,
            chr_dict,
            strand_dict,
            args.center_offset,
        ),
        convert_transcript_site_tsv(
            args.swarm_drach,
            os.path.join(args.outdir, f"SWARM_human_DRACH.{args.assembly_label}.centerPlus5.{suffix}"),
            "SWARM_human_DRACH",
            exon_dict,
            chr_dict,
            strand_dict,
            args.center_offset,
        ),
    ]

    summary_out = os.path.join(args.outdir, f"human_m6A.{args.assembly_label}.scorebed_standardization.summary.tsv")
    write_summary(summary_out, summaries)


if __name__ == "__main__":
    main()
