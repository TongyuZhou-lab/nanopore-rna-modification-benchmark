#!/usr/bin/env python3
import argparse
import csv
import gzip
import os


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def mouse_std_chrs():
    return {f"chr{i}" for i in range(1, 20)} | {"chrX", "chrY", "chrM"}


def is_std_chr(chrom, stdchr_only=True):
    if not stdchr_only:
        return True
    return chrom in mouse_std_chrs()


def ensure_outdir(outdir):
    os.makedirs(outdir, exist_ok=True)


def load_gtf_exons(gtf_path):
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

            transcript_id = None
            for part in attrs.split(";"):
                part = part.strip()
                if part.startswith('transcript_id "'):
                    transcript_id = part.split('"')[1]
                    break

            if transcript_id is None:
                continue

            start = int(start)
            end = int(end)

            if transcript_id not in exon_dict:
                exon_dict[transcript_id] = []
                chr_dict[transcript_id] = chrom
                strand_dict[transcript_id] = strand

            exon_dict[transcript_id].append((start, end))

    for tx in exon_dict:
        exon_dict[tx] = sorted(exon_dict[tx], key=lambda x: x[0])

    return exon_dict, chr_dict, strand_dict


def txpos_to_genome(tx_id, tx_pos_1based, exon_dict, chr_dict, strand_dict):
    if tx_id not in exon_dict:
        return None

    exons = exon_dict[tx_id]
    strand = strand_dict[tx_id]
    chrom = chr_dict[tx_id]

    if strand == "+":
        acc = 0
        for start, end in exons:
            exon_len = end - start + 1
            if tx_pos_1based > acc + exon_len:
                acc += exon_len
            else:
                genome_pos = start + (tx_pos_1based - acc) - 1
                return chrom, genome_pos, strand
    else:
        acc = 0
        for start, end in exons[::-1]:
            exon_len = end - start + 1
            if tx_pos_1based > acc + exon_len:
                acc += exon_len
            else:
                genome_pos = end - (tx_pos_1based - acc) + 1
                return chrom, genome_pos, strand

    return None


def write_score_bed(path, rows):
    with open(path, "w") as f:
        for chrom, start, end, score in sorted(rows, key=lambda x: (x[0], x[1], x[2])):
            f.write(f"{chrom}\t{start}\t{end}\t{score}\n")

    print("saved:", path, "rows =", len(rows))


def build_tandemmod(tandemmod_in, outdir, assembly_label, stdchr_only=True):
    rows = []

    with open_maybe_gz(tandemmod_in) as f:
        header = next(f).rstrip("\n").split("\t")
        idx = {x: i for i, x in enumerate(header)}

        if "chr" not in idx:
            raise ValueError(f"TandemMod header has no 'chr' column: {header}")
        chr_i = idx["chr"]

        if "genome_pos" in idx:
            pos_i = idx["genome_pos"]
        else:
            site_idxs = [i for i, col in enumerate(header) if col == "site"]
            if len(site_idxs) >= 2:
                pos_i = site_idxs[-1]
            elif len(site_idxs) == 1:
                pos_i = site_idxs[0]
            else:
                raise ValueError(f"TandemMod header has neither 'genome_pos' nor 'site': {header}")

        if "p_0.95" not in idx or "total" not in idx:
            raise ValueError(f"TandemMod header missing p_0.95/total: {header}")

        p95_i = idx["p_0.95"]
        total_i = idx["total"]

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            chrom = fields[chr_i]

            if not is_std_chr(chrom, stdchr_only):
                continue

            genome_pos = int(float(fields[pos_i]))
            p95 = float(fields[p95_i])
            total = float(fields[total_i])

            if total <= 0:
                continue

            score = p95 / total
            rows.append((chrom, genome_pos - 1, genome_pos, score))

    suffix = ".stdchr" if stdchr_only else ""
    out = os.path.join(outdir, f"TandemMod_all.{assembly_label}{suffix}.score.bed")
    write_score_bed(out, rows)
    return out


def build_cheui(cheui_in, outdir, assembly_label, exon_dict, chr_dict, strand_dict,
                output_prefix="CHEUI_all", stdchr_only=True):
    rows = []

    with open_maybe_gz(cheui_in) as f:
        header = next(f).rstrip("\n").split("\t")
        idx = {x: i for i, x in enumerate(header)}

        contig_i = idx["contig"]
        pos_i = idx["position"]
        prob_i = idx["probability"]

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            contig = fields[contig_i]
            tx_id = contig.split("|")[0]

            # CHEUI position is treated as the 0-based start of the 9-mer.
            # +5 gives the 1-based transcript coordinate of the center A.
            kmer_start = int(float(fields[pos_i]))
            tx_pos = kmer_start + 5
            score = float(fields[prob_i])

            res = txpos_to_genome(tx_id, tx_pos, exon_dict, chr_dict, strand_dict)
            if res is None:
                continue

            chrom, genome_pos, strand = res

            if not is_std_chr(chrom, stdchr_only):
                continue

            rows.append((chrom, genome_pos - 1, genome_pos, score))

    suffix = ".stdchr" if stdchr_only else ""
    out = os.path.join(outdir, f"{output_prefix}.{assembly_label}.centerPlus5{suffix}.score.bed")
    write_score_bed(out, rows)
    return out


def build_m6anet(m6anet_in, outdir, assembly_label, exon_dict, chr_dict, strand_dict,
                 stdchr_only=True):
    rows = []

    with open_maybe_gz(m6anet_in, newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            tx_full = row["transcript_id"]
            tx_id = tx_full.split("|")[0]

            # m6Anet transcript_position is treated as 0-based.
            # txpos_to_genome expects 1-based transcript coordinate.
            tx_pos0 = int(float(row["transcript_position"]))
            tx_pos = tx_pos0 + 1
            score = float(row["probability_modified"])

            res = txpos_to_genome(tx_id, tx_pos, exon_dict, chr_dict, strand_dict)
            if res is None:
                continue

            chrom, genome_pos, strand = res

            if not is_std_chr(chrom, stdchr_only):
                continue

            rows.append((chrom, genome_pos - 1, genome_pos, score))

    suffix = ".stdchr" if stdchr_only else ""
    out = os.path.join(outdir, f"m6Anet_all.{assembly_label}.txPlus1{suffix}.score.bed")
    write_score_bed(out, rows)
    return out


def build_swarm(swarm_in, outdir, assembly_label, output_prefix="SWARM_all",
                shift=5, stdchr_only=True):
    rows = []

    with open_maybe_gz(swarm_in) as f:
        header = next(f).rstrip("\n").split("\t")
        idx = {x.lower(): i for i, x in enumerate(header)}

        score_candidates = [
            "probability", "score", "prob", "pred", "pred_score",
            "site_probability", "site_prob", "modscore", "mod_score"
        ]

        score_col = None
        for c in score_candidates:
            if c in idx:
                score_col = c
                break

        chr_candidates = ["chr", "chrom", "chromosome", "contig"]
        pos_candidates = ["genome_pos", "genomic_position", "position", "site", "pos"]

        chr_col = None
        pos_col = None

        for c in chr_candidates:
            if c in idx:
                chr_col = c
                break

        for c in pos_candidates:
            if c in idx:
                pos_col = c
                break

        if score_col is None or chr_col is None or pos_col is None:
            raise ValueError(
                f"Cannot parse SWARM columns. Header={header}, "
                f"chr_col={chr_col}, pos_col={pos_col}, score_col={score_col}"
            )

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            chrom = fields[idx[chr_col]]
            if not is_std_chr(chrom, stdchr_only):
                continue

            pos = int(float(fields[idx[pos_col]]))
            score = float(fields[idx[score_col]])

            # Final mouse benchmark coordinate convention:
            # SWARM genomic position + 5.
            genome_pos = pos + shift
            rows.append((chrom, genome_pos - 1, genome_pos, score))

    suffix = ".stdchr" if stdchr_only else ""
    shift_tag = f"shiftPlus{shift}" if shift >= 0 else f"shiftMinus{abs(shift)}"
    out = os.path.join(outdir, f"{output_prefix}.{assembly_label}.{shift_tag}{suffix}.score.bed")
    write_score_bed(out, rows)
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Build mouse m6A genome-level score.bed files with final coordinate conventions."
    )

    parser.add_argument("--gtf", required=True, help="Mouse GTF annotation file.")
    parser.add_argument("--tandemmod", required=True, help="TandemMod site-level TSV file.")
    parser.add_argument("--cheui-all", required=True, help="CHEUI all-site prediction TSV file.")
    parser.add_argument("--m6anet", required=True, help="m6Anet data.site_proba.csv file.")
    parser.add_argument("--swarm-all", required=True, help="SWARM all-site prediction TSV file.")

    parser.add_argument("--cheui-drach", default=None, help="Optional CHEUI DRACH-only TSV file.")
    parser.add_argument("--swarm-drach", default=None, help="Optional SWARM DRACH-only TSV file.")

    parser.add_argument("--outdir", required=True, help="Output scorebed directory.")
    parser.add_argument("--assembly-label", default="mm39", help="Assembly label used in output file names.")
    parser.add_argument("--swarm-shift", type=int, default=5, help="SWARM coordinate shift. Default: +5.")

    parser.add_argument(
        "--keep-all-chromosomes",
        action="store_true",
        help="Do not filter to mouse standard chromosomes.",
    )

    args = parser.parse_args()

    stdchr_only = not args.keep_all_chromosomes
    ensure_outdir(args.outdir)

    exon_dict, chr_dict, strand_dict = load_gtf_exons(args.gtf)
    print("Loaded transcripts from GTF:", len(exon_dict))

    outputs = []

    outputs.append(
        build_tandemmod(
            tandemmod_in=args.tandemmod,
            outdir=args.outdir,
            assembly_label=args.assembly_label,
            stdchr_only=stdchr_only,
        )
    )

    outputs.append(
        build_cheui(
            cheui_in=args.cheui_all,
            outdir=args.outdir,
            assembly_label=args.assembly_label,
            exon_dict=exon_dict,
            chr_dict=chr_dict,
            strand_dict=strand_dict,
            output_prefix="CHEUI_all",
            stdchr_only=stdchr_only,
        )
    )

    outputs.append(
        build_m6anet(
            m6anet_in=args.m6anet,
            outdir=args.outdir,
            assembly_label=args.assembly_label,
            exon_dict=exon_dict,
            chr_dict=chr_dict,
            strand_dict=strand_dict,
            stdchr_only=stdchr_only,
        )
    )

    outputs.append(
        build_swarm(
            swarm_in=args.swarm_all,
            outdir=args.outdir,
            assembly_label=args.assembly_label,
            output_prefix="SWARM_all",
            shift=args.swarm_shift,
            stdchr_only=stdchr_only,
        )
    )

    if args.cheui_drach is not None:
        outputs.append(
            build_cheui(
                cheui_in=args.cheui_drach,
                outdir=args.outdir,
                assembly_label=args.assembly_label,
                exon_dict=exon_dict,
                chr_dict=chr_dict,
                strand_dict=strand_dict,
                output_prefix="CHEUI_DRACH_all",
                stdchr_only=stdchr_only,
            )
        )

    if args.swarm_drach is not None:
        outputs.append(
            build_swarm(
                swarm_in=args.swarm_drach,
                outdir=args.outdir,
                assembly_label=args.assembly_label,
                output_prefix="SWARM_DRACH_all",
                shift=args.swarm_shift,
                stdchr_only=stdchr_only,
            )
        )

    list_out = os.path.join(args.outdir, "mouse_m6A_scorebeds.list")
    with open(list_out, "w") as f:
        f.write("name\tpath\n")
        for path in outputs:
            name = os.path.basename(path).replace(".score.bed", "")
            f.write(f"{name}\t{path}\n")

    print("saved list:", list_out)


if __name__ == "__main__":
    main()