#!/usr/bin/env python3
import argparse
import os
import re
import gzip
import pandas as pd
from pathlib import Path
from collections import Counter


DRACH = re.compile(r"^[AGT][AG]AC[ACT]$")


def open_maybe_gz(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def get_standard_chromosomes(include_chrM=True):
    chrs = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    if include_chrM:
        chrs.add("chrM")
    return chrs


def load_genome(genome_fa, std_chrs):
    genome_fa = str(genome_fa)

    if not os.path.exists(genome_fa):
        raise FileNotFoundError(f"Genome FASTA not found: {genome_fa}")

    genome = {}
    chrom = None
    seqs = []

    with open_maybe_gz(genome_fa) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if chrom in std_chrs:
                    genome[chrom] = "".join(seqs).upper()

                chrom = line[1:].split()[0]
                seqs = []
            else:
                if chrom in std_chrs:
                    seqs.append(line.strip())

        if chrom in std_chrs:
            genome[chrom] = "".join(seqs).upper()

    print("Loaded genome:", genome_fa)
    print("Loaded chromosomes:", len(genome))
    return genome


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def is_drach_site(genome, chrom, start, end, strand):
    """
    Check whether a single-base BED site is centered in a DRACH motif.

    BED coordinate:
    chrom, start, end
    where end should equal start + 1.
    """
    if chrom not in genome:
        return False

    if end != start + 1:
        return False

    if start < 2:
        return False

    seq = genome[chrom]

    if start + 3 > len(seq):
        return False

    mer = seq[start - 2:start + 3].upper()

    if str(strand).strip() == "-":
        mer = revcomp(mer)

    return bool(DRACH.match(mer))


def is_drach_motif(motif):
    motif = str(motif).upper().replace("U", "T")
    return bool(DRACH.match(motif))


def write_bed(path, sites):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as out:
        for chrom, start, end in sorted(sites, key=lambda x: (x[0], x[1], x[2])):
            out.write(f"{chrom}\t{start}\t{end}\n")

    print("saved:", path, "n =", len(sites))


def read_bed(path):
    sites = set()

    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue

            fs = line.rstrip("\n").split("\t")
            if len(fs) < 3:
                continue

            sites.add((fs[0], int(fs[1]), int(fs[2])))

    return sites


def clean_columns(df):
    df.columns = df.columns.astype(str).str.strip().str.replace("\ufeff", "", regex=False)
    return df


def load_glori_rep(
    path,
    genome,
    std_chrs,
    min_agcov=10,
    max_padj=0.05,
    min_ratio=0.2,
    chr_col="Chr",
    site_col="Sites",
    strand_col="Strand",
    ratio_col="Ratio",
    padj_col="P_adjust",
    agcov_col=None,
):
    """
    Load one GLORI replicate and return DRACH-positive sites.

    This preserves the original human m6A logic:
    Ratio >= min_ratio, P_adjust < max_padj, AGcov >= min_agcov,
    then genome-based DRACH motif filtering.
    """
    df = pd.read_csv(path, sep=None, engine="python")
    df = clean_columns(df)

    required = [chr_col, site_col, strand_col, ratio_col, padj_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}. Available columns: {list(df.columns)}")

    if agcov_col is None:
        if "AGcov" in df.columns:
            agcov_col = "AGcov"
        elif "Acov" in df.columns:
            agcov_col = "Acov"
        else:
            raise ValueError(f"No AG coverage column found in {path}. Expected AGcov or Acov.")

    df[chr_col] = df[chr_col].astype(str)
    df[site_col] = pd.to_numeric(df[site_col], errors="coerce")
    df[ratio_col] = pd.to_numeric(df[ratio_col], errors="coerce")
    df[padj_col] = pd.to_numeric(df[padj_col], errors="coerce")
    df[agcov_col] = pd.to_numeric(df[agcov_col], errors="coerce")

    df = df[
        df[chr_col].isin(std_chrs)
        & df[site_col].notna()
        & df[ratio_col].notna()
        & df[padj_col].notna()
        & df[agcov_col].notna()
        & (df[padj_col] < max_padj)
        & (df[ratio_col] >= min_ratio)
        & (df[agcov_col] >= min_agcov)
    ]

    sites = set()

    for _, r in df.iterrows():
        chrom = r[chr_col]
        pos1 = int(r[site_col])
        start = pos1 - 1
        end = pos1
        strand = r[strand_col]

        if is_drach_site(genome, chrom, start, end, strand):
            sites.add((chrom, start, end))

    print(os.path.basename(str(path)), "positive DRACH sites:", len(sites))
    return sites


def load_m6aboost_bed(path, genome, std_chrs):
    """
    Load m6Aboost BED/BED.gz and keep genome-confirmed DRACH single-base sites.
    """
    sites = set()

    with open_maybe_gz(path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("track") or line.startswith("#") or line.startswith("browser"):
                continue

            fs = line.rstrip("\n").split()
            if len(fs) < 3:
                continue

            chrom = fs[0]
            if chrom not in std_chrs:
                continue

            start = int(fs[1])
            end = int(fs[2])
            strand = fs[5] if len(fs) >= 6 else "+"

            if end != start + 1:
                center = (start + end) // 2
                start, end = center, center + 1

            if is_drach_site(genome, chrom, start, end, strand):
                sites.add((chrom, start, end))

    print("m6Aboost positive DRACH sites:", len(sites))
    return sites


def load_m6ace_csv(
    path,
    genome,
    std_chrs,
    chr_col="Chr",
    start_col="Start",
    end_col="End",
    strand_col="Strand",
    motif_col="Motif",
):
    """
    Load m6ACE CSV/TSV and keep genome-confirmed DRACH single-base sites.
    If a motif column exists, it is used as an additional DRACH check.
    """
    df = pd.read_csv(path, sep=None, engine="python")
    df = clean_columns(df)

    required = [chr_col, start_col, end_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}. Available columns: {list(df.columns)}")

    df[chr_col] = df[chr_col].astype(str)
    df[start_col] = pd.to_numeric(df[start_col], errors="coerce")
    df[end_col] = pd.to_numeric(df[end_col], errors="coerce")

    sites = set()

    for _, r in df.iterrows():
        chrom = r[chr_col]
        if chrom not in std_chrs:
            continue

        if pd.isna(r[start_col]) or pd.isna(r[end_col]):
            continue

        start = int(r[start_col])
        end = int(r[end_col])

        strand = r[strand_col] if strand_col in df.columns else "+"

        if end != start + 1:
            center = (start + end) // 2
            start, end = center, center + 1

        if motif_col in df.columns and not is_drach_motif(r[motif_col]):
            continue

        if is_drach_site(genome, chrom, start, end, strand):
            sites.add((chrom, start, end))

    print("m6ACE positive DRACH sites:", len(sites))
    return sites


def build_reference_sets(source_out, ref_out, level):
    """
    Build final reference sets for one level.

    strict level:
    GLORI.strict + m6Aboost.strict + m6ACE.strict
    → strict_consensus_2of3

    loose level:
    GLORI.loose + m6Aboost.loose + m6ACE.loose
    → loose_union_3way
    """
    glori = read_bed(source_out / f"GLORI.{level}.bed")
    boost = read_bed(source_out / f"m6Aboost.{level}.bed")
    ace = read_bed(source_out / f"m6ACE.{level}.bed")

    counter = Counter()
    for s in [glori, boost, ace]:
        counter.update(s)

    consensus = {site for site, n in counter.items() if n >= 2}
    union = glori | boost | ace

    print("\n===", level, "===")
    print("GLORI:", len(glori))
    print("m6Aboost:", len(boost))
    print("m6ACE:", len(ace))
    print("consensus_2of3:", len(consensus))
    print("union_3way:", len(union))

    if level == "strict":
        write_bed(ref_out / "strict_consensus_2of3.bed", consensus)
        return {
            "strict_consensus_2of3": len(consensus),
            "strict_GLORI": len(glori),
            "strict_m6Aboost": len(boost),
            "strict_m6ACE": len(ace),
        }

    if level == "loose":
        write_bed(ref_out / "loose_union_3way.bed", union)
        return {
            "loose_union_3way": len(union),
            "loose_GLORI": len(glori),
            "loose_m6Aboost": len(boost),
            "loose_m6ACE": len(ace),
        }

    raise ValueError(f"Unexpected level: {level}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build human m6A DRACH reference sets from GLORI, m6Aboost, and m6ACE sources. "
            "This is a parameterized version of build_human_truth_DRACH_simple.py."
        )
    )

    parser.add_argument("--genome-fa", required=True, help="GRCh38 genome FASTA, plain or .gz.")
    parser.add_argument("--glori-rep1", required=True, help="GLORI HEK293T replicate 1 csv/tsv, plain or .gz.")
    parser.add_argument("--glori-rep2", required=True, help="GLORI HEK293T replicate 2 csv/tsv, plain or .gz.")
    parser.add_argument("--m6aboost-bed", required=True, help="m6Aboost BED/BED.gz file.")
    parser.add_argument("--m6ace-csv", required=True, help="m6ACE CSV/TSV file.")
    parser.add_argument("--outdir", required=True, help="Output directory.")

    parser.add_argument("--min-agcov", type=int, default=10)
    parser.add_argument("--max-padj", type=float, default=0.05)
    parser.add_argument("--glori-min-ratio", type=float, default=0.2)

    parser.add_argument("--glori-chr-col", default="Chr")
    parser.add_argument("--glori-site-col", default="Sites")
    parser.add_argument("--glori-strand-col", default="Strand")
    parser.add_argument("--glori-ratio-col", default="Ratio")
    parser.add_argument("--glori-padj-col", default="P_adjust")
    parser.add_argument("--glori-agcov-col", default=None)

    parser.add_argument("--m6ace-chr-col", default="Chr")
    parser.add_argument("--m6ace-start-col", default="Start")
    parser.add_argument("--m6ace-end-col", default="End")
    parser.add_argument("--m6ace-strand-col", default="Strand")
    parser.add_argument("--m6ace-motif-col", default="Motif")

    parser.add_argument(
        "--no-chrM",
        action="store_true",
        help="Exclude chrM from standard chromosomes.",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    source_out = outdir / "truth_sources_DRACH"
    ref_out = outdir / "reference_sets_DRACH"

    source_out.mkdir(parents=True, exist_ok=True)
    ref_out.mkdir(parents=True, exist_ok=True)

    std_chrs = get_standard_chromosomes(include_chrM=not args.no_chrM)
    genome = load_genome(args.genome_fa, std_chrs)

    print("\n=== GLORI reps ===")
    glori1 = load_glori_rep(
        path=args.glori_rep1,
        genome=genome,
        std_chrs=std_chrs,
        min_agcov=args.min_agcov,
        max_padj=args.max_padj,
        min_ratio=args.glori_min_ratio,
        chr_col=args.glori_chr_col,
        site_col=args.glori_site_col,
        strand_col=args.glori_strand_col,
        ratio_col=args.glori_ratio_col,
        padj_col=args.glori_padj_col,
        agcov_col=args.glori_agcov_col,
    )

    glori2 = load_glori_rep(
        path=args.glori_rep2,
        genome=genome,
        std_chrs=std_chrs,
        min_agcov=args.min_agcov,
        max_padj=args.max_padj,
        min_ratio=args.glori_min_ratio,
        chr_col=args.glori_chr_col,
        site_col=args.glori_site_col,
        strand_col=args.glori_strand_col,
        ratio_col=args.glori_ratio_col,
        padj_col=args.glori_padj_col,
        agcov_col=args.glori_agcov_col,
    )

    glori_strict = glori1 & glori2
    glori_loose = glori1 | glori2

    print("GLORI.strict rep intersection:", len(glori_strict))
    print("GLORI.loose rep union:", len(glori_loose))

    print("\n=== m6Aboost ===")
    boost = load_m6aboost_bed(args.m6aboost_bed, genome, std_chrs)

    print("\n=== m6ACE ===")
    ace = load_m6ace_csv(
        path=args.m6ace_csv,
        genome=genome,
        std_chrs=std_chrs,
        chr_col=args.m6ace_chr_col,
        start_col=args.m6ace_start_col,
        end_col=args.m6ace_end_col,
        strand_col=args.m6ace_strand_col,
        motif_col=args.m6ace_motif_col,
    )

    write_bed(source_out / "GLORI.strict.bed", glori_strict)
    write_bed(source_out / "GLORI.loose.bed", glori_loose)

    write_bed(source_out / "m6Aboost.strict.bed", boost)
    write_bed(source_out / "m6Aboost.loose.bed", boost)

    write_bed(source_out / "m6ACE.strict.bed", ace)
    write_bed(source_out / "m6ACE.loose.bed", ace)

    strict_counts = build_reference_sets(source_out, ref_out, "strict")
    loose_counts = build_reference_sets(source_out, ref_out, "loose")

    counts_out = ref_out / "reference_set_counts.tsv"
    with open(counts_out, "w") as f:
        f.write("set_name\tsite_count\n")
        for k, v in {**strict_counts, **loose_counts}.items():
            f.write(f"{k}\t{v}\n")

    print("\nSaved summary:", counts_out)

    print("\nFinal source BEDs:")
    for f in sorted(os.listdir(source_out)):
        if f.endswith(".bed"):
            p = source_out / f
            print(f, sum(1 for _ in open(p)))

    print("\nFinal reference BEDs:")
    for f in sorted(os.listdir(ref_out)):
        if f.endswith(".bed"):
            p = ref_out / f
            print(f, sum(1 for _ in open(p)))


if __name__ == "__main__":
    main()