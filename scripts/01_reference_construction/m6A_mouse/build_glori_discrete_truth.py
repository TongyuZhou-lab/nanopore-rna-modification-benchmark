#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path


def open_text(path):
    """
    Open plain text or gzip-compressed text file.
    """
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def normalize_delimiter(delim):
    """
    Allow users to pass '\\t' from command line.
    """
    if delim == "\\t":
        return "\t"
    return delim


def check_required_columns(reader, required_cols, path):
    """
    Make sure all required columns exist in the input file.
    """
    fieldnames = reader.fieldnames or []
    missing = [c for c in required_cols if c not in fieldnames]

    if missing:
        raise ValueError(
            f"Missing required columns in {path}: {missing}\n"
            f"Available columns: {fieldnames}"
        )


def read_glori(
    path,
    rep_name,
    delimiter,
    chr_col,
    site_col,
    strand_col,
    ratio_col,
    agcov_col,
    padj_col,
):
    """
    Read one GLORI replicate.

    The site key is:
    chrom, 1-based site position, strand

    This follows the same logic as the original mouse script.
    """
    data = {}

    required_cols = [
        chr_col,
        site_col,
        strand_col,
        ratio_col,
        agcov_col,
        padj_col,
    ]

    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        check_required_columns(reader, required_cols, path)

        for row in reader:
            chrom = row[chr_col]
            site = int(float(row[site_col]))
            strand = row[strand_col]

            key = (chrom, site, strand)

            data[key] = {
                f"Ratio_{rep_name}": float(row[ratio_col]),
                f"AGcov_{rep_name}": int(float(row[agcov_col])),
                f"Padj_{rep_name}": float(row[padj_col]),
            }

    return data


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build GLORI-derived discrete m6A truth sets from two replicates. "
            "The script outputs a mean-ratio table, a loose BED set, and a strict BED set."
        )
    )

    parser.add_argument(
        "--rep1",
        required=True,
        help="GLORI replicate 1 file, plain text or .gz.",
    )
    parser.add_argument(
        "--rep2",
        required=True,
        help="GLORI replicate 2 file, plain text or .gz.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "--prefix",
        default="GLORI",
        help="Output prefix, e.g. GLORI_mESC or GLORI_HEK293T.",
    )

    parser.add_argument(
        "--delimiter",
        default="\\t",
        help="Input delimiter. Default: tab. Use ',' for comma-separated files.",
    )

    parser.add_argument("--chr-col", default="Chr", help="Chromosome column name.")
    parser.add_argument("--site-col", default="Sites", help="1-based site position column name.")
    parser.add_argument("--strand-col", default="Strand", help="Strand column name.")
    parser.add_argument("--ratio-col", default="Ratio", help="Modification ratio column name.")
    parser.add_argument("--agcov-col", default="AGcov", help="AG coverage column name.")
    parser.add_argument("--padj-col", default="P_adjust", help="Adjusted P-value column name.")

    parser.add_argument(
        "--min-agcov",
        type=int,
        default=10,
        help="Minimum AG coverage required in both replicates. Default: 10.",
    )
    parser.add_argument(
        "--max-padj",
        type=float,
        default=0.05,
        help="Maximum adjusted P-value required in both replicates. Default: 0.05.",
    )
    parser.add_argument(
        "--loose-ratio",
        type=float,
        default=0.2,
        help="Mean ratio cutoff for loose GLORI positive BED. Default: 0.2.",
    )
    parser.add_argument(
        "--strict-ratio",
        type=float,
        default=0.95,
        help="Mean ratio cutoff for strict GLORI positive BED. Default: 0.95.",
    )

    args = parser.parse_args()

    delimiter = normalize_delimiter(args.delimiter)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    out_mean = outdir / f"{args.prefix}_mean_ratio.tsv"
    out_loose = outdir / f"{args.prefix}_positive_meanRatio{args.loose_ratio}.bed"
    out_strict = outdir / f"{args.prefix}_positive_meanRatio{args.strict_ratio}.bed"

    rep1 = read_glori(
        path=args.rep1,
        rep_name="rep1",
        delimiter=delimiter,
        chr_col=args.chr_col,
        site_col=args.site_col,
        strand_col=args.strand_col,
        ratio_col=args.ratio_col,
        agcov_col=args.agcov_col,
        padj_col=args.padj_col,
    )

    rep2 = read_glori(
        path=args.rep2,
        rep_name="rep2",
        delimiter=delimiter,
        chr_col=args.chr_col,
        site_col=args.site_col,
        strand_col=args.strand_col,
        ratio_col=args.ratio_col,
        agcov_col=args.agcov_col,
        padj_col=args.padj_col,
    )

    shared_keys = sorted(
        set(rep1.keys()) & set(rep2.keys()),
        key=lambda x: (x[0], x[1], x[2]),
    )

    base_filtered = 0
    loose_count = 0
    strict_count = 0

    with open(out_mean, "w", newline="") as fout_mean, \
         open(out_loose, "w", newline="") as fout_loose, \
         open(out_strict, "w", newline="") as fout_strict:

        mean_writer = csv.writer(fout_mean, delimiter="\t")
        loose_writer = csv.writer(fout_loose, delimiter="\t")
        strict_writer = csv.writer(fout_strict, delimiter="\t")

        mean_writer.writerow([
            "Chr", "Sites", "Strand",
            "Ratio_rep1", "Ratio_rep2", "mean_ratio",
            "AGcov_rep1", "AGcov_rep2",
            "Padj_rep1", "Padj_rep2",
        ])

        for key in shared_keys:
            chrom, site, strand = key

            r1 = rep1[key]
            r2 = rep2[key]

            ratio1 = r1["Ratio_rep1"]
            ratio2 = r2["Ratio_rep2"]
            ag1 = r1["AGcov_rep1"]
            ag2 = r2["AGcov_rep2"]
            p1 = r1["Padj_rep1"]
            p2 = r2["Padj_rep2"]

            mean_ratio = (ratio1 + ratio2) / 2.0

            mean_writer.writerow([
                chrom, site, strand,
                ratio1, ratio2, mean_ratio,
                ag1, ag2, p1, p2,
            ])

            if (
                ag1 >= args.min_agcov
                and ag2 >= args.min_agcov
                and p1 < args.max_padj
                and p2 < args.max_padj
            ):
                base_filtered += 1

                # GLORI site is treated as 1-based single-base position.
                # BED uses 0-based half-open coordinates.
                start = site - 1
                end = site

                if mean_ratio >= args.loose_ratio:
                    loose_writer.writerow([chrom, start, end, mean_ratio, mean_ratio, strand])
                    loose_count += 1

                if mean_ratio >= args.strict_ratio:
                    strict_writer.writerow([chrom, start, end, mean_ratio, mean_ratio, strand])
                    strict_count += 1

    print("shared_sites =", len(shared_keys))
    print("base_filtered_sites =", base_filtered)
    print(f"positive_sites_meanRatio{args.loose_ratio} =", loose_count)
    print(f"positive_sites_meanRatio{args.strict_ratio} =", strict_count)
    print("saved:", out_mean)
    print("saved:", out_loose)
    print("saved:", out_strict)


if __name__ == "__main__":
    main()