#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare low-pass VCF genotypes against WGS truth VCF and compute concordance metrics."
    )

    parser.add_argument(
        "--truth-vcf",
        required=True,
        help="WGS truth VCF/BCF, bgzipped and indexed.",
    )
    parser.add_argument(
        "--query-vcf",
        required=True,
        help="Pipeline output VCF/BCF to evaluate, bgzipped and indexed.",
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample name in query VCF. Also used for truth unless --truth-sample is provided.",
    )
    parser.add_argument(
        "--truth-sample",
        default=None,
        help="Sample name in truth VCF. Default: same as --sample.",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output prefix. Creates .metrics.tsv, .confusion.tsv, .sites.tsv.gz.",
    )
    parser.add_argument(
        "--coverage",
        type=float,
        required=True,
        help="Target coverage label as numeric value.",
    )
    parser.add_argument(
        "--replicate",
        type=int,
        required=True,
        help="Replicate number.",
    )
    parser.add_argument(
        "--method",
        default="pipeline",
        help="Method label, e.g. bcftools, glimpse2, varscan.",
    )
    parser.add_argument(
        "--regions",
        default=None,
        help="Optional BED file or region string is not implemented here. Placeholder.",
    )
    parser.add_argument(
        "--snps-only",
        action="store_true",
        help="Compare SNP sites only.",
    )
    parser.add_argument(
        "--biallelic-only",
        action="store_true",
        help="Compare biallelic sites only.",
    )
    parser.add_argument(
        "--require-pass",
        action="store_true",
        help="Use only PASS records in both VCFs.",
    )
    parser.add_argument(
        "--min-truth-gq",
        type=float,
        default=None,
        help="Optional minimum GQ in truth sample.",
    )
    parser.add_argument(
        "--min-query-gq",
        type=float,
        default=None,
        help="Optional minimum GQ in query sample.",
    )
    parser.add_argument(
        "--min-truth-dp",
        type=float,
        default=None,
        help="Optional minimum DP in truth sample.",
    )
    parser.add_argument(
        "--min-query-dp",
        type=float,
        default=None,
        help="Optional minimum DP in query sample.",
    )

    return parser.parse_args()


def normalize_gt(gt) -> Optional[Tuple[int, ...]]:
    """
    Convert pysam GT tuple to sorted allele tuple.
    Examples:
      (0, 0) -> (0, 0)
      (1, 0) -> (0, 1)
      (1, 1) -> (1, 1)
      (None, 1) -> None
    """
    if gt is None:
        return None
    if any(a is None for a in gt):
        return None
    if len(gt) == 0:
        return None
    return tuple(sorted(int(a) for a in gt))


def gt_class(gt: Optional[Tuple[int, ...]]) -> str:
    if gt is None:
        return "missing"

    if len(gt) != 2:
        return "non_diploid"

    a, b = gt

    if a == 0 and b == 0:
        return "hom_ref"

    if a != b:
        return "het"

    if a == b and a > 0:
        return "hom_alt"

    return "other"


def has_alt(gt: Optional[Tuple[int, ...]]) -> Optional[bool]:
    if gt is None:
        return None
    return any(a > 0 for a in gt)


def is_het(gt: Optional[Tuple[int, ...]]) -> Optional[bool]:
    if gt is None or len(gt) != 2:
        return None
    return gt[0] != gt[1]


def is_snp_record(rec) -> bool:
    if len(rec.ref) != 1:
        return False
    if rec.alts is None:
        return False
    return all(len(alt) == 1 for alt in rec.alts)


def is_biallelic_record(rec) -> bool:
    return rec.alts is not None and len(rec.alts) == 1


def pass_filter(rec) -> bool:
    # In pysam, rec.filter.keys() can be empty for PASS in some VCFs.
    keys = list(rec.filter.keys())
    return len(keys) == 0 or keys == ["PASS"]


def get_sample_value(sample_data, key: str):
    try:
        return sample_data.get(key)
    except Exception:
        return None


def sample_passes_thresholds(sample_data, min_gq=None, min_dp=None) -> bool:
    if min_gq is not None:
        gq = get_sample_value(sample_data, "GQ")
        if gq is None:
            return False
        if float(gq) < min_gq:
            return False

    if min_dp is not None:
        dp = get_sample_value(sample_data, "DP")
        if dp is None:
            return False
        if float(dp) < min_dp:
            return False

    return True


def make_site_key(rec) -> Tuple[str, int, str, Tuple[str, ...]]:
    alts = tuple(rec.alts) if rec.alts is not None else tuple()
    return rec.chrom, rec.pos, rec.ref, alts


def build_truth_index(
    truth_vcf_path: str,
    truth_sample: str,
    snps_only: bool,
    biallelic_only: bool,
    require_pass: bool,
    min_truth_gq,
    min_truth_dp,
) -> Dict[Tuple[str, int, str, Tuple[str, ...]], Tuple[Tuple[int, ...], str]]:
    truth_vcf = pysam.VariantFile(truth_vcf_path)

    if truth_sample not in truth_vcf.header.samples:
        raise ValueError(f"Truth sample '{truth_sample}' not found in {truth_vcf_path}")

    index = {}

    for rec in truth_vcf.fetch():
        if require_pass and not pass_filter(rec):
            continue
        if snps_only and not is_snp_record(rec):
            continue
        if biallelic_only and not is_biallelic_record(rec):
            continue

        sample_data = rec.samples[truth_sample]

        if not sample_passes_thresholds(
            sample_data,
            min_gq=min_truth_gq,
            min_dp=min_truth_dp,
        ):
            continue

        gt = normalize_gt(sample_data.get("GT"))
        if gt is None:
            continue

        key = make_site_key(rec)
        index[key] = (gt, gt_class(gt))

    truth_vcf.close()
    return index


def rate(num: int, den: int) -> Optional[float]:
    if den == 0:
        return None
    return num / den


def main():
    args = parse_args()

    truth_sample = args.truth_sample or args.sample
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    truth_index = build_truth_index(
        truth_vcf_path=args.truth_vcf,
        truth_sample=truth_sample,
        snps_only=args.snps_only,
        biallelic_only=args.biallelic_only,
        require_pass=args.require_pass,
        min_truth_gq=args.min_truth_gq,
        min_truth_dp=args.min_truth_dp,
    )

    query_vcf = pysam.VariantFile(args.query_vcf)

    if args.sample not in query_vcf.header.samples:
        raise ValueError(f"Query sample '{args.sample}' not found in {args.query_vcf}")

    counters = Counter()
    confusion = Counter()

    sites_path = Path(str(output_prefix) + ".sites.tsv.gz")
    with gzip.open(sites_path, "wt") as sites_out:
        sites_writer = csv.DictWriter(
            sites_out,
            fieldnames=[
                "sample",
                "coverage",
                "replicate",
                "method",
                "chrom",
                "pos",
                "ref",
                "alt",
                "truth_gt",
                "query_gt",
                "truth_class",
                "query_class",
                "exact_match",
                "alt_presence_match",
                "het_status_match",
            ],
            delimiter="\t",
        )
        sites_writer.writeheader()

        for rec in query_vcf.fetch():
            counters["query_records_total"] += 1

            if args.require_pass and not pass_filter(rec):
                counters["query_records_filter_fail"] += 1
                continue
            if args.snps_only and not is_snp_record(rec):
                counters["query_records_non_snp"] += 1
                continue
            if args.biallelic_only and not is_biallelic_record(rec):
                counters["query_records_non_biallelic"] += 1
                continue

            key = make_site_key(rec)

            if key not in truth_index:
                counters["query_records_not_in_truth"] += 1
                continue

            sample_data = rec.samples[args.sample]

            if not sample_passes_thresholds(
                sample_data,
                min_gq=args.min_query_gq,
                min_dp=args.min_query_dp,
            ):
                counters["query_records_sample_threshold_fail"] += 1
                continue

            query_gt = normalize_gt(sample_data.get("GT"))
            if query_gt is None:
                counters["query_gt_missing"] += 1
                continue

            truth_gt, truth_cls = truth_index[key]
            query_cls = gt_class(query_gt)

            counters["n_compared"] += 1

            exact = truth_gt == query_gt
            alt_match = has_alt(truth_gt) == has_alt(query_gt)
            het_match = is_het(truth_gt) == is_het(query_gt)

            if exact:
                counters["n_exact_match"] += 1
            if alt_match:
                counters["n_alt_presence_match"] += 1
            if het_match:
                counters["n_het_status_match"] += 1

            if truth_cls == "hom_ref":
                counters["truth_hom_ref"] += 1
                if query_cls == "hom_ref":
                    counters["truth_hom_ref_called_hom_ref"] += 1

            if truth_cls == "het":
                counters["truth_het"] += 1
                if query_cls == "het":
                    counters["truth_het_called_het"] += 1

            if truth_cls == "hom_alt":
                counters["truth_hom_alt"] += 1
                if query_cls == "hom_alt":
                    counters["truth_hom_alt_called_hom_alt"] += 1

            if has_alt(truth_gt):
                counters["truth_nonref"] += 1
                if exact:
                    counters["truth_nonref_exact"] += 1
                if has_alt(query_gt):
                    counters["truth_nonref_called_nonref"] += 1

            confusion[(truth_cls, query_cls)] += 1

            sites_writer.writerow(
                {
                    "sample": args.sample,
                    "coverage": args.coverage,
                    "replicate": args.replicate,
                    "method": args.method,
                    "chrom": rec.chrom,
                    "pos": rec.pos,
                    "ref": rec.ref,
                    "alt": ",".join(rec.alts or []),
                    "truth_gt": "/".join(map(str, truth_gt)),
                    "query_gt": "/".join(map(str, query_gt)),
                    "truth_class": truth_cls,
                    "query_class": query_cls,
                    "exact_match": int(exact),
                    "alt_presence_match": int(alt_match),
                    "het_status_match": int(het_match),
                }
            )

    query_vcf.close()

    n = counters["n_compared"]

    metrics = {
        "sample": args.sample,
        "truth_sample": truth_sample,
        "coverage": args.coverage,
        "replicate": args.replicate,
        "method": args.method,
        "truth_vcf": args.truth_vcf,
        "query_vcf": args.query_vcf,
        "n_truth_index_sites": len(truth_index),
        "n_query_records_total": counters["query_records_total"],
        "n_compared": n,
        "n_exact_match": counters["n_exact_match"],
        "exact_match_rate": rate(counters["n_exact_match"], n),
        "n_alt_presence_match": counters["n_alt_presence_match"],
        "alt_presence_match_rate": rate(counters["n_alt_presence_match"], n),
        "n_het_status_match": counters["n_het_status_match"],
        "het_status_match_rate": rate(counters["n_het_status_match"], n),
        "truth_hom_ref": counters["truth_hom_ref"],
        "hom_ref_match_rate": rate(
            counters["truth_hom_ref_called_hom_ref"],
            counters["truth_hom_ref"],
        ),
        "truth_het": counters["truth_het"],
        "het_recall": rate(
            counters["truth_het_called_het"],
            counters["truth_het"],
        ),
        "truth_hom_alt": counters["truth_hom_alt"],
        "hom_alt_recall": rate(
            counters["truth_hom_alt_called_hom_alt"],
            counters["truth_hom_alt"],
        ),
        "truth_nonref": counters["truth_nonref"],
        "nonref_exact_match_rate": rate(
            counters["truth_nonref_exact"],
            counters["truth_nonref"],
        ),
        "alt_recall": rate(
            counters["truth_nonref_called_nonref"],
            counters["truth_nonref"],
        ),
        "query_records_not_in_truth": counters["query_records_not_in_truth"],
        "query_gt_missing": counters["query_gt_missing"],
    }

    metrics_path = Path(str(output_prefix) + ".metrics.tsv")
    with metrics_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(metrics.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(metrics)

    confusion_path = Path(str(output_prefix) + ".confusion.tsv")
    with confusion_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["sample", "coverage", "replicate", "method", "truth_class", "query_class", "count"],
            delimiter="\t",
        )
        writer.writeheader()

        for (truth_cls, query_cls), count in sorted(confusion.items()):
            writer.writerow(
                {
                    "sample": args.sample,
                    "coverage": args.coverage,
                    "replicate": args.replicate,
                    "method": args.method,
                    "truth_class": truth_cls,
                    "query_class": query_cls,
                    "count": count,
                }
            )

    json_path = Path(str(output_prefix) + ".metrics.json")
    with json_path.open("w") as f:
        json.dump(metrics, f, indent=2)

    print(f"Metrics:   {metrics_path}")
    print(f"Confusion: {confusion_path}")
    print(f"Sites:     {sites_path}")


if __name__ == "__main__":
    main()
