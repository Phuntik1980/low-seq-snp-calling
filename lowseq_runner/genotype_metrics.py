import csv
import gzip
import json
import math
from collections import Counter
from pathlib import Path
from typing import Optional

import pysam


def normalize_gt(gt) -> Optional[tuple[int, ...]]:
    if gt is None:
        return None
    if any(a is None for a in gt):
        return None
    if len(gt) == 0:
        return None
    return tuple(sorted(int(a) for a in gt))


def gt_class(gt: Optional[tuple[int, ...]]) -> str:
    """
    Define classification rules for hom_ref/het/hom_alt
    """
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


def has_alt(gt: Optional[tuple[int, ...]]) -> Optional[bool]:
    if gt is None:
        return None
    return any(a > 0 for a in gt)


def alt_dosage(gt: Optional[tuple[int, ...]]) -> Optional[int]:
    """
    ALT dosage for biallelic/diploid comparison:
      0/0 -> 0
      0/1 -> 1
      1/1 -> 2

    For multiallelic sites, any allele > 0 is counted as ALT allele.
    But for strict interpretation use --biallelic-only.
    """
    if gt is None:
        return None
    if len(gt) != 2:
        return None
    return sum(1 for a in gt if a > 0)


def is_het(gt: Optional[tuple[int, ...]]) -> Optional[bool]:
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
        if gq is None or float(gq) < min_gq:
            return False

    if min_dp is not None:
        dp = get_sample_value(sample_data, "DP")
        if dp is None or float(dp) < min_dp:
            return False

    return True


def make_site_key(rec) -> tuple[str, int, str, tuple[str, ...]]:
    alts = tuple(rec.alts) if rec.alts is not None else tuple()
    return rec.chrom, rec.pos, rec.ref, alts


def build_truth_index(
    truth_vcf_path: str | Path,
    truth_sample: str,
    snps_only: bool,
    biallelic_only: bool,
    require_pass: bool,
    min_truth_gq,
    min_truth_dp,
) -> dict[tuple[str, int, str, tuple[str, ...]], tuple[tuple[int, ...], str]]:
    """
    Explain why sites without GT in truth are excluded
    """
    truth_vcf = pysam.VariantFile(str(truth_vcf_path))

    if truth_sample not in truth_vcf.header.samples:
        raise ValueError(
            f"Truth sample '{truth_sample}' not in {truth_vcf_path}")

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

        index[make_site_key(rec)] = (gt, gt_class(gt))

    truth_vcf.close()
    return index


def rate(num: int, den: int) -> Optional[float]:
    if den == 0:
        return None
    return num / den


def f1_score(precision: Optional[float],
             recall: Optional[float]) -> Optional[float]:
    if precision is None or recall is None:
        return None
    if precision + recall == 0:
        return None
    return 2.0 * precision * recall / (precision + recall)


def pearson_r(xs: list[float], ys: list[float]) -> Optional[float]:
    if len(xs) < 2 or len(ys) < 2:
        return None

    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)

    dx = [x - mean_x for x in xs]
    dy = [y - mean_y for y in ys]

    ss_x = sum(x * x for x in dx)
    ss_y = sum(y * y for y in dy)

    if ss_x == 0 or ss_y == 0:
        return None

    cov = sum(x * y for x, y in zip(dx, dy))
    return cov / math.sqrt(ss_x * ss_y)


def pearson_r2(xs: list[float], ys: list[float]) -> Optional[float]:
    r = pearson_r(xs, ys)
    if r is None:
        return None
    return r * r


def compute_genotype_metrics(
    truth_vcf: str | Path,
    query_vcf: str | Path,
    sample: str,
    output_prefix: str | Path,
    coverage: float,
    replicate: int,
    method: str = "pipeline",
    truth_sample: Optional[str] = None,
    snps_only: bool = False,
    biallelic_only: bool = False,
    require_pass: bool = False,
    min_truth_gq: Optional[float] = None,
    min_query_gq: Optional[float] = None,
    min_truth_dp: Optional[float] = None,
    min_query_dp: Optional[float] = None,
    write_sites: bool = False,
) -> dict:
    truth_sample = truth_sample or sample
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    truth_index = build_truth_index(
        truth_vcf_path=truth_vcf,
        truth_sample=truth_sample,
        snps_only=snps_only,
        biallelic_only=biallelic_only,
        require_pass=require_pass,
        min_truth_gq=min_truth_gq,
        min_truth_dp=min_truth_dp,
    )

    query_vcf_obj = pysam.VariantFile(str(query_vcf))

    if sample not in query_vcf_obj.header.samples:
        raise ValueError(f"Query sample '{sample}' not in {query_vcf}")

    counters = Counter()
    confusion = Counter()

    truth_dosages: list[float] = []
    query_dosages: list[float] = []

    sites_writer = None
    sites_out = None
    sites_path = Path(str(output_prefix) + ".sites.tsv.gz")

    if write_sites:
        sites_out = gzip.open(sites_path, "wt")
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
                "truth_dosage",
                "query_dosage",
                "truth_class",
                "query_class",
                "exact_match",
                "alt_presence_match",
                "het_status_match",
            ],
            delimiter="\t",
        )
        sites_writer.writeheader()

    try:
        for rec in query_vcf_obj.fetch():
            counters["query_records_total"] += 1

            if require_pass and not pass_filter(rec):
                counters["query_records_filter_fail"] += 1
                continue
            if snps_only and not is_snp_record(rec):
                counters["query_records_non_snp"] += 1
                continue
            if biallelic_only and not is_biallelic_record(rec):
                counters["query_records_non_biallelic"] += 1
                continue

            key = make_site_key(rec)

            if key not in truth_index:
                counters["query_records_not_in_truth"] += 1
                continue

            sample_data = rec.samples[sample]

            if not sample_passes_thresholds(
                    sample_data,
                    min_gq=min_query_gq,
                    min_dp=min_query_dp,
            ):
                counters["query_records_sample_threshold_fail"] += 1
                continue

            query_gt = normalize_gt(sample_data.get("GT"))
            if query_gt is None:
                counters["query_gt_missing"] += 1
                continue

            truth_gt, truth_cls = truth_index[key]
            query_cls = gt_class(query_gt)

            truth_has_alt = has_alt(truth_gt)
            query_has_alt = has_alt(query_gt)

            truth_dosage = alt_dosage(truth_gt)
            query_dosage = alt_dosage(query_gt)

            if truth_dosage is not None and query_dosage is not None:
                truth_dosages.append(float(truth_dosage))
                query_dosages.append(float(query_dosage))

            counters["n_compared"] += 1

            exact = truth_gt == query_gt
            alt_match = truth_has_alt == query_has_alt
            het_match = is_het(truth_gt) == is_het(query_gt)

            if exact:
                counters["n_exact_match"] += 1
            if alt_match:
                counters["n_alt_presence_match"] += 1
            if het_match:
                counters["n_het_status_match"] += 1

            # Binary ALT-presence confusion:
            # positive = truth has any non-reference allele.
            if truth_has_alt is True and query_has_alt is True:
                counters["binary_tp"] += 1
            elif truth_has_alt is False and query_has_alt is True:
                counters["binary_fp"] += 1
            elif truth_has_alt is False and query_has_alt is False:
                counters["binary_tn"] += 1
            elif truth_has_alt is True and query_has_alt is False:
                counters["binary_fn"] += 1

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

            if truth_has_alt:
                counters["truth_nonref"] += 1
                if exact:
                    counters["truth_nonref_exact"] += 1
                if query_has_alt:
                    counters["truth_nonref_called_nonref"] += 1

            confusion[(truth_cls, query_cls)] += 1

            if sites_writer is not None:
                sites_writer.writerow({
                    "sample": sample,
                    "coverage": coverage,
                    "replicate": replicate,
                    "method": method,
                    "chrom": rec.chrom,
                    "pos": rec.pos,
                    "ref": rec.ref,
                    "alt": ",".join(rec.alts or []),
                    "truth_gt": "/".join(map(str, truth_gt)),
                    "query_gt": "/".join(map(str, query_gt)),
                    "truth_dosage": truth_dosage,
                    "query_dosage": query_dosage,
                    "truth_class": truth_cls,
                    "query_class": query_cls,
                    "exact_match": int(exact),
                    "alt_presence_match": int(alt_match),
                    "het_status_match": int(het_match),
                })
    finally:
        query_vcf_obj.close()
        if sites_out is not None:
            sites_out.close()

    n = counters["n_compared"]

    tp = counters["binary_tp"]
    fp = counters["binary_fp"]
    tn = counters["binary_tn"]
    fn = counters["binary_fn"]

    precision = rate(tp, tp + fp)
    recall = rate(tp, tp + fn)
    specificity = rate(tn, tn + fp)
    f1 = f1_score(precision, recall)

    dosage_r = pearson_r(truth_dosages, query_dosages)
    dosage_r2 = pearson_r2(truth_dosages, query_dosages)

    metrics = {
        "sample":
        sample,
        "truth_sample":
        truth_sample,
        "coverage":
        coverage,
        "replicate":
        replicate,
        "method":
        method,
        "truth_vcf":
        str(truth_vcf),
        "query_vcf":
        str(query_vcf),
        "n_truth_index_sites":
        len(truth_index),
        "n_query_records_total":
        counters["query_records_total"],
        "n_compared":
        n,

        # Existing concordance metrics.
        "n_exact_match":
        counters["n_exact_match"],
        "exact_match_rate":
        rate(counters["n_exact_match"], n),
        "n_alt_presence_match":
        counters["n_alt_presence_match"],
        "alt_presence_match_rate":
        rate(counters["n_alt_presence_match"], n),
        "n_het_status_match":
        counters["n_het_status_match"],
        "het_status_match_rate":
        rate(counters["n_het_status_match"], n),

        # Class-specific metrics.
        "truth_hom_ref":
        counters["truth_hom_ref"],
        "hom_ref_match_rate":
        rate(
            counters["truth_hom_ref_called_hom_ref"],
            counters["truth_hom_ref"],
        ),
        "truth_het":
        counters["truth_het"],
        "het_recall":
        rate(
            counters["truth_het_called_het"],
            counters["truth_het"],
        ),
        "truth_hom_alt":
        counters["truth_hom_alt"],
        "hom_alt_recall":
        rate(
            counters["truth_hom_alt_called_hom_alt"],
            counters["truth_hom_alt"],
        ),
        "truth_nonref":
        counters["truth_nonref"],
        "nonref_exact_match_rate":
        rate(
            counters["truth_nonref_exact"],
            counters["truth_nonref"],
        ),
        "alt_recall":
        rate(
            counters["truth_nonref_called_nonref"],
            counters["truth_nonref"],
        ),

        # New binary ALT-presence classification metrics.
        "binary_tp":
        tp,
        "binary_fp":
        fp,
        "binary_tn":
        tn,
        "binary_fn":
        fn,
        "precision":
        precision,
        "positive_predictive_value":
        precision,
        "recall":
        recall,
        "sensitivity":
        recall,
        "specificity":
        specificity,
        "f1_score":
        f1,

        # New dosage correlation metrics.
        "n_dosage_pairs":
        len(truth_dosages),
        "dosage_pearson_r":
        dosage_r,
        "dosage_r2":
        dosage_r2,

        # Diagnostics.
        "query_records_not_in_truth":
        counters["query_records_not_in_truth"],
        "query_gt_missing":
        counters["query_gt_missing"],
    }

    metrics_path = Path(str(output_prefix) + ".metrics.tsv")
    with metrics_path.open("w", newline="") as f:
        writer = csv.DictWriter(f,
                                fieldnames=list(metrics.keys()),
                                delimiter="\t")
        writer.writeheader()
        writer.writerow(metrics)

    confusion_path = Path(str(output_prefix) + ".confusion.tsv")
    with confusion_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "sample",
                "coverage",
                "replicate",
                "method",
                "truth_class",
                "query_class",
                "count",
            ],
            delimiter="\t",
        )
        writer.writeheader()

        for (truth_cls, query_cls), count in sorted(confusion.items()):
            writer.writerow({
                "sample": sample,
                "coverage": coverage,
                "replicate": replicate,
                "method": method,
                "truth_class": truth_cls,
                "query_class": query_cls,
                "count": count,
            })

    binary_confusion_path = Path(str(output_prefix) + ".binary_confusion.tsv")
    with binary_confusion_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "sample", "coverage", "replicate", "method", "class", "count"
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for cls in ["binary_tp", "binary_fp", "binary_tn", "binary_fn"]:
            writer.writerow({
                "sample": sample,
                "coverage": coverage,
                "replicate": replicate,
                "method": method,
                "class": cls,
                "count": counters[cls],
            })

    json_path = Path(str(output_prefix) + ".metrics.json")
    with json_path.open("w") as f:
        json.dump(metrics, f, indent=2)

    return metrics
