#!/usr/bin/env python3

import argparse
import csv
import shlex
import subprocess
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Launch genotype metrics for all runs from launch_manifest.tsv."
    )

    parser.add_argument(
        "--manifest",
        required=True,
        help="launch_manifest.tsv from launch_experiments.py",
    )
    parser.add_argument(
        "--metrics-script",
        required=True,
        help="Path to compute_genotype_metrics.py",
    )
    parser.add_argument(
        "--truth-dir",
        required=True,
        help="Directory with truth VCFs.",
    )
    parser.add_argument(
        "--truth-template",
        default="{sample}.truth.vcf.gz",
        help="Truth VCF filename template relative to truth-dir.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for metrics output.",
    )
    parser.add_argument(
        "--vcf-kind",
        choices=["imputed", "bcftools"],
        default="imputed",
        help="Which query VCF to evaluate.",
    )
    parser.add_argument(
        "--method",
        default=None,
        help="Method label. Default: same as --vcf-kind.",
    )
    parser.add_argument(
        "--snps-only",
        action="store_true",
    )
    parser.add_argument(
        "--biallelic-only",
        action="store_true",
    )
    parser.add_argument(
        "--require-pass",
        action="store_true",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
    )

    return parser.parse_args()


def get_query_vcf(output_dir: Path, sample: str, vcf_kind: str) -> Path:
    if vcf_kind == "imputed":
        return output_dir / "9_imputed" / f"{sample}.imputed.vcf.gz"
    if vcf_kind == "bcftools":
        return output_dir / "8_vcf" / f"{sample}.bcftools.vcf.gz"
    raise ValueError(vcf_kind)


def safe_cov_label_from_manifest(row):
    return row.get("coverage_label") or f"{float(row['coverage']):g}x".replace(".", "p")


def main():
    args = parse_args()

    manifest = Path(args.manifest)
    metrics_script = Path(args.metrics_script).resolve()
    truth_dir = Path(args.truth_dir).resolve()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    method = args.method or args.vcf_kind

    commands = []

    with manifest.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            coverage = float(row["coverage"])
            coverage_label = safe_cov_label_from_manifest(row)
            replicate = int(row["replicate"])
            run_output_dir = Path(row["output_dir"])

            truth_vcf = truth_dir / args.truth_template.format(sample=sample)
            query_vcf = get_query_vcf(run_output_dir, sample, args.vcf_kind)

            if not truth_vcf.exists():
                print(f"WARN: truth VCF not found, skip: {truth_vcf}")
                continue

            if not query_vcf.exists():
                print(f"WARN: query VCF not found, skip: {query_vcf}")
                continue

            prefix = out_dir / f"{sample}_cov_{coverage_label}_rep_{replicate}_{method}"

            cmd = [
                "python",
                str(metrics_script),
                "--truth-vcf",
                str(truth_vcf),
                "--query-vcf",
                str(query_vcf),
                "--sample",
                sample,
                "--coverage",
                str(coverage),
                "--replicate",
                str(replicate),
                "--method",
                method,
                "--output-prefix",
                str(prefix),
            ]

            if args.snps_only:
                cmd.append("--snps-only")
            if args.biallelic_only:
                cmd.append("--biallelic-only")
            if args.require_pass:
                cmd.append("--require-pass")

            commands.append(cmd)

    print(f"Total metric jobs: {len(commands)}")

    for cmd in commands:
        print(" ".join(shlex.quote(x) for x in cmd))

        if not args.dry_run:
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()