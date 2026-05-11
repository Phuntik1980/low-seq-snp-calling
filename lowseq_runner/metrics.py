import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from lowseq_runner.config import LowseqConfig
from lowseq_runner.genotype_metrics import compute_genotype_metrics


def get_query_vcf(output_dir: Path, sample: str, vcf_kind: str) -> Path:
    if vcf_kind == "imputed":
        return output_dir / "8_imputed" / f"{sample}.imputed.vcf.gz"

    if vcf_kind == "bcftools":
        return output_dir / "7_vcf" / f"{sample}.bcftools.vcf.gz"

    raise ValueError(f"Unsupported vcf_kind: {vcf_kind}")


def safe_cov_label_from_manifest(row: dict) -> str:
    return row.get("coverage_label") or f"{float(row['coverage']):g}x".replace(
        ".", "p")


def build_metric_jobs(config: LowseqConfig,
                      manifest: Path | None = None) -> list[dict]:
    manifest = manifest or config.manifest
    jobs: list[dict] = []

    with manifest.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            sample = row["sample"]
            coverage = float(row["coverage"])
            coverage_label = safe_cov_label_from_manifest(row)
            replicate = int(row["replicate"])
            run_output_dir = Path(row["output_dir"])

            truth_vcf = config.truth_dir / config.truth_template.format(
                sample=sample)
            query_vcf = get_query_vcf(run_output_dir, sample, config.vcf_kind)

            if not truth_vcf.exists():
                print(f"WARN: truth VCF not found, skip: {truth_vcf}")
                continue

            if not query_vcf.exists():
                print(f"WARN: query VCF not found, skip: {query_vcf}")
                continue

            prefix = (
                config.metrics_dir /
                f"{sample}_cov_{coverage_label}_rep_{replicate}_{config.method}"
            )

            jobs.append({
                "truth_vcf": truth_vcf,
                "query_vcf": query_vcf,
                "sample": sample,
                "coverage": coverage,
                "replicate": replicate,
                "method": config.method,
                "output_prefix": prefix,
            })

    return jobs


def _run_metric_job(config: LowseqConfig, job: dict) -> dict:
    return compute_genotype_metrics(
        truth_vcf=job["truth_vcf"],
        query_vcf=job["query_vcf"],
        sample=job["sample"],
        coverage=job["coverage"],
        replicate=job["replicate"],
        method=job["method"],
        output_prefix=job["output_prefix"],
        snps_only=config.snps_only,
        biallelic_only=config.biallelic_only,
        require_pass=config.require_pass,
        min_truth_gq=config.min_truth_gq,
        min_query_gq=config.min_query_gq,
        min_truth_dp=config.min_truth_dp,
        min_query_dp=config.min_query_dp,
        write_sites=config.write_sites,
    )


def run_metrics(config: LowseqConfig, manifest: Path | None = None) -> Path:
    config.metrics_dir.mkdir(parents=True, exist_ok=True)
    jobs = build_metric_jobs(config, manifest)

    print(f"Total metric jobs: {len(jobs)}")
    print(f"Parallel metric jobs: {config.metrics_parallel_jobs}")

    if config.dry_run:
        for job in jobs:
            print(job)
        return config.metrics_dir

    failed: list[tuple[dict, Exception]] = []

    with ThreadPoolExecutor(
            max_workers=config.metrics_parallel_jobs) as executor:
        future_to_job = {
            executor.submit(_run_metric_job, config, job): job
            for job in jobs
        }

        for future in as_completed(future_to_job):
            job = future_to_job[future]
            try:
                future.result()
                print(
                    "[METRICS DONE] "
                    f"{job['sample']} cov={job['coverage']} rep={job['replicate']}"
                )
            except Exception as exc:
                print(
                    "[METRICS FAIL] "
                    f"{job['sample']} cov={job['coverage']} rep={job['replicate']}: {exc}"
                )
                failed.append((job, exc))

    if failed:
        raise RuntimeError(f"{len(failed)} metric jobs failed")

    return config.metrics_dir
