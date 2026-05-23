import os
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv
from pydantic import BaseModel, Field


def _bool(value: str | bool | None, default: bool = False) -> bool:
    if value is None or value == "":
        return default
    if isinstance(value, bool):
        return value
    return value.strip().lower() in {"1", "true", "yes", "y", "on"}


def _csv_str(value: str | None) -> list[str]:
    if not value:
        return []
    return [x.strip() for x in value.split(",") if x.strip()]


def _csv_float(value: str | None) -> list[float]:
    return [float(x) for x in _csv_str(value)]


def _optional_int(value: str | None) -> Optional[int]:
    if value is None or value.strip() == "":
        return None
    return int(value)


def _optional_float(value: str | None) -> Optional[float]:
    if value is None or value.strip() == "":
        return None
    return float(value)


class LowseqConfig(BaseModel):
    pipeline: Path
    input_root: Path
    output_root: Path
    reference: Path
    glimpse_panel: Optional[Path] = None
    glimpse_sites_vcf: Optional[Path] = None
    glimpse_sites_tsv: Optional[Path] = None
    glimpse_chunks: Optional[Path] = None

    truth_dir: Path
    truth_template: str = "{sample}.truth.vcf.gz"

    samples: list[str]
    coverages: list[float]
    replicates: int = 3
    read_length: int = 150
    genome_size: Optional[int] = None
    input_mode: str = "root"

    parallel_runs: int = 1
    threads_per_run: int = 1
    metrics_parallel_jobs: int = 1

    skip_varscan: bool = True
    skip_glimpse: bool = False
    skip_qualimap: bool = True

    vcf_kind: str = "imputed"
    method: str = "glimpse2"
    snps_only: bool = True
    biallelic_only: bool = True
    require_pass: bool = False

    min_truth_gq: Optional[float] = None
    min_query_gq: Optional[float] = None
    min_truth_dp: Optional[float] = None
    min_query_dp: Optional[float] = None

    write_sites: bool = False

    report_title_prefix: str = "Low-pass SNP genotyping"
    seed_base: int = 1000
    dry_run: bool = False

    manifest_path: Path = Field(default=Path("launch_manifest.tsv"))

    @property
    def metrics_dir(self) -> Path:
        return self.output_root / "genotype_metrics"

    @property
    def report_dir(self) -> Path:
        return self.output_root / "report"

    @property
    def manifest(self) -> Path:
        return self.output_root / "launch_manifest.tsv"


def load_config(env_file: str | Path = ".env") -> LowseqConfig:
    load_dotenv(env_file)

    return LowseqConfig(
        pipeline=Path(os.environ["LOWSEQ_PIPELINE"]),
        input_root=Path(os.environ["LOWSEQ_INPUT_ROOT"]),
        output_root=Path(os.environ["LOWSEQ_OUTPUT_ROOT"]),
        reference=Path(os.environ["LOWSEQ_REFERENCE"]),
        glimpse_panel=(Path(os.environ["LOWSEQ_GLIMPSE_PANEL"])
                       if os.getenv("LOWSEQ_GLIMPSE_PANEL") else None),
        glimpse_sites_vcf=(Path(os.environ["LOWSEQ_GLIMPSE_SITES_VCF"])
                           if os.getenv("LOWSEQ_GLIMPSE_SITES_VCF") else None),
        glimpse_sites_tsv=(Path(os.environ["LOWSEQ_GLIMPSE_SITES_TSV"])
                           if os.getenv("LOWSEQ_GLIMPSE_SITES_TSV") else None),
        glimpse_chunks=(Path(os.environ["LOWSEQ_GLIMPSE_CHUNKS"])
                        if os.getenv("LOWSEQ_GLIMPSE_CHUNKS") else None),
        truth_dir=Path(os.environ["LOWSEQ_TRUTH_DIR"]),
        truth_template=os.getenv("LOWSEQ_TRUTH_TEMPLATE"),
        samples=_csv_str(os.getenv("LOWSEQ_SAMPLES")),
        coverages=_csv_float(os.getenv("LOWSEQ_COVERAGES")),
        replicates=int(os.getenv("LOWSEQ_REPLICATES", "3")),
        read_length=int(os.getenv("LOWSEQ_READ_LENGTH", "150")),
        genome_size=_optional_int(os.getenv("LOWSEQ_GENOME_SIZE")),
        input_mode=os.getenv("LOWSEQ_INPUT_MODE", "root"),
        parallel_runs=int(os.getenv("LOWSEQ_PARALLEL_RUNS", "1")),
        threads_per_run=int(os.getenv("LOWSEQ_THREADS_PER_RUN", "1")),
        metrics_parallel_jobs=int(os.getenv("LOWSEQ_METRICS_PARALLEL_JOBS",
                                            "1")),
        skip_varscan=_bool(os.getenv("LOWSEQ_SKIP_VARSCAN"), True),
        skip_glimpse=_bool(os.getenv("LOWSEQ_SKIP_GLIMPSE"), False),
        skip_qualimap=_bool(os.getenv("LOWSEQ_SKIP_QUALIMAP"), True),
        vcf_kind=os.getenv("LOWSEQ_VCF_KIND", "imputed"),
        method=os.getenv("LOWSEQ_METHOD", "glimpse2"),
        snps_only=_bool(os.getenv("LOWSEQ_SNPS_ONLY"), True),
        biallelic_only=_bool(os.getenv("LOWSEQ_BIALLELIC_ONLY"), True),
        require_pass=_bool(os.getenv("LOWSEQ_REQUIRE_PASS"), False),
        min_truth_gq=_optional_float(os.getenv("LOWSEQ_MIN_TRUTH_GQ")),
        min_query_gq=_optional_float(os.getenv("LOWSEQ_MIN_QUERY_GQ")),
        min_truth_dp=_optional_float(os.getenv("LOWSEQ_MIN_TRUTH_DP")),
        min_query_dp=_optional_float(os.getenv("LOWSEQ_MIN_QUERY_DP")),
        write_sites=_bool(os.getenv("LOWSEQ_WRITE_SITES"), False),
        report_title_prefix=os.getenv("LOWSEQ_REPORT_TITLE_PREFIX"),
        seed_base=int(os.getenv("LOWSEQ_SEED_BASE", "1000")),
        dry_run=_bool(os.getenv("LOWSEQ_DRY_RUN"), False),
    )
