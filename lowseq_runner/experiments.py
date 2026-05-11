import csv
import json
from concurrent.futures import FIRST_COMPLETED, Future, wait
from pathlib import Path
from typing import Optional

from lowseq_runner.config import LowseqConfig
from lowseq_runner.shell import quote_cmd, run_to_log


def read_genome_size_from_fai(reference: Path) -> int:
    fai = Path(str(reference) + ".fai")
    if not fai.exists():
        raise FileNotFoundError(
            f"Genome size was not provided and FASTA index was not found: {fai}. "
            f"Run: samtools faidx {reference}")

    total = 0
    with fai.open() as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            total += int(fields[1])
    return total


def coverage_to_reads_per_file(coverage: float, genome_size: int,
                               read_length: int) -> int:
    reads = coverage * genome_size / (2.0 * read_length)
    return max(1, int(round(reads)))


def safe_cov_label(cov: float) -> str:
    return f"{str(f'{cov:g}').replace('.', 'p')}x"


def make_seed(seed_base: int, sample_index: int, coverage_index: int,
              replicate_index: int) -> int:
    return seed_base + sample_index * 10000 + coverage_index * 100 + replicate_index


def build_pipeline_command(
    pipeline: Path,
    input_dir: Path,
    output_dir: Path,
    sample: str,
    sample_reads: int,
    seed: int,
    reference: Path,
    glimpse_panel: Optional[Path],
    threads: int,
    skip_varscan: bool,
    skip_glimpse: bool,
    skip_qualimap: bool,
) -> list[str]:
    cmd = [
        "bash",
        str(pipeline),
        "--input-dir",
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--sample",
        sample,
        "--sample-reads",
        str(sample_reads),
        "--seed",
        str(seed),
        "--reference",
        str(reference),
        "--threads",
        str(threads),
    ]

    if glimpse_panel:
        cmd.extend(["--glimpse-panel", str(glimpse_panel)])

    if skip_varscan:
        cmd.append("--skip-varscan")
    if skip_glimpse:
        cmd.append("--skip-glimpse")
    if skip_qualimap:
        cmd.append("--skip-qualimap")

    return cmd


def build_experiment_plan(config: LowseqConfig) -> list[dict]:
    genome_size = config.genome_size or read_genome_size_from_fai(
        config.reference)
    rows: list[dict] = []

    for sample_idx, sample in enumerate(config.samples):
        input_dir = config.input_root if config.input_mode == "root" else config.input_root / sample

        for cov_idx, cov in enumerate(config.coverages):
            cov_label = safe_cov_label(cov)
            sample_reads = coverage_to_reads_per_file(cov, genome_size,
                                                      config.read_length)

            for rep in range(1, config.replicates + 1):
                seed = make_seed(config.seed_base, sample_idx, cov_idx, rep)
                run_dir = config.output_root / sample / f"cov_{cov_label}" / f"rep_{rep}"
                log_path = run_dir / "pipeline.log"

                cmd = build_pipeline_command(
                    pipeline=config.pipeline,
                    input_dir=input_dir,
                    output_dir=run_dir,
                    sample=sample,
                    sample_reads=sample_reads,
                    seed=seed,
                    reference=config.reference,
                    glimpse_panel=config.glimpse_panel,
                    threads=config.threads_per_run,
                    skip_varscan=config.skip_varscan,
                    skip_glimpse=config.skip_glimpse,
                    skip_qualimap=config.skip_qualimap,
                )

                rows.append({
                    "sample": sample,
                    "coverage": cov,
                    "coverage_label": cov_label,
                    "replicate": rep,
                    "seed": seed,
                    "sample_reads_per_file": sample_reads,
                    "input_dir": str(input_dir),
                    "output_dir": str(run_dir),
                    "log": str(log_path),
                    "command": quote_cmd(cmd),
                    "_cmd": cmd,
                    "_log": log_path,
                })

    return rows


def write_manifest(config: LowseqConfig, rows: list[dict]) -> Path:
    config.output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = config.manifest

    public_rows = [{
        k: v
        for k, v in row.items() if not k.startswith("_")
    } for row in rows]

    with manifest_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(public_rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(public_rows)

    return manifest_path


def write_launch_config(config: LowseqConfig, genome_size: int) -> Path:
    path = config.output_root / "launch_config.json"

    with path.open("w") as f:
        json.dump(
            {
                "pipeline":
                str(config.pipeline),
                "input_root":
                str(config.input_root),
                "output_root":
                str(config.output_root),
                "reference":
                str(config.reference),
                "glimpse_panel":
                str(config.glimpse_panel) if config.glimpse_panel else None,
                "samples":
                config.samples,
                "coverages":
                config.coverages,
                "replicates":
                config.replicates,
                "genome_size":
                genome_size,
                "read_length":
                config.read_length,
                "threads_per_run":
                config.threads_per_run,
                "parallel_runs":
                config.parallel_runs,
                "seed_base":
                config.seed_base,
            },
            f,
            indent=2,
        )

    return path


def run_experiments(config: LowseqConfig) -> Path:
    genome_size = config.genome_size or read_genome_size_from_fai(
        config.reference)
    rows = build_experiment_plan(config)

    manifest_path = write_manifest(config, rows)
    config_path = write_launch_config(config, genome_size)

    print(f"Manifest: {manifest_path}")
    print(f"Config:   {config_path}")
    print(f"Total runs: {len(rows)}")
    print(f"Parallel runs: {config.parallel_runs}")
    print(f"Threads per run: {config.threads_per_run}")

    if config.dry_run:
        for row in rows:
            print(row["command"])
        return manifest_path

    failed: list[tuple[int, dict]] = []
    pending_rows = list(rows)
    running: dict[Future[int], dict] = {}

    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=config.parallel_runs) as executor:
        while pending_rows or running:
            while pending_rows and len(running) < config.parallel_runs:
                row = pending_rows.pop(0)
                print(f"[START] {row['command']}")
                print(f"[LOG]   {row['_log']}")
                fut = executor.submit(run_to_log, row["_cmd"], row["_log"])
                running[fut] = row

            done, _ = wait(running.keys(), return_when=FIRST_COMPLETED)

            for fut in done:
                row = running.pop(fut)
                ret = fut.result()

                if ret == 0:
                    print(f"[DONE]  {row['_log']}")
                else:
                    print(f"[FAIL]  exit code {ret}: {row['_log']}")
                    failed.append((ret, row))

    if failed:
        message = "\n".join(
            f"exit={ret}, log={row['_log']}, cmd={row['command']}"
            for ret, row in failed)
        raise RuntimeError(f"{len(failed)} pipeline runs failed:\n{message}")

    return manifest_path
