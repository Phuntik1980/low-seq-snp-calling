import argparse
import csv
import json
import math
import os
import shlex
import subprocess
from pathlib import Path
from typing import Dict, List, Optional


def parse_args():
    parser = argparse.ArgumentParser(
        description="Launch low-pass sequencing pipeline for samples × coverages × replicates."
    )

    parser.add_argument(
        "--pipeline",
        required=True,
        help="Path to parameterized bash pipeline, e.g. pipeline_param.sh",
    )
    parser.add_argument(
        "--input-root",
        required=True,
        help="Root directory with FASTQ files or per-sample FASTQ directories.",
    )
    parser.add_argument(
        "--output-root",
        required=True,
        help="Root output directory for all runs.",
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA.",
    )
    parser.add_argument(
        "--glimpse-panel",
        default=None,
        help="Optional GLIMPSE2 reference panel VCF.gz.",
    )
    parser.add_argument(
        "--samples",
        required=True,
        nargs="+",
        help="Sample names, e.g. 208_S15 209_S16 210_S17.",
    )
    parser.add_argument(
        "--coverages",
        nargs="+",
        type=float,
        default=[0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0],
        help="Target coverages.",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=3,
        help="Number of downsampling replicates per sample and coverage.",
    )
    parser.add_argument(
        "--genome-size",
        type=int,
        default=None,
        help="Genome size in bp. If not provided, will be calculated from reference.fai if available.",
    )
    parser.add_argument(
        "--read-length",
        type=int,
        default=150,
        help="Read length after sequencing. Used for reads calculation.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=os.cpu_count() or 1,
        help="Threads per pipeline run.",
    )
    parser.add_argument(
        "--parallel-runs",
        type=int,
        default=1,
        help="How many pipeline runs to execute in parallel.",
    )
    parser.add_argument(
        "--seed-base",
        type=int,
        default=1000,
        help="Base seed. Actual seeds are deterministically generated from sample, coverage, replicate.",
    )
    parser.add_argument(
        "--input-mode",
        choices=["root", "per-sample"],
        default="root",
        help=(
            "root: FASTQs are directly in input-root. "
            "per-sample: FASTQs are in input-root/sample/."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands and manifest but do not run.",
    )
    parser.add_argument(
        "--skip-varscan",
        action="store_true",
        help="Pass --skip-varscan to pipeline.",
    )
    parser.add_argument(
        "--skip-glimpse",
        action="store_true",
        help="Pass --skip-glimpse to pipeline.",
    )
    parser.add_argument(
        "--skip-qualimap",
        action="store_true",
        help="Pass --skip-qualimap to pipeline.",
    )

    return parser.parse_args()


def read_genome_size_from_fai(reference: Path) -> int:
    fai = Path(str(reference) + ".fai")
    if not fai.exists():
        raise FileNotFoundError(
            f"Genome size was not provided and FASTA index was not found: {fai}. "
            f"Run: samtools faidx {reference}"
        )

    total = 0
    with fai.open() as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            total += int(fields[1])
    return total


def coverage_to_reads_per_file(coverage: float, genome_size: int, read_length: int) -> int:
    reads = coverage * genome_size / (2.0 * read_length)
    return max(1, int(round(reads)))


def safe_cov_label(cov: float) -> str:
    text = f"{cov:g}".replace(".", "p")
    return f"{text}x"


def make_seed(seed_base: int, sample_index: int, coverage_index: int, replicate_index: int) -> int:
    return seed_base + sample_index * 10000 + coverage_index * 100 + replicate_index


def build_command(
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
) -> List[str]:

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

    if glimpse_panel is not None:
        cmd.extend(["--glimpse-panel", str(glimpse_panel)])

    if skip_varscan:
        cmd.append("--skip-varscan")

    if skip_glimpse:
        cmd.append("--skip-glimpse")

    if skip_qualimap:
        cmd.append("--skip-qualimap")

    return cmd


def run_commands(commands: List[List[str]], logs: List[Path], parallel_runs: int, dry_run: bool):
    if dry_run:
        for cmd in commands:
            print(" ".join(shlex.quote(x) for x in cmd))
        return

    running = []
    command_iter = iter(zip(commands, logs))

    def start_next():
        try:
            cmd, log_path = next(command_iter)
        except StopIteration:
            return False

        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_file = log_path.open("w")

        print(f"[START] {' '.join(shlex.quote(x) for x in cmd)}")
        print(f"[LOG]   {log_path}")

        proc = subprocess.Popen(
            cmd,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
        )
        running.append((proc, log_file, cmd, log_path))
        return True

    for _ in range(parallel_runs):
        if not start_next():
            break

    failed = []

    while running:
        for item in list(running):
            proc, log_file, cmd, log_path = item
            ret = proc.poll()

            if ret is None:
                continue

            log_file.close()
            running.remove(item)

            if ret == 0:
                print(f"[DONE]  {log_path}")
            else:
                print(f"[FAIL]  exit code {ret}: {log_path}")
                failed.append((ret, cmd, log_path))

            start_next()

    if failed:
        message = "\n".join(
            f"exit={ret}, log={log_path}, cmd={' '.join(shlex.quote(x) for x in cmd)}"
            for ret, cmd, log_path in failed
        )
        raise RuntimeError(f"{len(failed)} pipeline runs failed:\n{message}")


def main():
    args = parse_args()

    pipeline = Path(args.pipeline).resolve()
    input_root = Path(args.input_root).resolve()
    output_root = Path(args.output_root).resolve()
    reference = Path(args.reference).resolve()
    glimpse_panel = Path(args.glimpse_panel).resolve() if args.glimpse_panel else None

    if args.genome_size is not None:
        genome_size = args.genome_size
    else:
        genome_size = read_genome_size_from_fai(reference)

    output_root.mkdir(parents=True, exist_ok=True)

    commands = []
    logs = []
    manifest_rows = []

    for sample_idx, sample in enumerate(args.samples):
        if args.input_mode == "root":
            input_dir = input_root
        else:
            input_dir = input_root / sample

        for cov_idx, cov in enumerate(args.coverages):
            cov_label = safe_cov_label(cov)
            sample_reads = coverage_to_reads_per_file(
                coverage=cov,
                genome_size=genome_size,
                read_length=args.read_length,
            )

            for rep in range(1, args.replicates + 1):
                seed = make_seed(
                    seed_base=args.seed_base,
                    sample_index=sample_idx,
                    coverage_index=cov_idx,
                    replicate_index=rep,
                )

                run_dir = output_root / sample / f"cov_{cov_label}" / f"rep_{rep}"
                log_path = run_dir / "pipeline.log"

                cmd = build_command(
                    pipeline=pipeline,
                    input_dir=input_dir,
                    output_dir=run_dir,
                    sample=sample,
                    sample_reads=sample_reads,
                    seed=seed,
                    reference=reference,
                    glimpse_panel=glimpse_panel,
                    threads=args.threads,
                    skip_varscan=args.skip_varscan,
                    skip_glimpse=args.skip_glimpse,
                    skip_qualimap=args.skip_qualimap,
                )

                commands.append(cmd)
                logs.append(log_path)

                manifest_rows.append(
                    {
                        "sample": sample,
                        "coverage": cov,
                        "coverage_label": cov_label,
                        "replicate": rep,
                        "seed": seed,
                        "sample_reads_per_file": sample_reads,
                        "input_dir": str(input_dir),
                        "output_dir": str(run_dir),
                        "log": str(log_path),
                        "command": " ".join(shlex.quote(x) for x in cmd),
                    }
                )

    manifest_tsv = output_root / "launch_manifest.tsv"
    with manifest_tsv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(manifest_rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    config_json = output_root / "launch_config.json"
    with config_json.open("w") as f:
        json.dump(
            {
                "pipeline": str(pipeline),
                "input_root": str(input_root),
                "output_root": str(output_root),
                "reference": str(reference),
                "glimpse_panel": str(glimpse_panel) if glimpse_panel else None,
                "samples": args.samples,
                "coverages": args.coverages,
                "replicates": args.replicates,
                "genome_size": genome_size,
                "read_length": args.read_length,
                "threads": args.threads,
                "parallel_runs": args.parallel_runs,
                "seed_base": args.seed_base,
            },
            f,
            indent=2,
        )

    print(f"Manifest: {manifest_tsv}")
    print(f"Config:   {config_json}")
    print(f"Total runs: {len(commands)}")

    run_commands(
        commands=commands,
        logs=logs,
        parallel_runs=args.parallel_runs,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()