import argparse
from pathlib import Path

from lowseq_runner.config import load_config
from lowseq_runner.experiments import run_experiments
from lowseq_runner.metrics import run_metrics
from lowseq_runner.orchestrator import run_all
from lowseq_runner.report import generate_report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Low-pass SNP calling experiment runner.", )

    parser.add_argument(
        "command",
        choices=[
            "all",
            "experiments",
            "metrics",
            "report",
        ],
        help="What to run.",
    )

    parser.add_argument(
        "--env-file",
        default=".env",
        help="Path to .env file.",
    )

    parser.add_argument(
        "--manifest",
        default=None,
        help="Optional manifest for metrics command.",
    )

    parser.add_argument(
        "--metrics-dir",
        default=None,
        help="Optional metrics directory for report command.",
    )

    parser.add_argument(
        "--report-dir",
        default=None,
        help="Optional report output directory for report command.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.env_file)

    if args.command == "all":
        run_all(config)

    elif args.command == "experiments":
        run_experiments(config)

    elif args.command == "metrics":
        manifest = Path(args.manifest) if args.manifest else config.manifest
        run_metrics(config, manifest=manifest)

    elif args.command == "report":
        metrics_dir = Path(
            args.metrics_dir) if args.metrics_dir else config.metrics_dir
        report_dir = Path(
            args.report_dir) if args.report_dir else config.report_dir

        generate_report(
            metrics_dir=metrics_dir,
            output_dir=report_dir,
            title_prefix=config.report_title_prefix,
        )

    else:
        raise ValueError(args.command)


if __name__ == "__main__":
    main()
