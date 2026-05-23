from lowseq_runner.config import LowseqConfig
from lowseq_runner.experiments import run_experiments
from lowseq_runner.metrics import run_metrics
from lowseq_runner.report import generate_report


def run_all(config: LowseqConfig) -> None:
    manifest = run_experiments(config)
    run_metrics(config, manifest=manifest)
    generate_report(
        metrics_dir=config.metrics_dir,
        output_dir=config.report_dir,
        title_prefix=config.report_title_prefix,
    )
