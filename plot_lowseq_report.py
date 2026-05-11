#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create plots and summary report for low-seq genotype concordance metrics."
    )

    parser.add_argument(
        "--metrics-dir",
        required=True,
        help="Directory with *.metrics.tsv files.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for plots and summary tables.",
    )
    parser.add_argument(
        "--title-prefix",
        default="Low-pass sequencing",
        help="Prefix for plot titles.",
    )

    return parser.parse_args()


def read_metrics(metrics_dir: Path) -> pd.DataFrame:
    files = sorted(metrics_dir.glob("*.metrics.tsv"))

    if not files:
        raise FileNotFoundError(f"No *.metrics.tsv files found in {metrics_dir}")

    frames = []
    for path in files:
        df = pd.read_csv(path, sep="\t")
        df["metrics_file"] = str(path)
        frames.append(df)

    data = pd.concat(frames, ignore_index=True)

    numeric_cols = [
        "coverage",
        "replicate",
        "n_compared",
        "exact_match_rate",
        "alt_presence_match_rate",
        "het_status_match_rate",
        "hom_ref_match_rate",
        "het_recall",
        "hom_alt_recall",
        "nonref_exact_match_rate",
        "alt_recall",
    ]

    for col in numeric_cols:
        if col in data.columns:
            data[col] = pd.to_numeric(data[col], errors="coerce")

    return data


def summarize(data: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "n_compared",
        "exact_match_rate",
        "alt_presence_match_rate",
        "het_status_match_rate",
        "hom_ref_match_rate",
        "het_recall",
        "hom_alt_recall",
        "nonref_exact_match_rate",
        "alt_recall",
    ]

    available = [m for m in metrics if m in data.columns]

    summary = (
        data
        .groupby(["sample", "method", "coverage"], as_index=False)
        .agg(
            **{
                f"{m}_mean": (m, "mean")
                for m in available
            },
            **{
                f"{m}_sd": (m, "std")
                for m in available
            },
            n_replicates=("replicate", "nunique"),
        )
    )

    return summary


def save_lineplot(
    data: pd.DataFrame,
    y: str,
    output_path: Path,
    title: str,
    ylabel: str,
):
    plt.figure(figsize=(9, 6))

    sns.lineplot(
        data=data,
        x="coverage",
        y=y,
        hue="sample",
        style="method",
        markers=True,
        dashes=False,
        errorbar="sd",
    )

    plt.title(title)
    plt.xlabel("Target coverage, x")
    plt.ylabel(ylabel)
    plt.ylim(0, 1.02)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def save_pooled_lineplot(
    data: pd.DataFrame,
    y: str,
    output_path: Path,
    title: str,
    ylabel: str,
):
    plt.figure(figsize=(9, 6))

    sns.lineplot(
        data=data,
        x="coverage",
        y=y,
        hue="method",
        markers=True,
        dashes=False,
        errorbar="sd",
    )

    plt.title(title)
    plt.xlabel("Target coverage, x")
    plt.ylabel(ylabel)
    plt.ylim(0, 1.02)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def save_bar_ncompared(data: pd.DataFrame, output_path: Path, title: str):
    plt.figure(figsize=(10, 6))

    sns.lineplot(
        data=data,
        x="coverage",
        y="n_compared",
        hue="sample",
        style="method",
        markers=True,
        dashes=False,
        errorbar="sd",
    )

    plt.title(title)
    plt.xlabel("Target coverage, x")
    plt.ylabel("Number of compared sites")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def make_markdown_report(
    data: pd.DataFrame,
    summary: pd.DataFrame,
    output_dir: Path,
    title_prefix: str,
):
    report_path = output_dir / "report.md"

    metrics_to_show = [
        ("exact_match_rate", "Exact genotype concordance"),
        ("alt_presence_match_rate", "Alternative allele presence concordance"),
        ("het_status_match_rate", "Heterozygosity status concordance"),
        ("het_recall", "WGS heterozygote recall"),
        ("hom_alt_recall", "WGS homozygous alternative recall"),
        ("alt_recall", "WGS non-reference recall"),
    ]

    with report_path.open("w") as f:
        f.write(f"# {title_prefix}: genotype concordance report\n\n")

        f.write("## Input summary\n\n")
        f.write(f"- Number of metric rows: {len(data)}\n")
        f.write(f"- Samples: {', '.join(map(str, sorted(data['sample'].unique())))}\n")
        f.write(f"- Methods: {', '.join(map(str, sorted(data['method'].unique())))}\n")
        f.write(
            f"- Coverages: {', '.join(map(lambda x: f'{x:g}', sorted(data['coverage'].dropna().unique())))}\n\n"
        )

        f.write("## Main plots\n\n")

        for metric, label in metrics_to_show:
            png = f"{metric}.png"
            if (output_dir / png).exists():
                f.write(f"### {label}\n\n")
                f.write(f"![{label}]({png})\n\n")

        if (output_dir / "n_compared.png").exists():
            f.write("### Number of compared sites\n\n")
            f.write("![Number of compared sites](n_compared.png)\n\n")

        f.write("## Summary table\n\n")
        selected_cols = [
            col for col in summary.columns
            if col in [
                "sample",
                "method",
                "coverage",
                "n_replicates",
                "n_compared_mean",
                "exact_match_rate_mean",
                "exact_match_rate_sd",
                "alt_presence_match_rate_mean",
                "alt_presence_match_rate_sd",
                "het_status_match_rate_mean",
                "het_status_match_rate_sd",
                "het_recall_mean",
                "het_recall_sd",
                "hom_alt_recall_mean",
                "hom_alt_recall_sd",
                "alt_recall_mean",
                "alt_recall_sd",
            ]
        ]

        f.write(summary[selected_cols].to_markdown(index=False))
        f.write("\n")

    return report_path


def main():
    args = parse_args()

    metrics_dir = Path(args.metrics_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    sns.set_theme(style="whitegrid")

    data = read_metrics(metrics_dir)
    summary = summarize(data)

    all_metrics_path = output_dir / "all_metrics.tsv"
    summary_path = output_dir / "summary_by_sample_coverage.tsv"

    data.to_csv(all_metrics_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)

    plot_specs = [
        (
            "exact_match_rate",
            "Exact genotype concordance",
            "Fraction of exactly matching genotypes",
        ),
        (
            "alt_presence_match_rate",
            "Alternative allele presence concordance",
            "Fraction with matching ALT presence",
        ),
        (
            "het_status_match_rate",
            "Heterozygosity status concordance",
            "Fraction with matching heterozygosity status",
        ),
        (
            "het_recall",
            "WGS heterozygote recall",
            "Fraction of WGS heterozygotes called as heterozygotes",
        ),
        (
            "hom_alt_recall",
            "WGS homozygous ALT recall",
            "Fraction of WGS hom-alt sites called as hom-alt",
        ),
        (
            "alt_recall",
            "WGS non-reference recall",
            "Fraction of WGS non-ref sites called as non-ref",
        ),
        (
            "nonref_exact_match_rate",
            "Exact concordance on WGS non-reference sites",
            "Exact match among WGS non-ref sites",
        ),
    ]

    for metric, title, ylabel in plot_specs:
        if metric not in data.columns:
            continue

        out = output_dir / f"{metric}.png"
        save_lineplot(
            data=data,
            y=metric,
            output_path=out,
            title=f"{args.title_prefix}: {title}",
            ylabel=ylabel,
        )

        pooled_out = output_dir / f"{metric}.pooled.png"
        save_pooled_lineplot(
            data=data,
            y=metric,
            output_path=pooled_out,
            title=f"{args.title_prefix}: {title}, pooled samples",
            ylabel=ylabel,
        )

    if "n_compared" in data.columns:
        save_bar_ncompared(
            data=data,
            output_path=output_dir / "n_compared.png",
            title=f"{args.title_prefix}: number of compared sites",
        )

    report_path = make_markdown_report(
        data=data,
        summary=summary,
        output_dir=output_dir,
        title_prefix=args.title_prefix,
    )

    print(f"All metrics: {all_metrics_path}")
    print(f"Summary:     {summary_path}")
    print(f"Report:      {report_path}")
    print(f"Plots:       {output_dir}")


if __name__ == "__main__":
    main()