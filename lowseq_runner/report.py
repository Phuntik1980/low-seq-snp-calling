from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def read_metrics(metrics_dir: Path) -> pd.DataFrame:
    files = sorted(metrics_dir.glob("*.metrics.tsv"))

    if not files:
        raise FileNotFoundError(
            f"No *.metrics.tsv files found in {metrics_dir}")

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
        "precision",
        "positive_predictive_value",
        "recall",
        "sensitivity",
        "specificity",
        "f1_score",
        "n_dosage_pairs",
        "dosage_pearson_r",
        "dosage_r2",
        "binary_tp",
        "binary_fp",
        "binary_tn",
        "binary_fn",
    ]

    for col in numeric_cols:
        if col in data.columns:
            data[col] = pd.to_numeric(data[col], errors="coerce")

    if "coverage" in data.columns:
        data["coverage_label"] = data["coverage"].map(lambda x: f"{x: g}x")

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
        "precision",
        "recall",
        "sensitivity",
        "specificity",
        "f1_score",
        "dosage_pearson_r",
        "dosage_r2",
    ]

    available = [m for m in metrics if m in data.columns]

    summary = data.groupby(["sample", "method", "coverage"],
                           as_index=False).agg(
                               **{f"{m}_mean": (m, "mean")
                                  for m in available},
                               **{f"{m}_sd": (m, "std")
                                  for m in available},
                               n_replicates=("replicate", "nunique"),
                           )

    return summary


def save_lineplot(
        data: pd.DataFrame,
        y: str,
        output_path: Path,
        title: str,
        ylabel: str,
        ylim: tuple[float, float] | None = (0, 1.02),
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

    if ylim is not None:
        plt.ylim(*ylim)

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
        ylim: tuple[float, float] | None = (0, 1.02),
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

    if ylim is not None:
        plt.ylim(*ylim)

    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def save_ncompared_plot(data: pd.DataFrame, output_path: Path, title: str):
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


def _coverage_order(data: pd.DataFrame) -> list[str]:
    covs = sorted(data["coverage"].dropna().unique())
    return [f"{x: g}x" for x in covs]


def save_metric_boxplot(
        data: pd.DataFrame,
        y: str,
        output_path: Path,
        title: str,
        ylabel: str,
        xlabel: str = "Imputed dataset",
        baseline: float | None = None,
        ylim: tuple[float, float] | None = (0, 1.02),
):
    plot_data = data.dropna(subset=[y, "coverage_label"]).copy()

    if plot_data.empty:
        return

    order = _coverage_order(plot_data)

    plt.figure(figsize=(10, 6))

    palette = sns.color_palette("Blues", n_colors=max(len(order), 3))

    ax = sns.boxplot(
        data=plot_data,
        x="coverage_label",
        y=y,
        order=order,
        palette=palette,
        width=0.65,
        showfliers=True,
        linewidth=1.0,
        medianprops={
            "color": "black",
            "linewidth": 1.8
        },
        boxprops={
            "edgecolor": "black",
            "linewidth": 0.8
        },
        whiskerprops={
            "color": "black",
            "linestyle": "--",
            "linewidth": 0.8
        },
        capprops={
            "color": "black",
            "linewidth": 0.8
        },
    )

    sns.stripplot(
        data=plot_data,
        x="coverage_label",
        y=y,
        order=order,
        color="black",
        alpha=0.35,
        size=3,
        jitter=0.18,
        ax=ax,
    )

    if baseline is not None:
        plt.axhline(
            baseline,
            color="gray",
            linestyle=":",
            linewidth=1.2,
            label=f"Baseline = {baseline: .3f}",
        )
        plt.legend(loc="lower right")

    plt.title(title, fontweight="bold")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if ylim is not None:
        plt.ylim(*ylim)

    plt.grid(True, axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def save_sample_faceted_boxplot(
        data: pd.DataFrame,
        y: str,
        output_path: Path,
        title: str,
        ylabel: str,
        ylim: tuple[float, float] | None = (0, 1.02),
):
    plot_data = data.dropna(subset=[y, "coverage_label"]).copy()

    if plot_data.empty:
        return

    order = _coverage_order(plot_data)

    g = sns.catplot(
        data=plot_data,
        x="coverage_label",
        y=y,
        col="sample",
        kind="box",
        order=order,
        col_wrap=3,
        height=4,
        aspect=1.0,
        palette=sns.color_palette("Blues", n_colors=max(len(order), 3)),
        showfliers=True,
        linewidth=1.0,
        medianprops={
            "color": "black",
            "linewidth": 1.6
        },
        whiskerprops={
            "color": "black",
            "linestyle": "--",
            "linewidth": 0.8
        },
    )

    g.set_axis_labels("Target coverage", ylabel)
    g.set_titles("{col_name}")

    if ylim is not None:
        for ax in g.axes.flat:
            ax.set_ylim(*ylim)
            ax.grid(True, axis="y", alpha=0.25)

    g.fig.suptitle(title, fontweight="bold", y=1.04)
    g.tight_layout()
    g.savefig(output_path, dpi=300)
    plt.close(g.fig)


def save_sensitivity_specificity_panel(
    data: pd.DataFrame,
    output_path: Path,
    title_prefix: str,
):
    required = {"sensitivity", "specificity", "coverage_label"}
    if not required.issubset(set(data.columns)):
        return

    plot_data = data.dropna(
        subset=["sensitivity", "specificity", "coverage_label"]).copy()

    if plot_data.empty:
        return

    order = _coverage_order(plot_data)
    palette = sns.color_palette("Blues", n_colors=max(len(order), 3))

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sns.boxplot(
        data=plot_data,
        x="coverage_label",
        y="sensitivity",
        order=order,
        palette=palette,
        ax=axes[0],
        width=0.65,
        medianprops={
            "color": "black",
            "linewidth": 1.8
        },
        whiskerprops={
            "color": "black",
            "linestyle": "--",
            "linewidth": 0.8
        },
    )
    sns.stripplot(
        data=plot_data,
        x="coverage_label",
        y="sensitivity",
        order=order,
        color="black",
        alpha=0.35,
        size=3,
        jitter=0.18,
        ax=axes[0],
    )
    axes[0].set_title("Sensitivity", fontweight="bold")
    axes[0].set_xlabel("Imputed dataset")
    axes[0].set_ylabel("Sensitivity relative to high coverage call set")
    axes[0].set_ylim(0, 1.02)
    axes[0].grid(True, axis="y", alpha=0.25)

    sns.boxplot(
        data=plot_data,
        x="coverage_label",
        y="specificity",
        order=order,
        palette=palette,
        ax=axes[1],
        width=0.65,
        medianprops={
            "color": "black",
            "linewidth": 1.8
        },
        whiskerprops={
            "color": "black",
            "linestyle": "--",
            "linewidth": 0.8
        },
    )
    sns.stripplot(
        data=plot_data,
        x="coverage_label",
        y="specificity",
        order=order,
        color="black",
        alpha=0.35,
        size=3,
        jitter=0.18,
        ax=axes[1],
    )
    axes[1].set_title("Specificity", fontweight="bold")
    axes[1].set_xlabel("Imputed dataset")
    axes[1].set_ylabel("Specificity relative to high coverage call set")
    axes[1].set_ylim(0.95, 1.0005)
    axes[1].grid(True, axis="y", alpha=0.25)

    fig.suptitle(f"{title_prefix}: sensitivity and specificity",
                 fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def make_markdown_report(
    data: pd.DataFrame,
    summary: pd.DataFrame,
    output_dir: Path,
    title_prefix: str,
) -> Path:
    report_path = output_dir / "report.md"

    metrics_to_show = [
        ("exact_match_rate", "Exact genotype concordance"),
        ("alt_presence_match_rate", "Alternative allele presence concordance"),
        ("het_status_match_rate", "Heterozygosity status concordance"),
        ("het_recall", "WGS heterozygote recall"),
        ("hom_alt_recall", "WGS homozygous alternative recall"),
        ("alt_recall", "WGS non-reference recall"),
        ("precision", "Precision / positive predictive value"),
        ("sensitivity", "Sensitivity / recall"),
        ("specificity", "Specificity"),
        ("f1_score", "F1-score"),
        ("dosage_pearson_r", "ALT dosage Pearson correlation"),
        ("dosage_r2", "ALT dosage R²"),
    ]

    boxplots_to_show = [
        ("sensitivity.boxplot.png", "Sensitivity boxplot"),
        ("specificity.boxplot.png", "Specificity boxplot"),
        ("f1_score.boxplot.png", "F1-score boxplot"),
        ("dosage_r2.boxplot.png", "Dosage R² boxplot"),
        ("dosage_pearson_r.boxplot.png", "Dosage correlation boxplot"),
        ("sensitivity_specificity.panel.png",
         "Sensitivity and specificity panel"),
        ("f1_score.by_sample.boxplot.png", "F1-score by sample"),
        ("dosage_r2.by_sample.boxplot.png", "Dosage R² by sample"),
    ]

    with report_path.open("w") as f:
        f.write(f"# {title_prefix}: genotype concordance report\n\n")

        f.write("## Input summary\n\n")
        f.write(f"- Number of metric rows: {len(data)}\n")
        f.write(
            f"- Samples: {', '.join(map(str, sorted(data['sample'].unique())))}"
            f"\n")
        f.write(
            f"- Methods: {', '.join(map(str, sorted(data['method'].unique())))}"
            f"\n")
        _coverages = (', '.join(
            map(lambda x: f'{x: g}',
                sorted(data['coverage'].dropna().unique()))))
        f.write(f"- Coverages: {_coverages}\n\n")

        if "n_compared" in data.columns:
            f.write(f"- Median compared sites: "
                    f"{data['n_compared'].median(): .0f}\n")
        if "n_dosage_pairs" in data.columns:
            f.write(
                f"- Median dosage pairs: {data['n_dosage_pairs'].median(): .0f}"
                f"\n")

        f.write("\n")

        f.write("## Additional metrics\n\n")
        f.write(
            "- **Sensitivity / recall**: TP / (TP + FN), where positive means "
            "truth genotype contains at least one ALT allele.\n")
        f.write("- **Specificity**: TN / (TN + FP), "
                "where negative means truth genotype "
                "is homozygous reference.\n")
        f.write("- **Precision**: TP / (TP + FP).\n")
        f.write("- **F1-score**: harmonic mean of precision and recall.\n")
        f.write("- **Dosage R²**: squared "
                "Pearson correlation between truth ALT dosage "
                "and query ALT dosage. For diploid biallelic "
                "SNPs dosage is 0, 1 or 2.\n\n")

        f.write("## Main line plots\n\n")

        for metric, label in metrics_to_show:
            png = f"{metric}.png"
            pooled_png = f"{metric}.pooled.png"

            if (output_dir / png).exists():
                f.write(f"### {label}\n\n")
                f.write(f"![{label}]({png})\n\n")

            if (output_dir / pooled_png).exists():
                f.write(f"![{label}, pooled]({pooled_png})\n\n")

        if (output_dir / "n_compared.png").exists():
            f.write("### Number of compared sites\n\n")
            f.write("![Number of compared sites](n_compared.png)\n\n")

        f.write("## Boxplot-style accuracy plots\n\n")

        for png, label in boxplots_to_show:
            if (output_dir / png).exists():
                f.write(f"### {label}\n\n")
                f.write(f"![{label}]({png})\n\n")

        f.write("## Summary table\n\n")

        selected_cols = [
            col for col in summary.columns if col in [
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
                "precision_mean",
                "precision_sd",
                "sensitivity_mean",
                "sensitivity_sd",
                "specificity_mean",
                "specificity_sd",
                "f1_score_mean",
                "f1_score_sd",
                "dosage_pearson_r_mean",
                "dosage_pearson_r_sd",
                "dosage_r2_mean",
                "dosage_r2_sd",
            ]
        ]

        f.write(summary[selected_cols].to_markdown(index=False))
        f.write("\n")

    return report_path


def generate_report(
    metrics_dir: Path,
    output_dir: Path,
    title_prefix: str,
) -> Path:
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
            (0, 1.02),
        ),
        (
            "alt_presence_match_rate",
            "Alternative allele presence concordance",
            "Fraction with matching ALT presence",
            (0, 1.02),
        ),
        (
            "het_status_match_rate",
            "Heterozygosity status concordance",
            "Fraction with matching heterozygosity status",
            (0, 1.02),
        ),
        (
            "het_recall",
            "WGS heterozygote recall",
            "Fraction of WGS heterozygotes called as heterozygotes",
            (0, 1.02),
        ),
        (
            "hom_alt_recall",
            "WGS homozygous ALT recall",
            "Fraction of WGS hom-alt sites called as hom-alt",
            (0, 1.02),
        ),
        (
            "alt_recall",
            "WGS non-reference recall",
            "Fraction of WGS non-ref sites called as non-ref",
            (0, 1.02),
        ),
        (
            "nonref_exact_match_rate",
            "Exact concordance on WGS non-reference sites",
            "Exact match among WGS non-ref sites",
            (0, 1.02),
        ),
        (
            "precision",
            "Precision / positive predictive value",
            "TP / (TP + FP)",
            (0, 1.02),
        ),
        (
            "sensitivity",
            "Sensitivity / recall",
            "TP / (TP + FN)",
            (0, 1.02),
        ),
        (
            "specificity",
            "Specificity",
            "TN / (TN + FP)",
            (0.95, 1.0005),
        ),
        (
            "f1_score",
            "F1-score",
            "Harmonic mean of precision and recall",
            (0, 1.02),
        ),
        (
            "dosage_pearson_r",
            "ALT dosage Pearson correlation",
            "Pearson correlation between truth/query ALT dosage",
            (-0.05, 1.02),
        ),
        (
            "dosage_r2",
            "ALT dosage R²",
            "Squared Pearson correlation between truth/query ALT dosage",
            (-0.05, 1.02),
        ),
    ]

    for metric, title, ylabel, ylim in plot_specs:
        if metric not in data.columns:
            continue

        save_lineplot(
            data=data,
            y=metric,
            output_path=output_dir / f"{metric}.png",
            title=f"{title_prefix}: {title}",
            ylabel=ylabel,
            ylim=ylim,
        )

        save_pooled_lineplot(
            data=data,
            y=metric,
            output_path=output_dir / f"{metric}.pooled.png",
            title=f"{title_prefix}: {title}, pooled samples",
            ylabel=ylabel,
            ylim=ylim,
        )

    if "n_compared" in data.columns:
        save_ncompared_plot(
            data=data,
            output_path=output_dir / "n_compared.png",
            title=f"{title_prefix}: number of compared sites",
        )

    # Boxplots similar to examples.
    if "sensitivity" in data.columns:
        save_metric_boxplot(
            data=data,
            y="sensitivity",
            output_path=output_dir / "sensitivity.boxplot.png",
            title=f"{title_prefix}: Sensitivity",
            ylabel="Sensitivity relative to high coverage call set",
            xlabel="Imputed dataset",
            baseline=data["sensitivity"].median(),
            ylim=(0, 1.02),
        )

    if "specificity" in data.columns:
        save_metric_boxplot(
            data=data,
            y="specificity",
            output_path=output_dir / "specificity.boxplot.png",
            title=f"{title_prefix}: Specificity",
            ylabel="Specificity relative to high coverage call set",
            xlabel="Imputed dataset",
            baseline=data["specificity"].median(),
            ylim=(0.95, 1.0005),
        )

    if "f1_score" in data.columns:
        save_metric_boxplot(
            data=data,
            y="f1_score",
            output_path=output_dir / "f1_score.boxplot.png",
            title=f"{title_prefix}: F1-score",
            ylabel="F1-score",
            xlabel="Imputed dataset",
            baseline=data["f1_score"].median(),
            ylim=(0, 1.02),
        )

        save_sample_faceted_boxplot(
            data=data,
            y="f1_score",
            output_path=output_dir / "f1_score.by_sample.boxplot.png",
            title=f"{title_prefix}: F1-score by sample",
            ylabel="F1-score",
            ylim=(0, 1.02),
        )

    if "dosage_r2" in data.columns:
        save_metric_boxplot(
            data=data,
            y="dosage_r2",
            output_path=output_dir / "dosage_r2.boxplot.png",
            title=f"{title_prefix}: ALT dosage R²",
            ylabel="R² relative to high coverage genotypes",
            xlabel="Imputed dataset",
            baseline=data["dosage_r2"].median(),
            ylim=(-0.05, 1.02),
        )

        save_sample_faceted_boxplot(
            data=data,
            y="dosage_r2",
            output_path=output_dir / "dosage_r2.by_sample.boxplot.png",
            title=f"{title_prefix}: ALT dosage R² by sample",
            ylabel="ALT dosage R²",
            ylim=(-0.05, 1.02),
        )

    if "dosage_pearson_r" in data.columns:
        save_metric_boxplot(
            data=data,
            y="dosage_pearson_r",
            output_path=output_dir / "dosage_pearson_r.boxplot.png",
            title=f"{title_prefix}: correlation with high-coverage genotypes",
            ylabel="Correlation with high-coverage genotypes",
            xlabel="Imputed dataset",
            baseline=data["dosage_pearson_r"].median(),
            ylim=(-0.05, 1.02),
        )

    save_sensitivity_specificity_panel(
        data=data,
        output_path=output_dir / "sensitivity_specificity.panel.png",
        title_prefix=title_prefix,
    )

    report_path = make_markdown_report(
        data=data,
        summary=summary,
        output_dir=output_dir,
        title_prefix=title_prefix,
    )

    print(f"All metrics: {all_metrics_path}")
    print(f"Summary:     {summary_path}")
    print(f"Report:      {report_path}")
    print(f"Plots:       {output_dir}")

    return report_path
