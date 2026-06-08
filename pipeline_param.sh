#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Parameterized low-pass sequencing pipeline
# ============================================================
#
# Example:
# bash pipeline_param.sh \
#   --input-dir /data/fastq \
#   --output-dir /data/results/208_S15/cov_0.5x/rep_1 \
#   --sample 208_S15 \
#   --sample-reads 2000000 \
#   --seed 101 \
#   --reference /data/ref/genome.fa \
#   --glimpse-panel /data/panel/panel.vcf.gz \
#   --threads 32
#
# Expected input files:
#   ${INPUT_DIR}/${SAMPLE}_R1_001.fastq.gz
#   ${INPUT_DIR}/${SAMPLE}_R2_001.fastq.gz
#
# ============================================================

usage() {
    cat <<EOF
Usage:
  bash pipeline_param.sh \\
    --input-dir DIR \\
    --output-dir DIR \\
    --sample SAMPLE \\
    --sample-reads N \\
    --seed SEED \\
    --reference REF.fa \\
    --glimpse-panel panel.fixed.vcf.gz \\
    --glimpse-sites-vcf panel.sites.vcf.gz \\
    --glimpse-sites-tsv panel.sites.tsv.gz \\
    --glimpse-chunks chunks.txt \\
    [--threads N] \\
    [--skip-varscan] \\
    [--skip-glimpse] \\
    [--skip-qualimap]

Required:
  --input-dir       Directory with FASTQ.gz files
  --output-dir      Output directory for this exact run
  --sample          Sample name, e.g. 208_S15
  --sample-reads    Number of reads per FASTQ file after downsampling
  --seed            Downsampling seed
  --reference       Reference genome FASTA

Optional:
  --glimpse-panel   Reference panel for GLIMPSE2
  --threads         Number of threads. Default: nproc
  --skip-varscan    Do not run VarScan
  --skip-glimpse    Do not run GLIMPSE2 imputation
  --skip-qualimap   Do not run Qualimap
EOF
}

INPUT_DIR=""
OUTPUT_DIR=""
SAMPLE=""
SAMPLE_READS=""
SAMPLE_SEED=""
REF_PATH=""
GLIMPSE2_REF_PANEL=""
THREADS="$(nproc)"

RUN_VARSCAN=1
RUN_GLIMPSE=1
RUN_QUALIMAP=1

BAM_QUALITY=20

# This flag means to exclude reads that^
# * are unmapped
# * are secondary or supplementary alignments
# * failed quality control
FILTER_FLAGS=2820

GLIMPSE2_REF_PANEL=""
GLIMPSE2_SITES_VCF=""
GLIMPSE2_SITES_TSV=""
GLIMPSE2_CHUNKS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --sample)
            SAMPLE="$2"
            shift 2
            ;;
        --sample-reads)
            SAMPLE_READS="$2"
            shift 2
            ;;
        --seed)
            SAMPLE_SEED="$2"
            shift 2
            ;;
        --reference)
            REF_PATH="$2"
            shift 2
            ;;
        --glimpse-panel)
            GLIMPSE2_REF_PANEL="$2"
            shift 2
            ;;
        --glimpse-sites-vcf)
            GLIMPSE2_SITES_VCF="$2"
            shift 2
            ;;
        --glimpse-sites-tsv)
            GLIMPSE2_SITES_TSV="$2"
            shift 2
            ;;
        --glimpse-chunks)
            GLIMPSE2_CHUNKS="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --skip-varscan)
            RUN_VARSCAN=0
            shift
            ;;
        --skip-glimpse)
            RUN_GLIMPSE=0
            shift
            ;;
        --skip-qualimap)
            RUN_QUALIMAP=0
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$SAMPLE" || -z "$SAMPLE_READS" || -z "$SAMPLE_SEED" || -z "$REF_PATH" ]]; then
    echo "ERROR: missing required arguments"
    usage
    exit 1
fi

INPUT_DIR="$(realpath "$INPUT_DIR")"
OUTPUT_DIR="$(realpath -m "$OUTPUT_DIR")"
REF_PATH="$(realpath "$REF_PATH")"

if [[ -n "$GLIMPSE2_REF_PANEL" ]]; then
    GLIMPSE2_REF_PANEL="$(realpath "$GLIMPSE2_REF_PANEL")"
fi

if [[ -n "$GLIMPSE2_SITES_VCF" ]]; then
    GLIMPSE2_SITES_VCF="$(realpath "$GLIMPSE2_SITES_VCF")"
fi

if [[ -n "$GLIMPSE2_SITES_TSV" ]]; then
    GLIMPSE2_SITES_TSV="$(realpath "$GLIMPSE2_SITES_TSV")"
fi

if [[ -n "$GLIMPSE2_CHUNKS" ]]; then
    GLIMPSE2_CHUNKS="$(realpath "$GLIMPSE2_CHUNKS")"
fi

BB_THREADS="$THREADS"

# ============================================================
# Tool paths
# ============================================================

TOOLS="${CONDA_PREFIX:-}/bin"

FASTQC="${TOOLS}/fastqc"
BBDUK="${TOOLS}/bbduk.sh"
REPAIR="${TOOLS}/repair.sh"
REFORMAT="${TOOLS}/reformat.sh"
BWA="${TOOLS}/bwa"
SAMTOOLS="${TOOLS}/samtools"
BCFTOOLS="${TOOLS}/bcftools"
TABIX="${TOOLS}/tabix"
BGZIP="${TOOLS}/bgzip"
GLIMPSE2_CHUNK="${TOOLS}/GLIMPSE2_chunk"
GLIMPSE2_PHASE="${TOOLS}/GLIMPSE2_phase"
GLIMPSE2_LIGATE="${TOOLS}/GLIMPSE2_ligate"
VARSCAN="${TOOLS}/varscan"
VCFTOOLS="${TOOLS}/vcftools"
QUALIMAP="${TOOLS}/qualimap"

BBDUK_REF_DIR="${CONDA_PREFIX}/opt/bbmap-39.26-0/resources"
export BCFTOOLS_PLUGINS="${CONDA_PREFIX}/libexec/bcftools"

# ============================================================
# Input FASTQ
# ============================================================

R1="${INPUT_DIR}/${SAMPLE}_R1_001.fastq.gz"
R2="${INPUT_DIR}/${SAMPLE}_R2_001.fastq.gz"

if [[ ! -f "$R1" ]]; then
    echo "ERROR: R1 not found: $R1"
    exit 1
fi

if [[ ! -f "$R2" ]]; then
    echo "ERROR: R2 not found: $R2"
    exit 1
fi

# ============================================================
# Directories
# ============================================================

SAMPLED_DIR="$OUTPUT_DIR/2_sampled"
ALIGNED_DIR="$OUTPUT_DIR/3_aligned"
MARKDUP_DIR="$OUTPUT_DIR/4_markdup"
FILTERED_DIR="$OUTPUT_DIR/5_filtered"
COVERAGE_DIR="$OUTPUT_DIR/6_coverage"
VCF_DIR="$OUTPUT_DIR/7_vcf"
IMPUTED_DIR="$OUTPUT_DIR/8_imputed"
METRICS_DIR="$OUTPUT_DIR/9_metrics"
RUN_META_DIR="$OUTPUT_DIR/00_run_metadata"

mkdir -p "$SAMPLED_DIR" "$ALIGNED_DIR" \
         "$MARKDUP_DIR" "$FILTERED_DIR" "$COVERAGE_DIR" \
         "$VCF_DIR" "$IMPUTED_DIR" "$METRICS_DIR" "$RUN_META_DIR"

cat > "$RUN_META_DIR/run_config.tsv" <<EOF
sample	$SAMPLE
sample_reads	$SAMPLE_READS
seed	$SAMPLE_SEED
input_dir	$INPUT_DIR
output_dir	$OUTPUT_DIR
reference	$REF_PATH
glimpse_panel	$GLIMPSE2_REF_PANEL
threads	$THREADS
run_varscan	$RUN_VARSCAN
run_glimpse	$RUN_GLIMPSE
run_qualimap	$RUN_QUALIMAP
date	$(date -Is)
EOF

echo "============================================================"
echo "Pipeline run"
echo "Sample:       $SAMPLE"
echo "Sample reads: $SAMPLE_READS"
echo "Seed:         $SAMPLE_SEED"
echo "Output:       $OUTPUT_DIR"
echo "Threads:      $THREADS"
echo "============================================================"


# ============================================================
# STAGE 3 - Downsampling
# ============================================================

echo "STAGE 3 - Downsampling"

SAMPLED_R1="$SAMPLED_DIR/${SAMPLE}_R1_sampled.fastq.gz"
SAMPLED_R2="$SAMPLED_DIR/${SAMPLE}_R2_sampled.fastq.gz"

if [[ -f "$SAMPLED_R1" && -f "$SAMPLED_R2" ]]; then
    echo "SKIP: downsampling already done for $SAMPLE"
else
    echo "INFO: Downsampling $SAMPLE to ${SAMPLE_READS} reads/file, seed=${SAMPLE_SEED}"

    "$REFORMAT" \
        in1="$R1" in2="$R2" \
        out1="$SAMPLED_R1" out2="$SAMPLED_R2" \
        samplereadstarget="$SAMPLE_READS" \
        sampleseed="$SAMPLE_SEED" \
        t="$BB_THREADS"
fi

# ============================================================
# STAGE 4 - Alignment
# ============================================================

echo "STAGE 4 - Alignment"

ALIGNED_BAM="$ALIGNED_DIR/${SAMPLE}.sorted.bam"
ALIGNED_BAI="$ALIGNED_DIR/${SAMPLE}.sorted.bam.bai"

if [[ -f "$ALIGNED_BAM" && -f "$ALIGNED_BAI" ]]; then
    echo "SKIP: alignment already done"
else
    if [[ ! -f "${REF_PATH}.bwt" ]]; then
        echo "INFO: BWA index not found, indexing reference"
        "$BWA" index "$REF_PATH"
    fi

    "$BWA" mem \
        -t "$THREADS" \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
        "$REF_PATH" \
        "$SAMPLED_R1" "$SAMPLED_R2" \
    | "$SAMTOOLS" sort \
        -@ "$THREADS" \
        -o "$ALIGNED_BAM" \
        -T "$ALIGNED_DIR/${SAMPLE}_tmp"

    "$SAMTOOLS" index -@ "$THREADS" "$ALIGNED_BAM"
    "$SAMTOOLS" flagstat "$ALIGNED_BAM" > "$ALIGNED_DIR/${SAMPLE}.flagstat.txt"
fi

# ============================================================
# STAGE 5 - Mark duplicates
# ============================================================

echo "STAGE 5 - Mark duplicates"

MARKDUP_BAM="$MARKDUP_DIR/${SAMPLE}.markdup.bam"
MARKDUP_BAI="$MARKDUP_DIR/${SAMPLE}.markdup.bam.bai"
MARKDUP_STATS="$MARKDUP_DIR/${SAMPLE}.markdup.stats.txt"

if [[ -f "$MARKDUP_BAM" && -f "$MARKDUP_BAI" ]]; then
    echo "SKIP: markdup already done"
else
    "$SAMTOOLS" sort -n \
        -@ "$THREADS" \
        -o "$MARKDUP_DIR/${SAMPLE}.namesort.bam" \
        "$ALIGNED_BAM"

    "$SAMTOOLS" fixmate -m \
        -@ "$THREADS" \
        "$MARKDUP_DIR/${SAMPLE}.namesort.bam" \
        "$MARKDUP_DIR/${SAMPLE}.fixmate.bam"

    "$SAMTOOLS" sort \
        -@ "$THREADS" \
        -o "$MARKDUP_DIR/${SAMPLE}.fixmate.sorted.bam" \
        "$MARKDUP_DIR/${SAMPLE}.fixmate.bam"

    "$SAMTOOLS" markdup \
        -@ "$THREADS" \
        -S \
        -f "$MARKDUP_STATS" \
        "$MARKDUP_DIR/${SAMPLE}.fixmate.sorted.bam" \
        "$MARKDUP_BAM"

    "$SAMTOOLS" index -@ "$THREADS" "$MARKDUP_BAM"

    rm -f "$MARKDUP_DIR/${SAMPLE}.namesort.bam" \
          "$MARKDUP_DIR/${SAMPLE}.fixmate.bam" \
          "$MARKDUP_DIR/${SAMPLE}.fixmate.sorted.bam"
fi

# ============================================================
# STAGE 6 - Filtering
# ============================================================

echo "STAGE 6 - Filtering"

FILTERED_BAM="$FILTERED_DIR/${SAMPLE}.filtered.bam"
FILTERED_BAI="$FILTERED_DIR/${SAMPLE}.filtered.bam.bai"

if [[ -f "$FILTERED_BAM" && -f "$FILTERED_BAI" ]]; then
    echo "SKIP: filtering already done"
else
    "$SAMTOOLS" view \
        -@ "$THREADS" \
        -b \
        -q "$BAM_QUALITY" \
        -F "$FILTER_FLAGS" \
        -o "$FILTERED_BAM" \
        "$MARKDUP_BAM"

    "$SAMTOOLS" index -@ "$THREADS" "$FILTERED_BAM"
    "$SAMTOOLS" flagstat "$FILTERED_BAM" > "$FILTERED_DIR/${SAMPLE}.filtered.flagstat.txt"
fi

# ============================================================
# STAGE 7 - Coverage
# ============================================================

echo "STAGE 7 - Coverage"

DEPTH_SUMMARY="$COVERAGE_DIR/${SAMPLE}.depth.summary.txt"
LOW_COV_REGIONS="$COVERAGE_DIR/${SAMPLE}.low_cov_regions.bed"
DEPTH_DONE_FLAG="$COVERAGE_DIR/${SAMPLE}.depth.done"
DEPTH_RAW="$COVERAGE_DIR/${SAMPLE}.depth.raw.gz"

if [[ -f "$DEPTH_DONE_FLAG" ]]; then
    echo "SKIP: coverage already done"
else
    "$SAMTOOLS" depth -@ "$THREADS" -a "$FILTERED_BAM" | gzip -c > "$DEPTH_RAW"

    zcat "$DEPTH_RAW" | awk '
        BEGIN { total=0; covered=0; sum=0 }
        {
            total++
            sum += $3
            if ($3 >= 1) covered++
        }
        END {
            printf "Total positions:\t%d\n", total
            printf "Mean depth:\t%.6f\n", (total>0 ? sum/total : 0)
            printf "Positions >= 1x:\t%d (%.2f%%)\n", covered, (total>0 ? covered/total*100 : 0)
        }
    ' > "$DEPTH_SUMMARY"

    zcat "$DEPTH_RAW" | awk '$3 < 1 { print $1"\t"($2-1)"\t"$2 }' > "$LOW_COV_REGIONS"

    rm -f "$DEPTH_RAW"
    touch "$DEPTH_DONE_FLAG"

    cat "$DEPTH_SUMMARY"
fi

# ============================================================
# STAGE 8 - SNP calling
# ============================================================

echo "STAGE 8 - SNP calling"

RAW_VCF="$VCF_DIR/${SAMPLE}.raw.vcf.gz"
BCFTOOLS_VCF="$VCF_DIR/${SAMPLE}.bcftools.vcf.gz"
VARSCAN_VCF="$VCF_DIR/${SAMPLE}.varscan.vcf"

echo "--- 8a. BCFtools mpileup + call ---"

if [[ -f "$BCFTOOLS_VCF" && -f "${BCFTOOLS_VCF}.tbi" ]]; then
    echo "SKIP: BCFtools calling already done"
else
    rm -f "$BCFTOOLS_VCF" "${BCFTOOLS_VCF}.tbi"

    "$BCFTOOLS" mpileup \
        --threads "$THREADS" \
        --fasta-ref "$REF_PATH" \
        --min-MQ 20 \
        --min-BQ 20 \
        --annotate FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
        --output-type u \
        "$FILTERED_BAM" \
    | "$BCFTOOLS" call \
        --threads "$THREADS" \
        --multiallelic-caller \
        --variants-only \
        --output-type z \
        --output "$BCFTOOLS_VCF"

    "$BCFTOOLS" index -t --threads "$THREADS" "$BCFTOOLS_VCF" \
        || "$TABIX" -p vcf "$BCFTOOLS_VCF"
fi

echo "--- 8b. VarScan ---"

if [[ "$RUN_VARSCAN" -eq 0 ]]; then
    echo "SKIP: VarScan disabled"
elif [[ ! -x "$VARSCAN" ]]; then
    echo "WARN: VarScan not found, skipping"
elif [[ -f "$VARSCAN_VCF" ]]; then
    echo "SKIP: VarScan already done"
else
    "$SAMTOOLS" mpileup \
        -f "$REF_PATH" \
        -q 20 \
        -Q 20 \
        "$FILTERED_BAM" \
    | "$VARSCAN" mpileup2snp \
        --min-coverage 1 \
        --min-avg-qual 20 \
        --p-value 0.05 \
        --output-vcf 1 \
        > "$VARSCAN_VCF"
fi

# ============================================================
# STAGE 9 - GLIMPSE2 imputation
# ============================================================

echo "STAGE 9 - GLIMPSE2 imputation"

IMPUTED_DIR_CHUNKS="$IMPUTED_DIR/chunks"
IMPUTED_DIR_PHASED="$IMPUTED_DIR/phased"
mkdir -p "$IMPUTED_DIR_CHUNKS" "$IMPUTED_DIR_PHASED"

IMPUTED_VCF="$IMPUTED_DIR/${SAMPLE}.imputed.vcf.gz"

if [[ "$RUN_GLIMPSE" -eq 0 ]]; then
    echo "SKIP: GLIMPSE2 disabled"

elif [[ -z "$GLIMPSE2_REF_PANEL" ]]; then
    echo "WARN: GLIMPSE2 fixed panel not provided, skipping imputation"

elif [[ -z "$GLIMPSE2_SITES_VCF" || -z "$GLIMPSE2_SITES_TSV" || -z "$GLIMPSE2_CHUNKS" ]]; then
    echo "ERROR: GLIMPSE2 prepared files are required:"
    echo "  --glimpse-panel      fixed panel VCF.gz"
    echo "  --glimpse-sites-vcf  sites VCF.gz"
    echo "  --glimpse-sites-tsv  sites TSV.gz"
    echo "  --glimpse-chunks     chunks.txt"
    exit 1

else
    if [[ ! -f "$GLIMPSE2_REF_PANEL" ]]; then
        echo "ERROR: fixed GLIMPSE2 panel not found: $GLIMPSE2_REF_PANEL"
        exit 1
    fi

    if [[ ! -f "${GLIMPSE2_REF_PANEL}.tbi" && ! -f "${GLIMPSE2_REF_PANEL}.csi" ]]; then
        echo "ERROR: fixed GLIMPSE2 panel index not found: ${GLIMPSE2_REF_PANEL}.tbi/.csi"
        exit 1
    fi

    if [[ ! -f "$GLIMPSE2_SITES_VCF" ]]; then
        echo "ERROR: GLIMPSE2 sites VCF not found: $GLIMPSE2_SITES_VCF"
        exit 1
    fi

    if [[ ! -f "${GLIMPSE2_SITES_VCF}.tbi" && ! -f "${GLIMPSE2_SITES_VCF}.csi" ]]; then
        echo "ERROR: GLIMPSE2 sites VCF index not found"
        exit 1
    fi

    if [[ ! -f "$GLIMPSE2_SITES_TSV" ]]; then
        echo "ERROR: GLIMPSE2 sites TSV not found: $GLIMPSE2_SITES_TSV"
        exit 1
    fi

    if [[ ! -f "${GLIMPSE2_SITES_TSV}.tbi" ]]; then
        echo "ERROR: GLIMPSE2 sites TSV tabix index not found: ${GLIMPSE2_SITES_TSV}.tbi"
        exit 1
    fi

    if [[ ! -s "$GLIMPSE2_CHUNKS" ]]; then
        echo "ERROR: GLIMPSE2 chunks file is empty or missing: $GLIMPSE2_CHUNKS"
        exit 1
    fi

    if [[ -f "$IMPUTED_VCF" && -f "${IMPUTED_VCF}.tbi" ]]; then
        echo "SKIP: imputation already done for $SAMPLE"
    else
        # ------------------------------------------------------------
        # 9a. Generate GL VCF for current sample/run
        # ------------------------------------------------------------
        echo "--- 9a. Generate GL VCF by prepared panel sites ---"

        SAMPLE_GL_VCF="$IMPUTED_DIR/${SAMPLE}.gl.vcf.gz"

        if [[ -f "$SAMPLE_GL_VCF" && -f "${SAMPLE_GL_VCF}.tbi" ]]; then
            echo "SKIP: GL VCF already exists: $SAMPLE_GL_VCF"
        else
            "$BCFTOOLS" mpileup \
                --threads "$THREADS" \
                --fasta-ref "$REF_PATH" \
                --min-MQ 20 \
                --min-BQ 20 \
                -I \
                -E \
                -a 'FORMAT/DP,FORMAT/AD' \
                -T "$GLIMPSE2_SITES_VCF" \
                -Ou \
                "$FILTERED_BAM" \
            | "$BCFTOOLS" call \
                --threads "$THREADS" \
                -Aim \
                -C alleles \
                -T "$GLIMPSE2_SITES_TSV" \
                -Oz -o "$SAMPLE_GL_VCF"

            EXIT_CODE=$?

            if [[ $EXIT_CODE -ne 0 || ! -f "$SAMPLE_GL_VCF" ]]; then
                echo "ERROR: bcftools mpileup|call failed with code $EXIT_CODE"
                exit 1
            fi

            "$BCFTOOLS" index -t --threads "$THREADS" "$SAMPLE_GL_VCF"

            N_GL=$("$BCFTOOLS" view -H "$SAMPLE_GL_VCF" | wc -l)
            echo "INFO: GL VCF sites: $N_GL"

            if [[ "$N_GL" -eq 0 ]]; then
                echo "ERROR: GL VCF is empty"
                exit 1
            fi
        fi

        # ------------------------------------------------------------
        # 9b. Determine adaptive GLIMPSE2 parameters
        # ------------------------------------------------------------
        echo "--- 9b. Determine GLIMPSE2 parameters ---"

        N_REF_HAPS=$(( $("$BCFTOOLS" query -l "$GLIMPSE2_REF_PANEL" | wc -l) * 2 ))

        if (( N_REF_HAPS <= 2000 )); then
            KPBWT=$(( N_REF_HAPS / 2 ))
            KINIT=$(( N_REF_HAPS / 2 ))

            if (( KPBWT < 100 )); then
                KPBWT=100
            fi

            if (( KINIT < 100 )); then
                KINIT=100
            fi

            echo "INFO: panel has $N_REF_HAPS haplotypes; using --Kinit $KINIT --Kpbwt $KPBWT"
        else
            KPBWT=2000
            KINIT=1000
        fi

        # ------------------------------------------------------------
        # 9c. GLIMPSE2_phase by prepared chunks
        # ------------------------------------------------------------
        echo "--- 9c. GLIMPSE2_phase by prepared chunks ---"

        PHASE_DONE_FLAG="$IMPUTED_DIR_PHASED/.phase.done"

        if [[ -f "$PHASE_DONE_FLAG" ]]; then
            echo "SKIP: phase chunks already done"
        else
            while IFS=$'\t' read -r ID CHR IRG ORG REST; do
                [[ -z "$ID" || "$ID" == \#* ]] && continue

                PHASED_CHUNK="$IMPUTED_DIR_PHASED/${SAMPLE}.${CHR}.chunk${ID}.bcf"

                if [[ -f "$PHASED_CHUNK" && -f "${PHASED_CHUNK}.csi" ]]; then
                    echo "SKIP: chunk $ID ($IRG) already phased"
                    continue
                fi

                echo "INFO: Phasing chunk ID=$ID CHR=$CHR IRG=$IRG ORG=$ORG"

                "$GLIMPSE2_PHASE" \
                    --input-gl "$SAMPLE_GL_VCF" \
                    --reference "$GLIMPSE2_REF_PANEL" \
                    --input-region "$IRG" \
                    --output-region "$ORG" \
                    --output "$PHASED_CHUNK" \
                    --Kinit "$KINIT" \
                    --Kpbwt "$KPBWT" \
                    --threads "$THREADS"

                EXIT_CODE=$?

                if [[ $EXIT_CODE -ne 0 || ! -f "$PHASED_CHUNK" ]]; then
                    echo "ERROR: GLIMPSE2_phase failed with code $EXIT_CODE"
                    echo "ERROR: ID=$ID IRG=$IRG ORG=$ORG"
                    exit 1
                fi

                "$BCFTOOLS" index -f "$PHASED_CHUNK"

            done < "$GLIMPSE2_CHUNKS"

            touch "$PHASE_DONE_FLAG"
        fi

        # ------------------------------------------------------------
        # 9d. GLIMPSE2_ligate
        # ------------------------------------------------------------
        echo "--- 9d. GLIMPSE2_ligate ---"

        PHASED_LIST="$IMPUTED_DIR_PHASED/phased_chunks.list"

        while IFS=$'\t' read -r ID CHR _; do
            [[ -z "$ID" || "$ID" == \#* ]] && continue

            f="$IMPUTED_DIR_PHASED/${SAMPLE}.${CHR}.chunk${ID}.bcf"

            if [[ -f "$f" ]]; then
                echo "$f"
            else
                echo "ERROR: missing phased chunk: $f"
                exit 1
            fi
        done < "$GLIMPSE2_CHUNKS" > "$PHASED_LIST"

        if [[ ! -s "$PHASED_LIST" ]]; then
            echo "ERROR: phased chunks list is empty"
            exit 1
        fi

        "$GLIMPSE2_LIGATE" \
            --input "$PHASED_LIST" \
            --output "$IMPUTED_VCF" \
            --threads 4

        if [[ $? -ne 0 || ! -f "$IMPUTED_VCF" ]]; then
            echo "ERROR: GLIMPSE2_ligate did not create $IMPUTED_VCF"
            exit 1
        fi

        "$BCFTOOLS" index -t --threads "$THREADS" "$IMPUTED_VCF"

        echo "INFO: Imputed VCF: $IMPUTED_VCF"
        "$BCFTOOLS" stats "$IMPUTED_VCF" | grep "^SN"
    fi
fi

# ============================================================
# STAGE 10 - Basic pipeline QC metrics
# ============================================================

echo "STAGE 10 - Basic pipeline QC metrics"

count_vcf_records() {
    local vcf="$1"
    [[ -f "$vcf" ]] || { echo 0; return; }
    "$BCFTOOLS" view -H "$vcf" 2>/dev/null | wc -l
}

if [[ -n "${IMPUTED_VCF:-}" && -f "$IMPUTED_VCF" && "$(count_vcf_records "$IMPUTED_VCF")" -gt 0 ]]; then
    METRICS_VCF="$IMPUTED_VCF"
elif [[ -f "$BCFTOOLS_VCF" && "$(count_vcf_records "$BCFTOOLS_VCF")" -gt 0 ]]; then
    METRICS_VCF="$BCFTOOLS_VCF"
else
    METRICS_VCF=""
fi

BCFTOOLS_STATS="$METRICS_DIR/${SAMPLE}.bcftools.stats.txt"

if [[ -n "$METRICS_VCF" ]]; then
    "$BCFTOOLS" stats "$METRICS_VCF" > "$BCFTOOLS_STATS"
fi

VCFTOOLS_PREFIX="$METRICS_DIR/${SAMPLE}.vcftools"

if [[ -n "$METRICS_VCF" ]]; then
    "$VCFTOOLS" --gzvcf "$METRICS_VCF" --TsTv-summary --out "$VCFTOOLS_PREFIX" 2> "$VCFTOOLS_PREFIX.TsTv.log" || true
    "$VCFTOOLS" --gzvcf "$METRICS_VCF" --het          --out "$VCFTOOLS_PREFIX" 2> "$VCFTOOLS_PREFIX.het.log" || true
    "$VCFTOOLS" --gzvcf "$METRICS_VCF" --site-depth   --out "$VCFTOOLS_PREFIX" 2> "$VCFTOOLS_PREFIX.depth.log" || true
fi

if [[ "$RUN_QUALIMAP" -eq 1 ]]; then
    QUALIMAP_DIR="$METRICS_DIR/${SAMPLE}.qualimap"
    QUALIMAP_BIN="$(command -v qualimap || true)"
    if [[ -n "$QUALIMAP_BIN" && -f "$FILTERED_BAM" ]]; then
        unset DISPLAY
        "$QUALIMAP_BIN" bamqc \
            -bam "$FILTERED_BAM" \
            -outdir "$QUALIMAP_DIR" \
            -outformat HTML \
            -nt "$THREADS" \
            --java-mem-size=4G \
            || echo "WARN: qualimap failed, continuing"
    fi
fi

SUMMARY_REPORT="$METRICS_DIR/${SAMPLE}.pipeline_summary.txt"

{
    echo "===== Pipeline Summary ====="
    echo "Sample: $SAMPLE"
    echo "Sample reads: $SAMPLE_READS"
    echo "Seed: $SAMPLE_SEED"
    echo "Output dir: $OUTPUT_DIR"
    echo "Date: $(date -Is)"
    echo ""

    echo "--- Coverage summary ---"
    if [[ -f "$DEPTH_SUMMARY" ]]; then
        cat "$DEPTH_SUMMARY"
    else
        echo "No coverage summary"
    fi
    echo ""

    echo "--- VCF used for basic stats ---"
    if [[ -n "$METRICS_VCF" ]]; then
        echo "$METRICS_VCF"
        echo "Records: $(count_vcf_records "$METRICS_VCF")"
    else
        echo "No VCF"
    fi
} > "$SUMMARY_REPORT"

echo "============================================================"
echo "Pipeline finished"
echo "Sample:  $SAMPLE"
echo "Output:  $OUTPUT_DIR"
echo "Summary: $SUMMARY_REPORT"
echo "============================================================"