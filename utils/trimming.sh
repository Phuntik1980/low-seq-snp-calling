#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
THREADS="$(nproc)"
JOBS=1
SAMPLES=()
SAMPLES_FILE=""

usage() {
    cat <<EOF
Usage:
  bash trimming.sh \\
    --input-dir DIR \\
    --output-dir DIR \\
    --sample SAMPLE [--sample SAMPLE2 ...] \\
    [--samples-file FILE] \\
    [--threads N] \\
    [--jobs N]

Required:
  --input-dir        Directory with FASTQ.gz files
  --output-dir       Output directory for this exact run

Samples:
  --sample           Sample name, e.g. 208_S15
                     Can be used multiple times
  --samples-file     Text file with one sample name per line

Optional:
  --threads          Number of threads per sample. Default: nproc
  --jobs             Number of samples to process simultaneously. Default: 1

Example:
  bash trimming.sh \\
    --input-dir ./fastq \\
    --output-dir ./results \\
    --sample 208_S15 \\
    --sample 209_S16 \\
    --threads 4 \\
    --jobs 2

Example with samples file:
  bash trimming.sh \\
    --input-dir ./fastq \\
    --output-dir ./results \\
    --samples-file samples.txt \\
    --threads 4 \\
    --jobs 3
EOF
}

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
            SAMPLES+=("$2")
            shift 2
            ;;
        --samples-file)
            SAMPLES_FILE="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --jobs)
            JOBS="$2"
            shift 2
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

# -------------------------------
# Checks
# -------------------------------

if [[ -z "$INPUT_DIR" ]]; then
    echo "ERROR: --input-dir is required"
    usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: --output-dir is required"
    usage
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: input directory does not exist: $INPUT_DIR"
    exit 1
fi

if [[ -n "$SAMPLES_FILE" ]]; then
    if [[ ! -f "$SAMPLES_FILE" ]]; then
        echo "ERROR: samples file does not exist: $SAMPLES_FILE"
        exit 1
    fi

    while IFS= read -r sample; do
        [[ -z "$sample" ]] && continue
        [[ "$sample" =~ ^# ]] && continue
        SAMPLES+=("$sample")
    done < "$SAMPLES_FILE"
fi

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    echo "ERROR: no samples specified. Use --sample or --samples-file"
    usage
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

TRIMMED_DIR="$OUTPUT_DIR/1_trimmed"
mkdir -p "$TRIMMED_DIR"

TOOLS="${CONDA_PREFIX:-}/bin"

BBDUK="${TOOLS}/bbduk.sh"
REPAIR="${TOOLS}/repair.sh"

if [[ ! -x "$BBDUK" ]]; then
    echo "ERROR: bbduk.sh not found or not executable: $BBDUK"
    exit 1
fi

if [[ ! -x "$REPAIR" ]]; then
    echo "ERROR: repair.sh not found or not executable: $REPAIR"
    exit 1
fi

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: CONDA_PREFIX is not set. Activate conda environment with bbmap."
    exit 1
fi

BBDUK_REF_DIR="${CONDA_PREFIX}/opt/bbmap-39.26-0/resources"

if [[ ! -d "$BBDUK_REF_DIR" ]]; then
    echo "ERROR: BBMap resources directory not found: $BBDUK_REF_DIR"
    exit 1
fi

BB_THREADS="$THREADS"

# ============================================================
# Function for one sample
# ============================================================

trim_sample() {
    local SAMPLE="$1"

    echo "=============================================="
    echo "Trimming sample: $SAMPLE"
    echo "=============================================="

    local R1=""
    local R2=""

    if compgen -G "$INPUT_DIR/${SAMPLE}_R1_001.fastq.gz" > /dev/null; then
        R1="$(ls "$INPUT_DIR/${SAMPLE}"_R1_001.fastq.gz | head -n 1)"
    fi

    if compgen -G "$INPUT_DIR/${SAMPLE}_R2_001.fastq.gz" > /dev/null; then
        R2="$(ls "$INPUT_DIR/${SAMPLE}"_R2_001.fastq.gz | head -n 1)"
    fi

    if [[ -z "$R1" || -z "$R2" ]]; then
        echo "ERROR: cannot find R1/R2 FASTQ files for sample: $SAMPLE"
        echo "Expected pattern:"
        echo "  $INPUT_DIR/${SAMPLE}_R1_001.fastq.gz"
        echo "  $INPUT_DIR/${SAMPLE}_R2_001.fastq.gz"
        return 1
    fi

    echo "Input R1: $R1"
    echo "Input R2: $R2"

    local TRIMMED_R1="$TRIMMED_DIR/${SAMPLE}_R1_trimmed.fastq.gz"
    local TRIMMED_R2="$TRIMMED_DIR/${SAMPLE}_R2_trimmed.fastq.gz"

    local TRIMMED_FR1="$TRIMMED_DIR/${SAMPLE}_R1_trimmed_fixed.fastq.gz"
    local TRIMMED_FR2="$TRIMMED_DIR/${SAMPLE}_R2_trimmed_fixed.fastq.gz"

    local SINGLETONS="$TRIMMED_DIR/${SAMPLE}_singletons.fastq.gz"
    local LOG="$TRIMMED_DIR/${SAMPLE}_bbduk.log"

    if [[ -f "$TRIMMED_R1" && -f "$TRIMMED_R2" ]]; then
        echo "SKIP: trimming already done for $SAMPLE"
        return 0
    fi

    "$BBDUK" \
        in1="$R1" \
        in2="$R2" \
        out1="$TRIMMED_R1" \
        out2="$TRIMMED_R2" \
        ref="$BBDUK_REF_DIR/sequencing_artifacts.fa.gz,$BBDUK_REF_DIR/phix174_ill.ref.fa.gz,$BBDUK_REF_DIR/adapters.fa" \
        k=31 \
        ordered \
        cardinality \
        qtrim=rl \
        trimq=20 \
        maq=25 \
        tbo \
        mink=11 \
        ktrim=r \
        minlen=50 \
        t="$BB_THREADS" \
        2> "$LOG"

    "$REPAIR" \
        in1="$TRIMMED_R1" \
        in2="$TRIMMED_R2" \
        out1="$TRIMMED_FR1" \
        out2="$TRIMMED_FR2" \
        outsingle="$SINGLETONS" \
        threads="$BB_THREADS"

    mv "$TRIMMED_FR1" "$TRIMMED_R1"
    mv "$TRIMMED_FR2" "$TRIMMED_R2"

    echo "DONE: $SAMPLE"
}

export -f trim_sample

export INPUT_DIR
export TRIMMED_DIR
export BBDUK
export REPAIR
export BBDUK_REF_DIR
export BB_THREADS

# ============================================================
# Run samples in parallel
# ============================================================

echo "Samples to process: ${#SAMPLES[@]}"
echo "Threads per sample: $THREADS"
echo "Parallel samples: $JOBS"

printf "%s\n" "${SAMPLES[@]}" | xargs -n 1 -P "$JOBS" -I {} bash -c 'trim_sample "$@"' _ {}

echo "All trimming jobs finished."