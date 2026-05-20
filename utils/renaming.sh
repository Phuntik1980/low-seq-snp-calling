#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage:
  bash prepare_glimpse_panel.sh \\
    --reference REF.fa \\
    --panel panel.vcf.gz \\
    --output-dir DIR \\
    [--threads N] \\
    [--rename-table rename.tsv] \\
    [--chroms chroms.txt] \\
    [--window-mb 2] \\
    [--buffer-mb 0.4]

Required:
  --reference       Reference FASTA
  --panel           Input reference panel VCF.gz
  --output-dir      Output directory for prepared GLIMPSE2 files

Optional:
  --threads         Threads, default: nproc
  --rename-table    TSV OLD_PANEL_CHR<TAB>REFERENCE_CHR
  --chroms          List of chromosomes to chunk, one per line.
                    If omitted, use chromosomes from fixed panel.
  --window-mb       GLIMPSE2 chunk window size, default 2
  --buffer-mb       GLIMPSE2 chunk buffer size, default 0.4
EOF
}

REF_PATH=""
PANEL_IN=""
OUTPUT_DIR=""
THREADS="$(nproc)"
RENAME_TABLE=""
CHROMS_FILE=""
WINDOW_MB="2"
BUFFER_MB="0.4"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference)
            REF_PATH="$2"
            shift 2
            ;;
        --panel)
            PANEL_IN="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --rename-table)
            RENAME_TABLE="$2"
            shift 2
            ;;
        --chroms)
            CHROMS_FILE="$2"
            shift 2
            ;;
        --window-mb)
            WINDOW_MB="$2"
            shift 2
            ;;
        --buffer-mb)
            BUFFER_MB="$2"
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

if [[ -z "$REF_PATH" || -z "$PANEL_IN" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: missing required arguments"
    usage
    exit 1
fi

REF_PATH="$(realpath "$REF_PATH")"
PANEL_IN="$(realpath "$PANEL_IN")"
OUTPUT_DIR="$(realpath -m "$OUTPUT_DIR")"

mkdir -p "$OUTPUT_DIR"

TOOLS="${CONDA_PREFIX:-}/bin"
BCFTOOLS="${TOOLS}/bcftools"
TABIX="${TOOLS}/tabix"
BGZIP="${TOOLS}/bgzip"
GLIMPSE2_CHUNK="${TOOLS}/GLIMPSE2_chunk"

export BCFTOOLS_PLUGINS="${CONDA_PREFIX}/libexec/bcftools"

if [[ ! -f "$REF_PATH" ]]; then
    echo "ERROR: reference not found: $REF_PATH"
    exit 1
fi

if [[ ! -f "$PANEL_IN" ]]; then
    echo "ERROR: panel not found: $PANEL_IN"
    exit 1
fi

if [[ ! -f "${REF_PATH}.fai" ]]; then
    echo "ERROR: reference index not found: ${REF_PATH}.fai"
    echo "Run: samtools faidx $REF_PATH"
    exit 1
fi

if [[ ! -f "${PANEL_IN}.tbi" && ! -f "${PANEL_IN}.csi" ]]; then
    echo "INFO: indexing input panel"
    "$BCFTOOLS" index -t --threads "$THREADS" "$PANEL_IN" || "$TABIX" -p vcf "$PANEL_IN"
fi

echo "============================================================"
echo "Prepare GLIMPSE2 panel"
echo "Reference:  $REF_PATH"
echo "Panel:      $PANEL_IN"
echo "Output dir: $OUTPUT_DIR"
echo "Threads:    $THREADS"
echo "============================================================"

# ------------------------------------------------------------
# 1. Optional chromosome rename
# ------------------------------------------------------------

PANEL_WORK="$PANEL_IN"

if [[ -n "$RENAME_TABLE" ]]; then
    RENAME_TABLE="$(realpath "$RENAME_TABLE")"

    if [[ ! -s "$RENAME_TABLE" ]]; then
        echo "ERROR: rename table is empty or missing: $RENAME_TABLE"
        exit 1
    fi

    PANEL_RENAMED="$OUTPUT_DIR/panel.renamed.vcf.gz"
    PANEL_RENAMED_UNSORTED="$OUTPUT_DIR/panel.renamed.unsorted.vcf.gz"

    if [[ -f "$PANEL_RENAMED" && -f "${PANEL_RENAMED}.tbi" ]]; then
        echo "SKIP: renamed panel exists: $PANEL_RENAMED"
    else
        echo "INFO: renaming chromosomes by $RENAME_TABLE"

        "$BCFTOOLS" annotate \
            --threads "$THREADS" \
            --rename-chrs "$RENAME_TABLE" \
            "$PANEL_IN" \
            -Oz -o "$PANEL_RENAMED_UNSORTED"

        "$BCFTOOLS" sort \
            -Oz -o "$PANEL_RENAMED" \
            -T "$OUTPUT_DIR/bcftools_sort_tmp_rename" \
            "$PANEL_RENAMED_UNSORTED"

        rm -f "$PANEL_RENAMED_UNSORTED"

        "$BCFTOOLS" index -t -f --threads "$THREADS" "$PANEL_RENAMED"
    fi

    PANEL_WORK="$PANEL_RENAMED"
fi

# ------------------------------------------------------------
# 2. Fix panel for GLIMPSE2
# ------------------------------------------------------------

PANEL_FIXED="$OUTPUT_DIR/panel.fixed.vcf.gz"

if [[ -f "$PANEL_FIXED" && -f "${PANEL_FIXED}.tbi" ]]; then
    echo "SKIP: fixed panel exists: $PANEL_FIXED"
else
    N_SAMPLES=$("$BCFTOOLS" query -l "$PANEL_WORK" | wc -l)
    EXPECTED_AN=$((N_SAMPLES * 2))

    echo "INFO: panel samples: $N_SAMPLES"
    echo "INFO: expected AN: $EXPECTED_AN"
    echo "INFO: norm -> biallelic SNP -> fill missing GT -> phase -> fill AC/AN -> drop monomorphic"

    EXP_AN="$EXPECTED_AN"

    "$BCFTOOLS" norm --threads "$THREADS" -f "$REF_PATH" -m -any "$PANEL_WORK" -Ou \
    | "$BCFTOOLS" view --threads "$THREADS" -m 2 -M 2 -v snps -Ou \
    | "$BCFTOOLS" +setGT -Ou -- -t . -n 0 \
    | "$BCFTOOLS" +setGT -Ou -- -t ./x -n 0 \
    | "$BCFTOOLS" +setGT -Ou -- -t a -n p \
    | "$BCFTOOLS" annotate --threads "$THREADS" -x 'INFO,^FORMAT/GT' -Ou \
    | "$BCFTOOLS" +fill-tags -Ou -- -t AC,AN,AF,NS \
    | "$BCFTOOLS" view --threads "$THREADS" \
        -e "AC==0 || AC==AN || INFO/AN!=$EXP_AN" \
        -Oz -o "$PANEL_FIXED"

    "$BCFTOOLS" index -t -f --threads "$THREADS" "$PANEL_FIXED"

    N_BAD=$("$BCFTOOLS" query -f '%INFO/AN\n' "$PANEL_FIXED" \
        | awk -v want="$EXPECTED_AN" '$1 != want' \
        | wc -l)

    if [[ "$N_BAD" -gt 0 ]]; then
        echo "ERROR: $N_BAD sites have AN != $EXPECTED_AN"
        "$BCFTOOLS" query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' "$PANEL_FIXED" \
            | awk -v want="$EXPECTED_AN" '$4 != want' \
            | head -5
        exit 1
    fi

    N_UNPHASED=$("$BCFTOOLS" query -f '[%GT\n]' "$PANEL_FIXED" \
        | head -10000 \
        | grep -c '/' || true)

    if [[ "$N_UNPHASED" -gt 0 ]]; then
        echo "WARN: unphased GT remain in first 10000 GT: $N_UNPHASED"
    fi
fi

# ------------------------------------------------------------
# 3. Create sites VCF and TSV
# ------------------------------------------------------------

SITES_VCF="$OUTPUT_DIR/panel.sites.vcf.gz"
SITES_TSV="$OUTPUT_DIR/panel.sites.tsv.gz"

if [[ -f "$SITES_VCF" && -f "${SITES_VCF}.tbi" ]]; then
    echo "SKIP: sites VCF exists: $SITES_VCF"
else
    echo "INFO: creating sites VCF"

    "$BCFTOOLS" view \
        -G \
        "$PANEL_FIXED" \
        -Oz -o "$SITES_VCF"

    "$BCFTOOLS" index -t -f --threads "$THREADS" "$SITES_VCF"
fi

if [[ -f "$SITES_TSV" && -f "${SITES_TSV}.tbi" ]]; then
    echo "SKIP: sites TSV exists: $SITES_TSV"
else
    echo "INFO: creating sites TSV"

    "$BCFTOOLS" query \
        -f '%CHROM\t%POS\t%REF,%ALT\n' \
        "$SITES_VCF" \
    | "$BGZIP" -c > "$SITES_TSV"

    "$TABIX" -f -s1 -b2 -e2 "$SITES_TSV"
fi

# ------------------------------------------------------------
# 4. Create chunks
# ------------------------------------------------------------

CHUNKS_FILE="$OUTPUT_DIR/chunks.txt"
CHUNKS_DIR="$OUTPUT_DIR/chunks"
mkdir -p "$CHUNKS_DIR"

if [[ -f "$CHUNKS_FILE" && -s "$CHUNKS_FILE" ]]; then
    echo "SKIP: chunks file exists: $CHUNKS_FILE"
else
    echo "INFO: creating chunks"

    : > "$CHUNKS_FILE"

    if [[ -n "$CHROMS_FILE" ]]; then
        CHROMS_FILE="$(realpath "$CHROMS_FILE")"
        if [[ ! -s "$CHROMS_FILE" ]]; then
            echo "ERROR: chroms file empty or missing: $CHROMS_FILE"
            exit 1
        fi
        CHROMS_SOURCE="$CHROMS_FILE"
    else
        CHROMS_SOURCE="$OUTPUT_DIR/panel.chroms.txt"
        "$BCFTOOLS" index -s "$PANEL_FIXED" | cut -f1 > "$CHROMS_SOURCE"
    fi

    while IFS= read -r CHR; do
        [[ -z "$CHR" ]] && continue

        echo "INFO: chunking chromosome: $CHR"
        CHUNK_TMP="$CHUNKS_DIR/chunks_${CHR}.txt"

        "$GLIMPSE2_CHUNK" \
            --input "$PANEL_FIXED" \
            --region "$CHR" \
            --window-mb "$WINDOW_MB" \
            --buffer-mb "$BUFFER_MB" \
            --sequential \
            --output "$CHUNK_TMP"

        if [[ ! -s "$CHUNK_TMP" ]]; then
            echo "WARN: no chunks for chromosome $CHR"
            continue
        fi

        cat "$CHUNK_TMP" >> "$CHUNKS_FILE"

    done < "$CHROMS_SOURCE"

    if [[ ! -s "$CHUNKS_FILE" ]]; then
        echo "ERROR: no chunks created"
        exit 1
    fi
fi

# ------------------------------------------------------------
# 5. Manifest
# ------------------------------------------------------------

MANIFEST="$OUTPUT_DIR/panel_manifest.tsv"

{
    echo -e "key\tvalue"
    echo -e "reference\t$REF_PATH"
    echo -e "input_panel\t$PANEL_IN"
    echo -e "work_panel\t$PANEL_WORK"
    echo -e "fixed_panel\t$PANEL_FIXED"
    echo -e "sites_vcf\t$SITES_VCF"
    echo -e "sites_tsv\t$SITES_TSV"
    echo -e "chunks\t$CHUNKS_FILE"
    echo -e "window_mb\t$WINDOW_MB"
    echo -e "buffer_mb\t$BUFFER_MB"
    echo -e "threads\t$THREADS"
    echo -e "date\t$(date -Is)"
} > "$MANIFEST"

echo "============================================================"
echo "Prepared GLIMPSE2 panel files:"
echo "Fixed panel: $PANEL_FIXED"
echo "Sites VCF:   $SITES_VCF"
echo "Sites TSV:   $SITES_TSV"
echo "Chunks:      $CHUNKS_FILE"
echo "Manifest:    $MANIFEST"
echo "============================================================"