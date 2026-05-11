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
    [--glimpse-panel panel.vcf.gz] \\
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

TRIMMED_DIR="$OUTPUT_DIR/1_trimmed"
SAMPLED_DIR="$OUTPUT_DIR/2_sampled"
ALIGNED_DIR="$OUTPUT_DIR/3_aligned"
MARKDUP_DIR="$OUTPUT_DIR/4_markdup"
FILTERED_DIR="$OUTPUT_DIR/5_filtered"
COVERAGE_DIR="$OUTPUT_DIR/6_coverage"
VCF_DIR="$OUTPUT_DIR/7_vcf"
IMPUTED_DIR="$OUTPUT_DIR/8_imputed"
METRICS_DIR="$OUTPUT_DIR/9_metrics"
RUN_META_DIR="$OUTPUT_DIR/00_run_metadata"

mkdir -p "$TRIMMED_DIR" "$SAMPLED_DIR" "$ALIGNED_DIR" \
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
# STAGE 2 - Trimming
# ============================================================

echo "STAGE 2 - Trimming"

TRIMMED_R1="$TRIMMED_DIR/${SAMPLE}_R1_trimmed.fastq.gz"
TRIMMED_R2="$TRIMMED_DIR/${SAMPLE}_R2_trimmed.fastq.gz"
TRIMMED_FR1="$TRIMMED_DIR/${SAMPLE}_R1_trimmed_fixed.fastq.gz"
TRIMMED_FR2="$TRIMMED_DIR/${SAMPLE}_R2_trimmed_fixed.fastq.gz"
SINGLETONE="$TRIMMED_DIR/${SAMPLE}_singletones.fastq.gz"

if [[ -f "$TRIMMED_R1" && -f "$TRIMMED_R2" ]]; then
    echo "SKIP: trimming already done for $SAMPLE"
else
    "$BBDUK" \
        in1="$R1" in2="$R2" \
        out1="$TRIMMED_R1" out2="$TRIMMED_R2" \
        ref="$BBDUK_REF_DIR/sequencing_artifacts.fa.gz,$BBDUK_REF_DIR/phix174_ill.ref.fa.gz,$BBDUK_REF_DIR/adapters.fa" \
        k=31 ordered cardinality \
        qtrim=rl trimq=20 maq=25 tbo mink=11 ktrim=r \
        minlen=50 \
        t="$BB_THREADS"

    "$REPAIR" \
        in1="$TRIMMED_R1" \
        in2="$TRIMMED_R2" \
        out1="$TRIMMED_FR1" \
        out2="$TRIMMED_FR2" \
        outsingle="$SINGLETONE" \
        threads="$BB_THREADS"

    mv "$TRIMMED_FR1" "$TRIMMED_R1"
    mv "$TRIMMED_FR2" "$TRIMMED_R2"
fi

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
        in1="$TRIMMED_R1" in2="$TRIMMED_R2" \
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
        -q 20 \
        -F 2820 \
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
CHUNKS_FILE="$IMPUTED_DIR_CHUNKS/chunks.txt"

if [[ "$RUN_GLIMPSE" -eq 0 ]]; then
    echo "SKIP: GLIMPSE2 disabled"
elif [[ -z "$GLIMPSE2_REF_PANEL" ]]; then
    echo "WARN: GLIMPSE2 panel not provided, skipping imputation"
else
    echo "STAGE 9 - Импутация (GLIMPSE2)"


  if [[ -f "$IMPUTED_VCF" && -f "${IMPUTED_VCF}.tbi" ]]; then
      echo "SKIP: импутация уже выполнена для $SAMPLE"
  else
      # ------------------------------------------------------------------
      # Проверяем наличие референсной панели
      # ------------------------------------------------------------------
      if [[ -z "$GLIMPSE2_REF_PANEL" || ! -f "$GLIMPSE2_REF_PANEL" ]]; then
          echo "ERROR: Для GLIMPSE2 необходима референсная панель (4-й аргумент)."
          exit 1
      fi

      if [[ ! -f "${GLIMPSE2_REF_PANEL}.tbi" ]]; then
          echo "INFO: Индекс референсной панели не найден, создаём..."
          "$TABIX" -p vcf "$GLIMPSE2_REF_PANEL"
      fi

      # ------------------------------------------------------------------
      # ДИАГНОСТИКА имён хромосом
      # ------------------------------------------------------------------
      echo "--- Диагностика имён хромосом ---"

      BAM_CHROMS_FILE="$IMPUTED_DIR/bam_chroms.txt"
      PANEL_CHROMS_FILE="$IMPUTED_DIR/panel_chroms.txt"

      "$SAMTOOLS" view -H "$FILTERED_BAM" \
          | grep "^@SQ" \
          | sed 's/.*SN:\([^\t]*\).*/\1/' \
          | sort > "$BAM_CHROMS_FILE"

      "$BCFTOOLS" index -s "$GLIMPSE2_REF_PANEL" \
          | cut -f1 \
          | sort > "$PANEL_CHROMS_FILE"

      echo "INFO: Хромосомы в BAM:"
      cat "$BAM_CHROMS_FILE"
      echo "---"
      echo "INFO: Хромосомы в панели:"
      cat "$PANEL_CHROMS_FILE"
      echo "---"

      # ------------------------------------------------------------------
      # Автоматическое определение и создание таблицы переименования
      # ------------------------------------------------------------------
      RENAME_TABLE="$IMPUTED_DIR/chrom_rename.txt"
      PANEL_RENAMED="$IMPUTED_DIR/ref_panel_renamed.vcf.gz"

      # Проверяем прямое совпадение
      DIRECT_COMMON=$(comm -12 \
          "$BAM_CHROMS_FILE" \
          "$PANEL_CHROMS_FILE" | wc -l)

      echo "INFO: Прямых совпадений хромосом: $DIRECT_COMMON"

      if [[ "$DIRECT_COMMON" -eq 0 ]]; then
          echo "INFO: Прямых совпадений нет — пробуем автоматическое переименование..."

          # Берём первую хромосому из каждого файла для определения паттерна
          BAM_FIRST=$(head -1 "$BAM_CHROMS_FILE")
          PANEL_FIRST=$(head -1 "$PANEL_CHROMS_FILE")

          echo "INFO: Пример BAM:   '$BAM_FIRST'"
          echo "INFO: Пример панели: '$PANEL_FIRST'"

          # Определяем направление переименования
          # Случай 1: панель = "Chr01", BAM = "glyma.Wm82.gnm4.Gm01" (соя)
          # Случай 2: панель = "1", BAM = "chr1"
          # Случай 3: панель = "chr1", BAM = "1"
          # Случай 4: панель = "Gm01", BAM = "glyma.Wm82.gnm4.Gm01"

          > "$RENAME_TABLE"

      while IFS= read -r PANEL_CHR; do
        BAM_MATCH=""

        # Попытка 1: точное совпадение последнего компонента (после последней точки)
        # glyma.Wm82.gnm4.Gm01 → суффикс "Gm01"
        PANEL_SUFFIX=$(echo "$PANEL_CHR" | rev | cut -d'.' -f1 | rev)

        if [[ -n "$PANEL_SUFFIX" ]]; then
          # Ищем строку где суффикс после последней точки совпадает точно
          while IFS= read -r B; do
            B_SUFFIX=$(echo "$B" | rev | cut -d'.' -f1 | rev)
            if [[ "$B_SUFFIX" == "$PANEL_SUFFIX" ]]; then
              BAM_MATCH="$B"
              break
            fi
          done < "$BAM_CHROMS_FILE"
        fi

        # Попытка 2: числовое ядро совпадает (для chr1 vs 1 vs Gm01)
        if [[ -z "$BAM_MATCH" ]]; then
          PANEL_NUM=$(echo "$PANEL_CHR" | grep -oP '\d+' | tail -1)
          if [[ -n "$PANEL_NUM" ]]; then
            # Нормализуем: убираем ведущие нули для сравнения
            PANEL_NUM_NORM=$(echo "$PANEL_NUM" | sed 's/^0*//')
            while IFS= read -r B; do
              B_NUM=$(echo "$B" | grep -oP '\d+' | tail -1)
              B_NUM_NORM=$(echo "$B_NUM" | sed 's/^0*//')
              if [[ -n "$B_NUM_NORM" && "$B_NUM_NORM" == "$PANEL_NUM_NORM" ]]; then
                BAM_MATCH="$B"
                break
              fi
            done < "$BAM_CHROMS_FILE"
          fi
        fi

        # Попытка 3: панель содержится как точный суффикс BAM (chr → chromosome)
        if [[ -z "$BAM_MATCH" ]]; then
          while IFS= read -r B; do
            if [[ "$B" == *".$PANEL_CHR" || "$B" == *"_${PANEL_CHR}" ]]; then
              BAM_MATCH="$B"
              break
            fi
          done < "$BAM_CHROMS_FILE"
        fi

        if [[ -n "$BAM_MATCH" ]]; then
          echo -e "${PANEL_CHR}\t${BAM_MATCH}" >> "$RENAME_TABLE"
        else
          echo "WARN: Не найдено соответствие для панельной хромосомы: $PANEL_CHR"
        fi

      done < "$PANEL_CHROMS_FILE"

          echo "INFO: Таблица переименования ($RENAME_TABLE):"
          cat "$RENAME_TABLE"

          if [[ ! -s "$RENAME_TABLE" ]]; then
              echo "ERROR: Не удалось автоматически сопоставить хромосомы."
              echo "ERROR: Создайте вручную файл $RENAME_TABLE"
              echo "ERROR: Формат (TSV): OLD_PANEL_NAME<TAB>BAM_NAME"
              echo "ERROR: Пример:"
              echo "ERROR:   Chr01<TAB>glyma.Wm82.gnm4.Gm01"
              echo "ERROR:   Chr02<TAB>glyma.Wm82.gnm4.Gm02"
              exit 1
          fi

          RENAMED_COUNT=$(wc -l < "$RENAME_TABLE")
          echo "INFO: Сопоставлено $RENAMED_COUNT хромосом"

          # Переименовываем хромосомы в референсной панели
      if [[ ! -f "$PANEL_RENAMED" || ! -f "${PANEL_RENAMED}.tbi" ]]; then
        echo "INFO: Переименование хромосом в референсной панели..."

        PANEL_RENAMED_UNSORTED="${PANEL_RENAMED%.vcf.gz}.unsorted.vcf.gz"

        "$BCFTOOLS" annotate \
          --threads "$THREADS" \
          --rename-chrs "$RENAME_TABLE" \
          "$GLIMPSE2_REF_PANEL" \
          -Oz -o "$PANEL_RENAMED_UNSORTED"

        if [[ $? -ne 0 || ! -f "$PANEL_RENAMED_UNSORTED" ]]; then
          echo "ERROR: bcftools annotate --rename-chrs завершился с ошибкой"
          exit 1
        fi

        # ВАЖНО: сортируем после переименования, т.к. порядок записей может нарушиться
        echo "INFO: Сортировка переименованной панели..."
        "$BCFTOOLS" sort \
          -Oz -o "$PANEL_RENAMED" \
          -T "$IMPUTED_DIR/bcftools_sort_tmp" \
          "$PANEL_RENAMED_UNSORTED"

        if [[ $? -ne 0 || ! -f "$PANEL_RENAMED" ]]; then
          echo "ERROR: bcftools sort завершился с ошибкой"
          exit 1
        fi

        rm -f "$PANEL_RENAMED_UNSORTED"

        "$BCFTOOLS" index -t --threads "$THREADS" "$PANEL_RENAMED"

        if [[ $? -ne 0 ]]; then
          echo "ERROR: bcftools index завершился с ошибкой для $PANEL_RENAMED"
          exit 1
        fi

        echo "INFO: Переименованная панель → $PANEL_RENAMED"

        # Проверяем результат
        echo "INFO: Хромосомы в переименованной панели:"
        "$BCFTOOLS" index -s "$PANEL_RENAMED" | cut -f1
      fi

          GLIMPSE2_REF_PANEL="$PANEL_RENAMED"

          # Обновляем список хромосом панели после переименования
          "$BCFTOOLS" index -s "$PANEL_RENAMED" \
              | cut -f1 \
              | sort > "$PANEL_CHROMS_FILE"

          # Финальная проверка совпадения
          FINAL_COMMON=$(comm -12 \
              "$BAM_CHROMS_FILE" \
              "$PANEL_CHROMS_FILE" | wc -l)

          echo "INFO: Совпадений после переименования: $FINAL_COMMON"

          if [[ "$FINAL_COMMON" -eq 0 ]]; then
              echo "ERROR: После переименования совпадений всё равно нет."
              echo "ERROR: Проверьте таблицу $RENAME_TABLE вручную."
              exit 1
          fi
      fi

      # ------------------------------------------------------------------
      # 9a. Подготовка референсной панели для GLIMPSE2
      # ------------------------------------------------------------------
      # Корневая причина ошибки "AC/AN INFO fields are inconsistent with GT":
      # GLIMPSE2 требует AN == 2*Nsamples на КАЖДОМ сайте панели (полная диплоидная
      # панель без missing). bcftools +fill-tags считает AN только по присутствующим
      # аллелям, поэтому даже один './.' даёт AN < 2N → mismatch → фатал.
      # Решение: сначала отбросить все сайты, где у любого образца есть './.',
      # затем фазировать, пересчитать AC/AN, удалить мономорфные сайты.
      # ------------------------------------------------------------------
      echo "--- 9a. Подготовка референсной панели для GLIMPSE2 ---"
      GLIMPSE2_REF_PANEL_FIXED="${GLIMPSE2_REF_PANEL%.vcf.gz}.fixed.vcf.gz"

      if [[ -f "$GLIMPSE2_REF_PANEL_FIXED" && -f "${GLIMPSE2_REF_PANEL_FIXED}.tbi" ]]; then
          echo "SKIP: исправленная панель уже существует → $GLIMPSE2_REF_PANEL_FIXED"
      else
          N_SAMPLES=$("$BCFTOOLS" query -l "$GLIMPSE2_REF_PANEL" | wc -l)
          EXPECTED_AN=$((N_SAMPLES * 2))
          echo "INFO: Образцов в панели: $N_SAMPLES (ожидаемое AN = $EXPECTED_AN)"

          echo "INFO: norm → biallelic SNP → fill missing GT with REF → phase → strip INFO → fill AC/AN → drop monomorphic"

          # 1. norm     — split multi-allelic, левосдвиг (важно для AC/AN consistency)
          # 2. view     — оставить только биаллельные SNP
          # 3. setGT . 0 — заменить ВСЕ missing-аллели (./. и .|.) на REF (0)
          #               (стандартный приём для imputation-панелей; не теряем сайты,
          #               где у некоторых образцов нет покрытия)
          # 4. setGT a p — все генотипы → phased (| вместо /)
          # 5. annotate — удалить все INFO и все FORMAT кроме GT
          # 6. fill-tags — пересчитать AC/AN/AF/NS с нуля; теперь AN == 2*Nsamples везде
          # 7. view     — убрать мономорфные сайты (AC==0 || AC==AN)
          # Стратегия:
          #   bcftools norm -m -any      — split multi-allelic. ПОБОЧКА: создаёт
          #     haploid '.' для сэмплов, у которых аллели не попали в данный
          #     биаллельный split. setGT '-n 0' тогда выдаёт haploid '0' → AN<2N.
          #   setGT -t . -n 0             — заполнить полностью missing diploid (./.)
          #   setGT -t ./x -n 0           — заполнить partially missing (./0, .|1)
          #   setGT -t a -n p             — phase all
          #   final view AN==EXP_AN       — отбросить редкие сайты с остаточным
          #     haploid (норм-побочка) и сайты, где даже после fill-tags AN<2N.
          EXP_AN="$EXPECTED_AN"
          "$BCFTOOLS" norm --threads "$THREADS" -f "$REF_PATH" -m -any "$GLIMPSE2_REF_PANEL" -Ou \
          | "$BCFTOOLS" view --threads "$THREADS" -m 2 -M 2 -v snps -Ou \
          | "$BCFTOOLS" +setGT -Ou -- -t . -n 0 \
          | "$BCFTOOLS" +setGT -Ou -- -t ./x -n 0 \
          | "$BCFTOOLS" +setGT -Ou -- -t a -n p \
          | "$BCFTOOLS" annotate --threads "$THREADS" -x 'INFO,^FORMAT/GT' -Ou \
          | "$BCFTOOLS" +fill-tags -Ou -- -t AC,AN,AF,NS \
          | "$BCFTOOLS" view --threads "$THREADS" \
              -e "AC==0 || AC==AN || INFO/AN!=$EXP_AN" \
              -Oz -o "$GLIMPSE2_REF_PANEL_FIXED"

          if [[ $? -ne 0 || ! -f "$GLIMPSE2_REF_PANEL_FIXED" ]]; then
              echo "ERROR: Не удалось создать исправленную панель"
              exit 1
          fi

          "$BCFTOOLS" index -t -f --threads "$THREADS" "$GLIMPSE2_REF_PANEL_FIXED"
          N_SITES=$("$BCFTOOLS" index -n "$GLIMPSE2_REF_PANEL_FIXED" 2>/dev/null \
                    || "$BCFTOOLS" view -H "$GLIMPSE2_REF_PANEL_FIXED" | wc -l)
          echo "INFO: Исправленная панель → $GLIMPSE2_REF_PANEL_FIXED ($N_SITES сайтов)"

          # ---- Sanity check: убеждаемся, что AN == 2*Nsamples везде ----
          echo "INFO: Проверка консистентности AC/AN..."
          # ВАЖНО: имя awk-переменной должно НЕ совпадать со встроенной функцией (exp, log, ...)
          N_BAD=$("$BCFTOOLS" query -f '%INFO/AN\n' "$GLIMPSE2_REF_PANEL_FIXED" \
                  | awk -v want="$EXPECTED_AN" '$1 != want' | wc -l)
          if [[ "$N_BAD" -gt 0 ]]; then
              echo "ERROR: $N_BAD сайтов имеют AN != $EXPECTED_AN"
              echo "ERROR: Примеры:"
              "$BCFTOOLS" query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' "$GLIMPSE2_REF_PANEL_FIXED" \
                  | awk -v want="$EXPECTED_AN" '$4 != want' | head -5
              exit 1
          fi
          echo "INFO: OK — все $N_SITES сайтов имеют AN == $EXPECTED_AN"

          # ---- Проверка фазированности (первые 10000 GT) ----
          N_UNPHASED=$("$BCFTOOLS" query -f '[%GT\n]' "$GLIMPSE2_REF_PANEL_FIXED" \
                       | head -10000 | grep -c '/' || true)
          if [[ "$N_UNPHASED" -gt 0 ]]; then
              echo "WARN: В панели осталось $N_UNPHASED нефазированных GT (из первых 10000)"
          else
              echo "INFO: Все генотипы успешно фазированы (проверка первых 10000)"
          fi
      fi

      GLIMPSE2_REF_PANEL_WORK="$GLIMPSE2_REF_PANEL_FIXED"

      # ------------------------------------------------------------------
      # 9b. Генерация GL по сайтам ref-панели
      # ------------------------------------------------------------------
      echo "--- 9b. Генерация GL по сайтам ref-панели ---"
      SAMPLE_GL_VCF="$IMPUTED_DIR/${SAMPLE}.gl.vcf.gz"

      if [[ -f "$SAMPLE_GL_VCF" && -f "${SAMPLE_GL_VCF}.tbi" ]]; then
          echo "SKIP: GL VCF уже создан для $SAMPLE"
      else
          SITES_VCF="$IMPUTED_DIR/ref_panel_sites.vcf.gz"
          SITES_TSV="$IMPUTED_DIR/ref_panel_sites.tsv.gz"

          if [[ ! -f "$SITES_VCF" || ! -f "${SITES_VCF}.tbi" ]]; then
              echo "INFO: Извлечение сайтов из ref-панели..."
              "$BCFTOOLS" view \
                  -G \
                  "$GLIMPSE2_REF_PANEL_WORK" \
                  -Oz -o "$SITES_VCF"
              "$BCFTOOLS" index -t --threads "$THREADS" "$SITES_VCF"
          fi

          if [[ ! -f "$SITES_TSV" || ! -f "${SITES_TSV}.tbi" ]]; then
              echo "INFO: Создание TSV сайтов..."
              "$BCFTOOLS" query \
                  -f '%CHROM\t%POS\t%REF,%ALT\n' \
                  "$SITES_VCF" \
              | "$BGZIP" -c > "$SITES_TSV"
              "$TABIX" -s1 -b2 -e2 "$SITES_TSV"
          fi

          # Добавляем ##contig строки из BAM в заголовок SITES_VCF
          SITES_VCF_REHEADERED="$IMPUTED_DIR/ref_panel_sites.reheadered.vcf.gz"

          if [[ ! -f "$SITES_VCF_REHEADERED" || ! -f "${SITES_VCF_REHEADERED}.tbi" ]]; then
              echo "INFO: Добавление ##contig строк из BAM в заголовок SITES_VCF..."

              CONTIG_HEADER="$IMPUTED_DIR/contig_header.txt"
              "$SAMTOOLS" view -H "$FILTERED_BAM" \
                  | grep "^@SQ" \
                  | awk '{
                      name=""; length=""
                      for(i=1;i<=NF;i++){
                          if($i ~ /^SN:/) name=substr($i,4)
                          if($i ~ /^LN:/) length=substr($i,4)
                      }
                      if(name!="" && length!="")
                          print "##contig=<ID="name",length="length">"
                  }' > "$CONTIG_HEADER"

              echo "INFO: Контигов из BAM: $(wc -l < "$CONTIG_HEADER")"

              # Собираем новый заголовок: мета-строки + новые контиги + #CHROM
              NEW_HEADER="$IMPUTED_DIR/new_header.txt"
              {
                  "$BCFTOOLS" view -h "$SITES_VCF" \
                      | grep -v "^##contig" \
                      | grep -v "^#CHROM"
                  cat "$CONTIG_HEADER"
                  "$BCFTOOLS" view -h "$SITES_VCF" | grep "^#CHROM"
              } > "$NEW_HEADER"

              "$BCFTOOLS" reheader \
                  -h "$NEW_HEADER" \
                  "$SITES_VCF" \
                  -o "$SITES_VCF_REHEADERED"

              if [[ $? -ne 0 || ! -f "$SITES_VCF_REHEADERED" ]]; then
                  echo "ERROR: bcftools reheader завершился с ошибкой"
                  exit 1
              fi

              "$BCFTOOLS" index -t --threads "$THREADS" "$SITES_VCF_REHEADERED"
              echo "INFO: Reheadered SITES_VCF → $SITES_VCF_REHEADERED"
          fi

          echo "INFO: Запуск bcftools mpileup | call..."

          "$BCFTOOLS" mpileup \
              --threads "$THREADS" \
              --fasta-ref "$REF_PATH" \
              --min-MQ 20 \
              --min-BQ 20 \
              -I \
              -E \
              -a 'FORMAT/DP,FORMAT/AD' \
              -T "$SITES_VCF_REHEADERED" \
              -Ou \
              "$FILTERED_BAM" \
          | "$BCFTOOLS" call \
              --threads "$THREADS" \
              -Aim \
              -C alleles \
              -T "$SITES_TSV" \
              -Oz -o "$SAMPLE_GL_VCF"

          EXIT_CODE=$?

          if [[ $EXIT_CODE -ne 0 || ! -f "$SAMPLE_GL_VCF" ]]; then
              echo "ERROR: bcftools mpileup|call завершился с кодом $EXIT_CODE"
              exit 1
          fi

          "$BCFTOOLS" index -t --threads "$THREADS" "$SAMPLE_GL_VCF"
          N_GL=$("$BCFTOOLS" view -H "$SAMPLE_GL_VCF" | wc -l)
          echo "INFO: GL VCF содержит $N_GL сайтов → $SAMPLE_GL_VCF"

          if [[ "$N_GL" -eq 0 ]]; then
              echo "ERROR: GL VCF пустой!"
              exit 1
          fi
      fi

      # ------------------------------------------------------------------
      # 9c. GLIMPSE2_chunk
      # ------------------------------------------------------------------
      echo "--- 9c. GLIMPSE2_chunk ---"

      if [[ -f "$CHUNKS_FILE" && -s "$CHUNKS_FILE" ]]; then
          echo "SKIP: чанки уже созданы → $CHUNKS_FILE"
      else
          > "$CHUNKS_FILE"

          # Используем только хромосомы общие для BAM и (переименованной) панели
          while IFS= read -r CHR; do
              echo "INFO: Чанкинг хромосомы $CHR"
              CHUNK_TMP="$IMPUTED_DIR_CHUNKS/chunks_${CHR}.txt"

              "$GLIMPSE2_CHUNK" \
                  --input "$GLIMPSE2_REF_PANEL_WORK" \
                  --region "$CHR" \
                  --window-mb 2 \
                  --buffer-mb 0.4 \
                  --sequential \
                  --output "$CHUNK_TMP"

              if [[ $? -ne 0 || ! -s "$CHUNK_TMP" ]]; then
                  echo "WARN: GLIMPSE2_chunk не создал чанки для $CHR, пропускаем"
                  continue
              fi

              cat "$CHUNK_TMP" >> "$CHUNKS_FILE"
              echo "INFO: $CHR: $(wc -l < "$CHUNK_TMP") чанков"

          done < <(comm -12 "$BAM_CHROMS_FILE" "$PANEL_CHROMS_FILE")

          if [[ ! -s "$CHUNKS_FILE" ]]; then
              echo "ERROR: Ни одного чанка не создано"
              exit 1
          fi

          echo "INFO: Всего чанков: $(wc -l < "$CHUNKS_FILE")"
      fi

      # ------------------------------------------------------------------
      # 9d. GLIMPSE2_phase
      # ------------------------------------------------------------------
      echo "--- 9d. GLIMPSE2_phase (по чанкам) ---"

      PHASE_DONE_FLAG="$IMPUTED_DIR_PHASED/.phase.done"

      # Адаптивные параметры PBWT под маленькую панель.
      # Default Kpbwt=2000 рассчитан на панель типа 1KG (5008 hap).
      # Если в нашей панели меньше 2000 гаплотипов, PBWT-селекция короткозамыкается
      # (Kpbwt >= n_ref_haps → "No PBWT selection" → "States for individual 0 are zero").
      # Решение: выставить Kpbwt и Kinit ниже числа гаплотипов в панели.
      N_REF_HAPS=$(( $("$BCFTOOLS" query -l "$GLIMPSE2_REF_PANEL_WORK" | wc -l) * 2 ))
      if (( N_REF_HAPS <= 2000 )); then
          KPBWT=$(( N_REF_HAPS / 2 ))
          KINIT=$(( N_REF_HAPS / 2 ))
          if (( KPBWT < 100 )); then KPBWT=100; fi
          if (( KINIT < 100 )); then KINIT=100; fi
          echo "INFO: панель имеет $N_REF_HAPS гаплотипов; используем --Kinit $KINIT --Kpbwt $KPBWT"
      else
          KPBWT=2000
          KINIT=1000
      fi

      if [[ -f "$PHASE_DONE_FLAG" ]]; then
          echo "SKIP: фазирование чанков уже выполнено"
      else
          while IFS=$'\t' read -r ID CHR IRG ORG REST; do
              [[ -z "$ID" || "$ID" == \#* ]] && continue

              # Chunk ID нумеруется ПЕР-хромосомно (0,1,...) → имя файла должно
              # включать хромосому, иначе chunk0 из Gm02 затрёт chunk0 из Gm01.
              PHASED_CHUNK="$IMPUTED_DIR_PHASED/${SAMPLE}.${CHR}.chunk${ID}.bcf"

              if [[ -f "$PHASED_CHUNK" && -f "${PHASED_CHUNK}.csi" ]]; then
                  echo "SKIP: чанк $ID ($IRG) уже обработан"
                  continue
              fi

              echo "INFO: Фазирование чанка $ID — input: $IRG  output: $ORG"

              "$GLIMPSE2_PHASE" \
                  --input-gl "$SAMPLE_GL_VCF" \
                  --reference "$GLIMPSE2_REF_PANEL_WORK" \
                  --input-region "$IRG" \
                  --output-region "$ORG" \
                  --output "$PHASED_CHUNK" \
                  --Kinit "$KINIT" \
                  --Kpbwt "$KPBWT" \
                  --threads "$THREADS"

              EXIT_CODE=$?

              if [[ $EXIT_CODE -ne 0 || ! -f "$PHASED_CHUNK" ]]; then
                  echo "ERROR: GLIMPSE2_phase завершился с кодом $EXIT_CODE"
                  echo "ERROR: ID=$ID  IRG=$IRG  ORG=$ORG"
                  exit 1
              fi

              "$BCFTOOLS" index -f "$PHASED_CHUNK"
              echo "INFO: Чанк $ID → $PHASED_CHUNK"

          done < "$CHUNKS_FILE"

          touch "$PHASE_DONE_FLAG"
          echo "INFO: Все чанки успешно обработаны"
      fi

      # ------------------------------------------------------------------
      # 9e. GLIMPSE2_ligate
      # ------------------------------------------------------------------
      echo "--- 9e. GLIMPSE2_ligate ---"

      PHASED_LIST="$IMPUTED_DIR_PHASED/phased_chunks.list"
      # Сортируем по (CHR, ID) — нужно для GLIMPSE2_ligate (порядок чанков должен
      # совпадать с порядком в исходном CHUNKS_FILE).
      while IFS=$'\t' read -r ID CHR _ ; do
          [[ -z "$ID" || "$ID" == \#* ]] && continue
          f="$IMPUTED_DIR_PHASED/${SAMPLE}.${CHR}.chunk${ID}.bcf"
          [[ -f "$f" ]] && echo "$f"
      done < "$CHUNKS_FILE" > "$PHASED_LIST"

      if [[ ! -s "$PHASED_LIST" ]]; then
          echo "ERROR: Не найдены фазированные чанки"
          exit 1
      fi

      echo "INFO: Сборка $(wc -l < "$PHASED_LIST") чанков (GLIMPSE2_ligate)"

      "$GLIMPSE2_LIGATE" \
          --input "$PHASED_LIST" \
          --output "$IMPUTED_VCF" \
          --threads 4

      if [[ $? -ne 0 || ! -f "$IMPUTED_VCF" ]]; then
          echo "ERROR: GLIMPSE2_ligate не создал $IMPUTED_VCF"
          exit 1
      fi

      "$BCFTOOLS" index -t --threads "$THREADS" "$IMPUTED_VCF"
      echo "INFO: Импутированный VCF → $IMPUTED_VCF"
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