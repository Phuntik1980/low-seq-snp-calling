#!/bin/bash
# Скрипт для извлечения образца из панели и приведения к формату truth VCF

set -e  # остановиться при любой ошибке

SAMPLE_IN_PANEL="211"           # имя образца в панели
SAMPLE_FINAL="211_S16"          # желаемое имя образца в выходном файле
PANEL="/mnt/low/annotation_soybean.vcf.gz"
REF="/mnt/low/genome/glyma.Wm82.gnm4.4PTR.genome_main.fna"
OUT="211_S16.truth.vcf.gz"

echo "=== Извлечение образца $SAMPLE_IN_PANEL из панели ==="
bcftools view -s "$SAMPLE_IN_PANEL" "$PANEL" -Oz -o temp_211.vcf.gz

echo "=== Переименование образца в $SAMPLE_FINAL ==="
bcftools reheader -s <(echo "$SAMPLE_FINAL") temp_211.vcf.gz -o temp_211_renamed.vcf.gz

echo "=== Нормализация и фильтрация (только биаллельные SNPs) ==="
bcftools norm -m -any -f "$REF" temp_211_renamed.vcf.gz -Ou | \
bcftools view -m2 -M2 -v snps -Oz -o "$OUT"

echo "=== Индексация ==="
bcftools index -t "$OUT"

echo "=== Очистка временных файлов ==="
rm -f temp_211.vcf.gz temp_211_renamed.vcf.gz

echo "=== Готово: $OUT ==="
echo "Число SNP в выходном файле:"
bcftools view -H "$OUT" | wc -l
