set -euo pipefail
PANEL="/mnt/low/soybean_vcf_without_208_211_45_samples.vcf.gz"
REGION="glyma.Wm82.gnm4.Gm01:1-4354142"
BCFTOOLS="$(command -v bcftools || true)"
if [[ -z "$BCFTOOLS" ]]; then
  if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/bcftools" ]]; then
    BCFTOOLS="${CONDA_PREFIX}/bin/bcftools"
  else
    echo "ERROR: bcftools not found"
    exit 127
  fi
fi

echo "BCFTOOLS=$BCFTOOLS"
echo "PANEL=$PANEL"
echo "REGION=$REGION"

echo "== panel sample count =="
"$BCFTOOLS" query -l "$PANEL" | wc -l

echo "== AC/AN header definitions =="
"$BCFTOOLS" view -h "$PANEL" | grep -E '^##INFO=<ID=(AC|AN),' || true

echo "== first 10 records with stored AC/AN and recomputed AC/AN from GT =="
TMP="$(mktemp --suffix=.bcf)"
"$BCFTOOLS" +fill-tags "$PANEL" -r "$REGION" -Ou -- -t AC,AN \
| "$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' \
| head -10

echo "== first 20 sites where stored AC/AN differ after recomputing from GT =="
"$BCFTOOLS" view -r "$REGION" -Ou "$PANEL" \
| "$BCFTOOLS" annotate -x INFO/AC,INFO/AN -Ou \
| "$BCFTOOLS" +fill-tags -Ou -- -t AC,AN \
| "$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' \
| head -20

echo "== original first 20 sites for comparison =="
"$BCFTOOLS" query -r "$REGION" -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' "$PANEL" | head -20
