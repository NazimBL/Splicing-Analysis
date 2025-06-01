cd ~/macca_ref/GSE_samples

for FQ in *_paired.fq; do
  RAW=$(head -n 2 "$FQ" | tail -n 1 | wc -c)
  ACTUAL_BASES=$((RAW - 1))
  printf "%s\t%d bases\n" "$FQ" "$ACTUAL_BASES"
done | sort -k2n

