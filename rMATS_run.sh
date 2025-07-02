#!/bin/bash
#SBATCH --job-name=rmats_baseline_vs_peak
#SBATCH --output=rmats_baseline_vs_peak.out
#SBATCH --error=rmats_baseline_vs_peak.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G


eval "$(conda shell.bash hook)"
conda activate rmats_env


BASE="/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SARS_Covid"
GTF="/home/nazimahmed_belabbaci_student_uml_edu/macca_ref/Macaca_mulatta.Mmul_10.108.gtf"
BAM_BASELINE="${BASE}/Baseline/trimmed_output/star_bam"
BAM_PEAK="${BASE}/Peak/trimmed_output/star_bam"
OUTDIR="${BASE}/rMATS_baseline_vs_peak"
TMPDIR="${OUTDIR}/tmp"


rm -rf "$OUTDIR"
mkdir -p "$OUTDIR" "$TMPDIR"

# generate comma-separated BAM lists and save to files
find "$BAM_BASELINE" -name '*Aligned.sortedByCoord.out.bam' | sort | paste -sd, - > "$OUTDIR/b1.txt"
find "$BAM_PEAK" -name '*Aligned.sortedByCoord.out.bam' | sort | paste -sd, - > "$OUTDIR/b2.txt"

# Log preview
echo "🧪 Baseline BAMs:"; cat "$OUTDIR/b1.txt"
echo "🧪 Peak BAMs:"; cat "$OUTDIR/b2.txt"

# Run rMATS
python ~/rmats-turbo/run_rmats.py \
  --b1 "$OUTDIR/b1.txt" \
  --b2 "$OUTDIR/b2.txt" \
  --gtf "$GTF" \
  --od "$OUTDIR" \
  --tmp "$TMPDIR" \
  -t single \
  --readLength 100 \
  --libType fr-unstranded \
  --nthread 8 \
  --bi 8
