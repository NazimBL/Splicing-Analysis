#!/bin/bash
#SBATCH --job-name=base_vs_24
#SBATCH --output=rmats_base_vs_24.out
#SBATCH --error=rmats_base_vs_24.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G

# Conda setup
eval "$(conda shell.bash hook)"
conda activate rmats_env

# Define directories
BASE="/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SIV"
GTF="/home/nazimahmed_belabbaci_student_uml_edu/macca_ref/Macaca_mulatta.Mmul_10.108.gtf"
BAM_BASELINE="${BASE}/Baseline/star_alignments"
BAM_PEAK="${BASE}/Week24/star_alignments"
OUTDIR="${BASE}/rMATS_Base_vs_W24"
TMPDIR="${OUTDIR}/tmp"

# Clean previous output if exists
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR" "$TMPDIR"

# Generate comma-separated BAM lists
find "$BAM_BASELINE" -name '*Aligned.sortedByCoord.out.bam' | sort | paste -sd, - > "$OUTDIR/b1.txt"
find "$BAM_PEAK" -name '*Aligned.sortedByCoord.out.bam' | sort | paste -sd, - > "$OUTDIR/b2.txt"

# Preview files
echo "🧪 Baseline BAMs:"; cat "$OUTDIR/b1.txt"
echo "🧪 Peak BAMs:"; cat "$OUTDIR/b2.txt"

# Run rMATS for paired-end
python ~/rmats-turbo/run_rmats.py \
  --b1 "$OUTDIR/b1.txt" \
  --b2 "$OUTDIR/b2.txt" \
  --gtf "$GTF" \
  --od "$OUTDIR" \
  --tmp "$TMPDIR" \
  -t paired \
  --readLength 75 \
  --libType fr-unstranded \
  --nthread 8
