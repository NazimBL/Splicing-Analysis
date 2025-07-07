#!/bin/bash
#SBATCH --job-name=star_align_se
#SBATCH --output=star_align_se.log
#SBATCH --error=star_align_se.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G


# STAR index directory
STAR_INDEX="/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SARS_Covid/STARindex"

# directory for all timepoints
BASE_DIR="/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SARS_Covid"
TIMEPOINTS=("Baseline" "Peak" "Resolution")

for TIMEPOINT in "${TIMEPOINTS[@]}"; do
    echo "🔄 Aligning $TIMEPOINT samples..."

    INPUT_DIR="${BASE_DIR}/${TIMEPOINT}/trimmed_output"
    OUTPUT_DIR="${INPUT_DIR}/star_bam"
    mkdir -p "$OUTPUT_DIR"

    for FQ in "${INPUT_DIR}"/*_trimmed.fq; do
        [ -e "$FQ" ] || continue  # skip if no matching files

        SAMPLE=$(basename "$FQ" _trimmed.fq)
        OUT_PREFIX="${OUTPUT_DIR}/${SAMPLE}_"

        STAR --runThreadN 8 \
            --genomeDir "$STAR_INDEX" \
            --readFilesIn "$FQ" \
            --outFileNamePrefix "$OUT_PREFIX" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM MD \
            --outSAMunmapped Within \
            --outFilterMultimapNmax 1 \
            --alignEndsType EndToEnd
    done

    echo "Finished $TIMEPOINT"
done

echo "All alignments completed."
