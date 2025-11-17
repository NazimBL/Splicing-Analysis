#!/bin/bash
#SBATCH --job-name=6h_star_align_pe
#SBATCH --output=pe_6h.log
#SBATCH --error=pe_6h.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G


THREADS=8 GENOME_DIR=/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SIV/STARindex 
GTF=/home/nazimahmed_belabbaci_student_uml_edu/macca_ref/Macaca_mulatta.Mmul_10.108.gtf 
SAMPLE_DIR=/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SIV/6hours/trimmed_output 
ALIGN_DIR=/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SIV/6hours/star_alignments 
SJDB_OVERHANG=74

mkdir -p "$ALIGN_DIR"

for FQ1 in ${SAMPLE_DIR}/*_1_paired.fq; do
  SAMPLE=$(basename "$FQ1" _1_paired.fq)
  FQ2="${SAMPLE_DIR}/${SAMPLE}_2_paired.fq"
  echo "=== Aligning $SAMPLE ==="

  mkdir -p ${ALIGN_DIR}/${SAMPLE}
  cd ${ALIGN_DIR}/${SAMPLE}

  # First pass
  STAR \
    --runThreadN ${THREADS} \
    --genomeDir ${GENOME_DIR} \
    --readFilesIn "$FQ1" "$FQ2" \
    --readFilesCommand cat \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang ${SJDB_OVERHANG} \
    --outFileNamePrefix ${SAMPLE}.pass1. \
    --outSAMtype BAM Unsorted \
    --outSAMstrandField intronMotif

  PASS1_SJ="${ALIGN_DIR}/${SAMPLE}/${SAMPLE}.pass1.SJ.out.tab"

  # Second pass
  STAR \
    --runThreadN ${THREADS} \
    --genomeDir ${GENOME_DIR} \
    --readFilesIn "$FQ1" "$FQ2" \
    --readFilesCommand cat \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang ${SJDB_OVERHANG} \
    --sjdbFileChrStartEnd "$PASS1_SJ" \
    --outFileNamePrefix ${SAMPLE}. \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI AS nM MD

  # Index the BAM
  samtools index ${SAMPLE}.Aligned.sortedByCoord.out.bam

  echo "=== Done with $SAMPLE ==="
done

