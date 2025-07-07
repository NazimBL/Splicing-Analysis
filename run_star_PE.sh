#!/usr/bin/env bash

THREADS=8
GENOME_DIR=/work/pi_mohammadariful_alam_uml_edu/nazim_macca_ref/STARindex
GTF=/home/nazimahmed_belabbaci_student_uml_edu/macca_ref/Macaca_mulatta.Mmul_10.108.gtf
SAMPLE_DIR=/home/nazimahmed_belabbaci_student_uml_edu/macca_ref/GSE_samples
ALIGN_DIR=/work/pi_mohammadariful_alam_uml_edu/nazim_macca_ref/star_alignments
SJDB_OVERHANG=74

for SAMPLE in SRR19357323 SRR19357324 SRR19357325 SRR19357344 SRR19357346 SRR19357347; do
  echo "=== Aligning $SAMPLE ==="

  mkdir -p ${ALIGN_DIR}/${SAMPLE}
  cd ${ALIGN_DIR}/${SAMPLE}

  # 1st-pass alignment:
  STAR \
    --runThreadN ${THREADS} \
    --genomeDir ${GENOME_DIR} \
    --readFilesIn ${SAMPLE_DIR}/${SAMPLE}_1_paired.fq \
                  ${SAMPLE_DIR}/${SAMPLE}_2_paired.fq \
    --readFilesCommand cat \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang ${SJDB_OVERHANG} \
    --outFileNamePrefix ${SAMPLE}.pass1. \
    --outSAMtype BAM Unsorted \
    --outSAMstrandField intronMotif

  PASS1_SJ="${ALIGN_DIR}/${SAMPLE}/${SAMPLE}.pass1.SJ.out.tab"

  # 2nd-pass alignment:
  STAR \
    --runThreadN ${THREADS} \
    --genomeDir ${GENOME_DIR} \
    --readFilesIn ${SAMPLE_DIR}/${SAMPLE}_1_paired.fq \
                  ${SAMPLE_DIR}/${SAMPLE}_2_paired.fq \
    --readFilesCommand cat \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang ${SJDB_OVERHANG} \
    --sjdbFileChrStartEnd ${PASS1_SJ} \
    --outFileNamePrefix ${SAMPLE}. \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI AS nM

  samtools index ${SAMPLE}.Aligned.sortedByCoord.out.bam

  echo "=== Done with $SAMPLE ==="
done

