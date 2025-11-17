#!/bin/bash
#SBATCH --job-name=star_run
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=60G            # <-- increase RAM significantly (e.g., 60GB)
#SBATCH --time=2:00:00
#SBATCH --output=star_run.log


THREADS=8
GENOME_DIR=/work/pi_mohammadariful_alam_uml_edu/nazim_macca_ref/STARindex
GTF=/home/nazimahmed_belabbaci_student_uml_edu/macca_ref/Macaca_mulatta.Mmul_10.108.gtf
SAMPLE_DIR=/work/pi_mohammadariful_alam_uml_edu/nazim_gse_samples
ALIGN_DIR=/work/pi_mohammadariful_alam_uml_edu/nazim_macca_ref/star_alignments
SJDB_OVERHANG=74

for SAMPLE in SRR19357304 SRR19357306 SRR19357307 SRR19357310 SRR19357321 SRR19357322 SRR19357323 SRR19357324 SRR19357325 SRR19357344 SRR19357346 SRR19357347 SRR19357350 SRR19357351 SRR19357352 SRR19357353 SRR19357354 SRR19357355; do
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
