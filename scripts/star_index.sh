#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=60G
#SBATCH --time=02:00:00
#SBATCH --output=star_index.log
#SBATCH --error=star_index.err

# Activate STAR environment
conda activate sra_download_env

# Set working directory for the index
index_dir="/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SIV/STARindex"
mkdir -p "$index_dir"

# Run STAR index generation
STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir "$index_dir" \
  --genomeFastaFiles /home/nazimahmed_belabbaci_student_uml_edu/macca_ref/Macaca_mulatta.Mmul_10.dna.toplevel.fa \
  --sjdbGTFfile /home/nazimahmed_belabbaci_student_uml_edu/macca_ref/Macaca_mulatta.Mmul_10.108.gtf \
  --sjdbOverhang 74 \
  --limitGenomeGenerateRAM 50000000000 \
  --genomeSAsparseD 2

