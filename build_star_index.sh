#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=60G            # <-- increase RAM significantly (e.g., 60GB)
#SBATCH --time=02:00:00
#SBATCH --output=star_index.log


cd /work/pi_mohammadariful_alam_uml_edu
mkdir -p nazim_macca_ref
cd nazim_macca_ref

STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir STARindex \
  --genomeFastaFiles /home/nazimahmed_belabbaci_student_uml_edu/macca_ref/Macaca_mulatta.Mmul_10.dna.toplevel.fa \
  --sjdbGTFfile /home/nazimahmed_belabbaci_student_uml_edu/macca_ref/Macaca_mulatta.Mmul_10.108.gtf \
  --sjdbOverhang 74
