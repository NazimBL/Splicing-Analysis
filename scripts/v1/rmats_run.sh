#!/bin/bash
#SBATCH --job-name=rmats_run
#SBATCH --output=rmats_run.out
#SBATCH --error=rmats_run.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=60G
#SBATCH --time=4:00:00

# Manually activate conda environm  # avoid `conda activate` to prevent CondaError

# Change working directory to where you want outputs
cd /work/pi_mohammadariful_alam_uml_edu/nazim_macca_ref/rmats_output

# Run rMATS
python ~/rmats-turbo/run_rmats.py \
--b1 /work/pi_mohammadariful_alam_uml_edu/nazim_macca_ref/baseline_bams.txt \
--b2 /work/pi_mohammadariful_alam_uml_edu/nazim_macca_ref/week24_bams.txt \
  --gtf /home/nazimahmed_belabbaci_student_uml_edu/macca_ref/Macaca_mulatta.Mmul_10.108.gtf \
  --od . \
  --tmp ./tmp_rmats \
  --readLength 65 \
  --libType fr-firststrand \
  --nthread 8 \
  --bi 8
