#!/bin/bash
#SBATCH --job-name=download_sra_3d
#SBATCH --output=download_sra_3d.out
#SBATCH --error=download_sra_3d.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G            # Slightly less than total available to accommodate overhead


conda activate sra_download_env

# Directory to save downloaded samples
output_dir="/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SIV/Day3"
mkdir -p "$output_dir"
cd "$output_dir"

# Path to SRR list file
srr_list="/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SIV/download_list_day3.txt"


threads=8

# loop over SRR list and download each sample
while read srr; do
    echo "Downloading $srr..."
    fasterq-dump "$srr" --split-files --threads "$threads"
done < "$srr_list"

echo "All downloads complete."

