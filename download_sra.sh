#!/bin/bash
#SBATCH --job-name=download_sra_peak
#SBATCH --output=download_sra_peak.out
#SBATCH --error=download_sra_peak.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G            # Slightly less than total available to accommodate overhead


conda activate sra_download_env

# Directory to save downloaded samples
output_dir="/work/pi_mohammadariful_alam_uml_edu/nazim_SARS_Covid/Peak"
mkdir -p "$output_dir"
cd "$output_dir"

# Path to SRR list file
srr_list="/work/pi_mohammadariful_alam_uml_edu/nazim_SARS_Covid/download_list_peak.txt"


threads=8

# loop over SRR list and download each sample
while read srr; do
    echo "Downloading $srr..."
    fasterq-dump "$srr" --split-files --threads "$threads"
done < "$srr_list"

echo "All downloads complete."
