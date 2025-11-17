#!/bin/bash
#SBATCH --job-name=trim_qc_peak
#SBATCH --output=trim_qc_peak.out
#SBATCH --error=trim_qc_peak.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G

conda activate rnaseq_env

# Set paths
baseline_dir="/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SARS_Covid/Peak"
trimmed_dir="${baseline_dir}/trimmed_output"
qc_dir="${trimmed_dir}/qc_results"

# Make output directories if they don't exist
mkdir -p "$trimmed_dir"
mkdir -p "$qc_dir"

# Path to adapter file
adapter_file="${HOME}/macca_ref/Scripts/TruSeq3-PE.fa"

# Download adapter file if not already present
if [ ! -f "$adapter_file" ]; then
    curl -o "$adapter_file" https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa
fi

# Loop through all *_1.fastq files in the Baseline folder
for fq in ${baseline_dir}/*_1.fastq; do
    sample=$(basename "$fq" _1.fastq)
    echo "Processing $sample..."

    # Trim using Trimmomatic SE mode
    trimmomatic SE -threads 8 \
        "$fq" "${trimmed_dir}/${sample}_trimmed.fq" \
        ILLUMINACLIP:${adapter_file}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # Run FastQC on trimmed file
    fastqc "${trimmed_dir}/${sample}_trimmed.fq" -o "$qc_dir"
done

echo "Trimming and quality check completed."
