#!/bin/bash
#SBATCH --job-name=week24_trim_qc_pe
#SBATCH --output=week24_trim_qc_pe.out
#SBATCH --error=week24_trim_qc_pe.err
#SBATCH --time=28:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G

conda activate rnaseq_env

baseline_dir="/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SIV/Week24"
trimmed_dir="${baseline_dir}/trimmed_output"
qc_dir="${trimmed_dir}/qc_results"

# output directories
mkdir -p "$trimmed_dir"
mkdir -p "$qc_dir"

# Path to adapter file
adapter_file="${HOME}/macca_ref/Scripts/TruSeq3-PE.fa"

# Download adapter file if not present
if [ ! -f "$adapter_file" ]; then
    curl -o "$adapter_file" https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa
fi

# Loop through all *_1.fastq files
for fq1 in ${baseline_dir}/*_1.fastq; do
    sample=$(basename "$fq1" _1.fastq)
    fq2="${baseline_dir}/${sample}_2.fastq"

    echo "Processing $sample..."

    # Run Trimmomatic in PE mode
    trimmomatic PE -threads 8 \
        "$fq1" "$fq2" \
        "${trimmed_dir}/${sample}_1_paired.fq" "${trimmed_dir}/${sample}_1_unpaired.fq" \
        "${trimmed_dir}/${sample}_2_paired.fq" "${trimmed_dir}/${sample}_2_unpaired.fq" \
        ILLUMINACLIP:${adapter_file}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # Run FastQC on the paired output files
    fastqc "${trimmed_dir}/${sample}_1_paired.fq" -o "$qc_dir"
    fastqc "${trimmed_dir}/${sample}_2_paired.fq" -o "$qc_dir"
done

echo "Paired-end trimming and QC completed."


