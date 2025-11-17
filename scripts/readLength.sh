#!/bin/bash
#SBATCH --job-name=read_length_mode
#SBATCH --output=read_length_mode.out
#SBATCH --error=read_length_mode.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

conda activate sra_download_env

# base path
base_dir="/scratch/workspace/nazimahmed_belabbaci_student_uml_edu-SplicingWorkspace/nazim_SIV/"
# folders  ("Baseline" "Peak" "Resolution")
folders=("Baseline" "Day1" "Day2" "Day3" "6hours" "Week24")
echo "Calculating mode read lengths for each folder..."

for folder in "${folders[@]}"; do
    echo "Processing $folder..."
	#trimmed folder
    fq_path="${base_dir}/${folder}/trimmed_output"

    if [ ! -d "$fq_path" ]; then
        echo "  Folder not found: $fq_path"
        continue
    fi

    # Temp file to hold all read lengths
    length_file=$(mktemp)

    # Loop over all *_1.fastq files
    for fq in "$fq_path"/*.fq; do
        if [ ! -f "$fq" ]; then continue; fi

        # Read length of the first sequence
        length=$(head -n 2 "$fq" | tail -n 1 | wc -c)
        length=$((length - 1))  # remove newline char

        echo "$length" >> "$length_file"
    done

    echo -n " Mode read length in $folder: "

    # Compute mode using sort + uniq + awk
    sort "$length_file" | uniq -c | sort -nr | head -n 1 | awk '{print $2 " bp (" $1 " occurrences)"}'

    rm "$length_file"
done

echo "Done."
