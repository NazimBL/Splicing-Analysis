import pandas as pd
import os

# Path to your rMATS output folder
rmats_output_dir = 'rMATS_baseline_vs_resolution'

# Event types you want to process
event_types = ['SE', 'RI', 'A5SS', 'A3SS', 'MXE']

# Output folder
filtered_output_dir = os.path.join(rmats_output_dir, 'filtered_events')
os.makedirs(filtered_output_dir, exist_ok=True)

# Filtering thresholds
fdr_threshold = 0.05
delta_psi_threshold = 0

for event in event_types:
    input_file = os.path.join(rmats_output_dir, f'{event}.MATS.JC.txt')
    output_file = os.path.join(filtered_output_dir, f'{event}_filtered.csv')

    # Load the rMATS file
    df = pd.read_csv(input_file, sep='\t')

    # Filter based on FDR and IncLevelDifference (ΔΨ)
    filtered_df = df[(df['FDR'] < fdr_threshold) & (df['IncLevelDifference'].abs() > delta_psi_threshold)]

    # Sort by absolute splicing change
    filtered_df = filtered_df.reindex(filtered_df['IncLevelDifference'].abs().sort_values(ascending=False).index)

    # Save the filtered file
    filtered_df.to_csv(output_file, index=False)

    print(f'{event}: {len(filtered_df)} events passed filters and saved to {output_file}')