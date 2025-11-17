import pandas as pd
import os

# Path to your rMATS output folder
rmats_output_dir = 'rmats_output'

# Event types you want to process
event_types = ['SE', 'RI', 'A5SS', 'A3SS', 'MXE']

# Output folder
filtered_output_dir = os.path.join(rmats_output_dir, 'filtered_events')
os.makedirs(filtered_output_dir, exist_ok=True)

# Filtering thresholds
fdr_threshold = 0.5
delta_psi_threshold = 0.01

for event in event_types:
    input_file = os.path.join(rmats_output_dir, f'{event}.MATS.JC.txt')
    output_file = os.path.join(filtered_output_dir, f'{event}_filtered.csv')

    # Load the rMATS file
    df = pd.read_csv(input_file, sep='\t')

    print(f"\n=== {event} ===")
    print(f"Total events: {len(df)}")
    print("Columns:", df.columns.tolist())

    # Basic FDR and ΔΨ stats
    fdr_pass = df[df['FDR'] < fdr_threshold]
    psi_pass = df[df['IncLevelDifference'].abs() > delta_psi_threshold]

    print(f"Events with FDR < {fdr_threshold}: {len(fdr_pass)}")
    print(f"Events with |ΔΨ| > {delta_psi_threshold}: {len(psi_pass)}")
    print(f"Events passing both filters: {len(df[(df['FDR'] < fdr_threshold) & (df['IncLevelDifference'].abs() > delta_psi_threshold)])}")

    # Optional: Show top 5 entries with lowest FDR
    print(df[['FDR', 'IncLevelDifference']].sort_values(by='FDR').head())

    # Apply filtering
    filtered_df = df[(df['FDR'] < fdr_threshold) & (df['IncLevelDifference'].abs() > delta_psi_threshold)]
    filtered_df = filtered_df.reindex(filtered_df['IncLevelDifference'].abs().sort_values(ascending=False).index)

    # Save filtered results
    filtered_df.to_csv(output_file, index=False)
    print(f"Saved {len(filtered_df)} events to {output_file}")
