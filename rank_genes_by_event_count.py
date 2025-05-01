import pandas as pd
import os

# Path to your filtered events folder
filtered_dir = 'rmats_output/filtered_events'
event_types = ['SE', 'RI', 'A5SS', 'A3SS', 'MXE']

for event in event_types:
    input_file = os.path.join(filtered_dir, f'{event}_filtered.csv')
    output_file = os.path.join(filtered_dir, f'{event}_gene_event_counts.csv')

    if os.path.exists(input_file):
        df = pd.read_csv(input_file)

        if 'GeneID' in df.columns:
            # Count how many events per gene
            gene_counts = df['GeneID'].value_counts()

            # Save to CSV
            gene_counts.to_csv(output_file, header=["event_count"])

            # Preview top 5
            print(f"\nTop genes with most {event} events:")
            print(gene_counts.head(5))
        else:
            print(f"{event}: 'GeneID' column not found in {input_file}")
    else:
        print(f"{event}: file not found")
