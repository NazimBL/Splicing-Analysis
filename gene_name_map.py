import pandas as pd
import os

# === Setup ===
filtered_dir = 'rmats_output/filtered_events'
gtf_file = ('Macaca_mulatta.Mmul_10.108.gtf')
event_types = ['SE', 'RI', 'A5SS', 'A3SS', 'MXE']

# === Step 1: Build GeneID → GeneName map from GTF ===
gene_id_to_name = {}

with open(gtf_file, 'r') as gtf:
    for line in gtf:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[2] != 'gene':
            continue
        attributes = fields[8]
        gene_id = ''
        gene_name = ''
        for attr in attributes.split(';'):
            attr = attr.strip()
            if attr.startswith('gene_id'):
                gene_id = attr.split(' ')[1].strip('"')
            elif attr.startswith('gene_name'):
                gene_name = attr.split(' ')[1].strip('"')
        if gene_id and gene_name:
            gene_id_to_name[gene_id] = gene_name

print(f"Loaded {len(gene_id_to_name)} gene name mappings from GTF.")

# === Step 2: Add GeneName column to each *_gene_event_counts.csv ===
for event in event_types:
    input_file = os.path.join(filtered_dir, f'{event}_gene_event_counts.csv')
    output_file = os.path.join(filtered_dir, f'{event}_gene_event_counts_named.csv')

    if os.path.exists(input_file):
        df = pd.read_csv(input_file)
        df['GeneName'] = df['GeneID'].map(gene_id_to_name)
        df.to_csv(output_file, index=False)
        print(f"{event}: saved named version to {output_file}")
    else:
        print(f"{event}: file not found")
