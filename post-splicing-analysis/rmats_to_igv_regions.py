import pandas as pd

# Change this to any of your fromGTF files
input_file = 'rmats_output/fromGTF.SE.txt'
output_file = 'SE_igv_regions.bed'

# Load the rMATS fromGTF file
df = pd.read_csv(input_file, sep='\t')

# Calculate the full region covering the exon skipping event:
# Add 100bp upstream/downstream buffer for IGV context
def region_with_buffer(row, buffer=100):
    chrom = row['chr']
    start = max(0, min(row['upstreamES'], row['downstreamEE']) - buffer)
    end = max(row['upstreamES'], row['downstreamEE']) + buffer
    return f"{chrom}:{start}-{end}"

# Apply and save
df['region'] = df.apply(region_with_buffer, axis=1)
df[['region']].to_csv(output_file, index=False, header=False)

print(f"✅ Saved IGV region list to: {output_file}")
