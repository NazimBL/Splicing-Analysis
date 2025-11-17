import pandas as pd
from collections import defaultdict
import gzip

def parse_gtf_for_exons(gtf_path):
    exon_map = defaultdict(list)  # {(chr, start, end, strand): [transcript_id]}

    open_func = gzip.open if gtf_path.endswith('.gz') else open
    with open_func(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) != 9: continue
            seqname, source, feature, start, end, score, strand, frame, attributes = fields
            if feature != 'exon': continue

            # Extract transcript_id
            attrs = {item.strip().split(' ')[0]: item.strip().split(' ')[1].replace('"', '') 
                     for item in attributes.strip().split(';') if item.strip()}
            transcript_id = attrs.get('transcript_id', None)
            if transcript_id:
                exon_map[(seqname, int(start), int(end), strand)].append(transcript_id)

    return exon_map

def match_rmats_to_transcripts(rmats_path, exon_map, output_csv='rmats_with_transcripts.csv'):
    df = pd.read_csv(rmats_path)
    df.columns = df.columns.str.strip().str.lower()

    results = []

    results = []

    for _, row in df.iterrows():
        # Then access columns like this:
        chrom = row['chr'].replace('chr', '')
        strand = row['strand']
        start = int(row['exonstart_0base']) + 1
        end = int(row['exonend'])

        key = (chrom, start, end, strand)
        transcripts = exon_map.get(key, [])
        results.append({
            'chr': row['chr'],
            'start': start,
            'end': end,
            'strand': strand,
            'geneSymbol': row.get('geneSymbol', ''),
            'transcripts': ';'.join(transcripts),
            'IncLevelDifference': row.get('IncLevelDifference', ''),
            'FDR': row.get('FDR', '')
        })

    out_df = pd.DataFrame(results)
    out_df.to_csv(output_csv, index=False)
    print(f"Saved: {output_csv}")

# Example usage:
exon_map = parse_gtf_for_exons("Macaca_mulatta.Mmul_10.108.gtf")
match_rmats_to_transcripts("SE.MATS.JC.txt", exon_map)
