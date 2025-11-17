
import pandas as pd
from collections import defaultdict
import gzip

def parse_gtf_exons(gtf_file):
    exon_lookup = defaultdict(list)  # key: (chr, start, end, strand), value: list of transcript_ids

    open_func = gzip.open if gtf_file.endswith('.gz') else open
    with open_func(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) != 9: continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            if feature != 'exon': continue

            attrs = {}
            for attr in attributes.strip().split(';'):
                if attr.strip():
                    key, value = attr.strip().split(' ', 1)
                    attrs[key] = value.replace('"', '').strip()

            transcript_id = attrs.get('transcript_id')
            if transcript_id:
                exon_lookup[(chrom.replace('chr', ''), int(start), int(end), strand)].append(transcript_id)

    return exon_lookup

def find_transcripts_for_mxe(mxe_file, gtf_file, output_file='mxe_transcript_mapping.csv'):
    mxe_df = pd.read_csv(mxe_file)
    exon_map = parse_gtf_exons(gtf_file)

    records = []

    for _, row in mxe_df.iterrows():
        chrom = row['chr'].replace('chr', '')
        strand = row['strand']

        exon1_key = (chrom, int(row['1stExonStart_0base']) + 1, int(row['1stExonEnd']), strand)
        exon2_key = (chrom, int(row['2ndExonStart_0base']) + 1, int(row['2ndExonEnd']), strand)

        transcripts1 = exon_map.get(exon1_key, [])
        transcripts2 = exon_map.get(exon2_key, [])

        records.append({
            'chr': row['chr'],
            'strand': strand,
            'Exon1_start': exon1_key[1],
            'Exon1_end': exon1_key[2],
            'Exon2_start': exon2_key[1],
            'Exon2_end': exon2_key[2],
            'transcripts_exon1': ';'.join(transcripts1),
            'transcripts_exon2': ';'.join(transcripts2),
            'geneSymbol': row.get('geneSymbol', ''),
            'FDR': row.get('FDR', ''),
            'IncLevelDifference': row.get('IncLevelDifference', '')
        })

    result_df = pd.DataFrame(records)
    result_df.to_csv(output_file, index=False)
    print(f"Saved transcript mapping to: {output_file}")

# Example usage:
find_transcripts_for_mxe("MXE_filtered.csv", "Macaca_mulatta.Mmul_10.108.gtf")
