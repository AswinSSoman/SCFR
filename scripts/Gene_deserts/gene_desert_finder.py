import pandas as pd
from scipy.stats import zscore
import argparse

def read_bed(file_path):
    # Read a BED file with chromosome, start, end, gene
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'gene'])
    return df

def merge_exons(df):
    # Merge overlapping/adjacent exons on the same chromosome
    merged = []
    for chrom, group in df.groupby('chr'):
        sorted_group = group.sort_values('start')
        current_start = None
        current_end = None
        for _, row in sorted_group.iterrows():
            if current_start is None:
                current_start, current_end = row['start'], row['end']
            elif row['start'] <= current_end:
                current_end = max(current_end, row['end'])
            else:
                merged.append((chrom, current_start, current_end))
                current_start, current_end = row['start'], row['end']
        if current_start is not None:
            merged.append((chrom, current_start, current_end))
    return pd.DataFrame(merged, columns=['chr', 'start', 'end'])

def find_intergenic_regions(merged_df):
    intergenic = []
    for chrom, group in merged_df.groupby('chr'):
        sorted_exons = group.sort_values('start').reset_index(drop=True)
        for i in range(1, len(sorted_exons)):
            prev_end = sorted_exons.loc[i-1, 'end']
            curr_start = sorted_exons.loc[i, 'start']
            if curr_start > prev_end:
                intergenic.append((chrom, prev_end, curr_start,
                                   curr_start - prev_end))
    return pd.DataFrame(intergenic, columns=['chr', 'start', 'end', 'length'])

def identify_deserts(intergenic_df, z_thresh=2):
    intergenic_df['zscore'] = zscore(intergenic_df['length'])
    deserts = intergenic_df[intergenic_df['zscore'] >= z_thresh]
    return deserts

def main(bed_file, z_threshold, output_prefix):
    df = read_bed(bed_file)
    merged_exons = merge_exons(df)
    intergenic = find_intergenic_regions(merged_exons)
    deserts = identify_deserts(intergenic, z_threshold)

    intergenic.to_csv(f"{output_prefix}_intergenic_regions.tsv", sep='\t', index=False)
    deserts.to_csv(f"{output_prefix}_gene_deserts.tsv", sep='\t', index=False)
    print(f"Found {len(deserts)} gene deserts (Z = {z_threshold}).")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify gene deserts based on coding exon BED file.")
    parser.add_argument("bed_file", help="Input BED file with coding exons")
    parser.add_argument("--z", type=float, default=2.0, help="Z-score threshold for gene deserts")
    parser.add_argument("--out", default="output", help="Prefix for output files")
    args = parser.parse_args()
    main(args.bed_file, args.z, args.out)

#python scripts/gene_desert_finder.py genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.bed --z 2 --out gene_deserts
