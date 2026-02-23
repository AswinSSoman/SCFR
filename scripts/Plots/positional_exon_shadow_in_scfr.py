import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def main():
    # 1. Setup Command Line Arguments
    parser = argparse.ArgumentParser(description='Convert SCFR length into bins and count feature distribution.')
    parser.add_argument('-i', '--input', required=True, help='Path to input text file (tab-separated)')
    parser.add_argument('-t', '--table', default='bin_counts.csv', help='Path for output CSV table')
    parser.add_argument('-p', '--pdf', default='scfr_distribution.pdf', help='Path for output PDF plot')
    parser.add_argument('-b', '--bins', type=int, default=20, help='Number of equal bins (default: 20)')
    
    args = parser.parse_args()

    # 2. Load the data
    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    # 3. Calculate SCFR total length and relative positions
    df['SCFR_total_len'] = df['SCFR_end'] - df['SCFR_start']
    
    # Position of Upstream end (normalized 0 to 1)
    df['upstream_rel'] = df['upstream_len_in_SCFR'] / df['SCFR_total_len']
    
    # Position of Downstream start (normalized 0 to 1)
    df['downstream_rel'] = (df['SCFR_total_len'] - df['downstream_len_in_SCFR']) / df['SCFR_total_len']

    # 4. Define Bins and Labels
    bins = np.linspace(0, 1, args.bins + 1)
    bin_ids = range(1, args.bins + 1)
    
    # Create percentage labels like (0,5], (5,10], etc.
    tick_labels = []
    for i in range(len(bins)-1):
        start_pct = int(bins[i] * 100)
        end_pct = int(bins[i+1] * 100)
        tick_labels.append(f"({start_pct},{end_pct}]")

    # 5. Categorize and Count
    df['up_bin'] = pd.cut(df['upstream_rel'], bins=bins, labels=bin_ids, include_lowest=True)
    df['down_bin'] = pd.cut(df['downstream_rel'], bins=bins, labels=bin_ids, include_lowest=True)

    up_counts = df['up_bin'].value_counts().reindex(bin_ids).fillna(0).astype(int)
    down_counts = df['down_bin'].value_counts().reindex(bin_ids).fillna(0).astype(int)

    # 6. Create and Save Summary Table
    bin_summary = pd.DataFrame({
        'bin_id': bin_ids,
        'bin_label': tick_labels,
        'upstream_count': up_counts.values,
        'downstream_count': down_counts.values
    }).set_index('bin_id')
    
    bin_summary.to_csv(args.table)
    print(f"Summary table saved to: {args.table}")

    # 7. Visualization
    plt.figure(figsize=(12, 7))
    x_indices = np.arange(len(tick_labels))
    bar_width = 0.35

    # Bar Plot
    plt.bar(x_indices - bar_width/2, bin_summary['upstream_count'], bar_width, 
            label='Upstream Counts', color='skyblue', alpha=0.7)
    plt.bar(x_indices + bar_width/2, bin_summary['downstream_count'], bar_width, 
            label='Downstream Counts', color='salmon', alpha=0.7)

    # Line Plot at the top
    plt.plot(x_indices, bin_summary['upstream_count'], marker='o', color='blue', 
             label='Upstream Trend', linewidth=2)
    plt.plot(x_indices, bin_summary['downstream_count'], marker='s', color='red', 
             label='Downstream Trend', linewidth=2)

    # Labels and Formatting
    plt.xlabel('Normalized position in SCFR', fontsize=12)
    plt.ylabel('Frequency of Exon Shadow', fontsize=12)
    plt.title(f'Distribution across {args.bins} SCFR Bins', fontsize=14)
    plt.xticks(x_indices, tick_labels, rotation=45, ha='right')
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(args.pdf)
    print(f"Plot saved to: {args.pdf}")

if __name__ == "__main__":
    main()
