#!/usr/bin/env python3
import pandas as pd
import sys
import os
import numpy as np

# Helper function to format numbers: max 1 decimal, no .0 for integers
def format_stat(value):
    """Formats a number to a maximum of 1 decimal point, suppressing .0."""
    if pd.isna(value):
        return None
    
    # Round to the nearest 1 decimal place
    rounded_value = round(value, 1)
    
    # Check if the rounded value is an integer (i.e., has no decimal part)
    if rounded_value == int(rounded_value):
        return int(rounded_value)
    else:
        # Format to one decimal place if it has a decimal part
        return f"{rounded_value:.1f}"

def summarize_lengths(lengths):
    """Calculates descriptive statistics for a series of lengths."""
    # Find the mode
    mode_val = lengths.mode()
    # Get the first mode value, or None if the series is empty
    mode_val = mode_val.iloc[0] if not mode_val.empty else np.nan

    return {
        "N": len(lengths),
        "min": lengths.min(),
        "q1": lengths.quantile(0.25),
        "median": lengths.median(),
        "mean": lengths.mean(),
        "mode": mode_val,
        "q3": lengths.quantile(0.75),
        "max": lengths.max(),
        "sd": lengths.std()
    }

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 gene_desert_stats.py gene_deserts/*.bed")
        sys.exit(1)

    files = sys.argv[1:]

    all_rows = []
    all_lengths = []

    for f in files:
        species = os.path.basename(f).replace(".gene_deserts.bed", "")

        try:
            # Added error handling for missing file
            df = pd.read_csv(f, sep="\t", header=None, names=["chr", "start", "end"])
        except FileNotFoundError:
            print(f"Error: File not found: {f}", file=sys.stderr)
            continue
        except pd.errors.EmptyDataError:
            print(f"Warning: File is empty or malformed: {f}. Skipping.", file=sys.stderr)
            continue

        df["length"] = df["end"] - df["start"]

        if df["length"].empty:
            print(f"Warning: No lengths calculated for {species}. Skipping.", file=sys.stderr)
            continue

        stats_dict = summarize_lengths(df["length"])
        stats_dict["species"] = species
        all_rows.append(stats_dict)

        temp = pd.DataFrame({"species": species, "length": df["length"]})
        all_lengths.append(temp)

    # --- Summary Table Generation and Formatting ---
    if not all_rows:
        print("Error: No data processed. Exiting.")
        sys.exit(1)

    summary_df = pd.DataFrame(all_rows)
    summary_df = summary_df[
        ["species", "N", "min", "q1", "median", "mean", "mode", "q3", "max", "sd"]
    ].copy() # Use .copy() to avoid SettingWithCopyWarning
    
    # Define columns to format (all numerical stats except 'N')
    stat_cols = ["min", "q1", "median", "mean", "mode", "q3", "max", "sd"]

    # Apply the custom formatting function to the numerical columns
    # We use applymap for element-wise application
    summary_df[stat_cols] = summary_df[stat_cols].applymap(format_stat)
    
    # Ensure 'N' is integer and 'species' is string (just in case)
    summary_df['N'] = summary_df['N'].astype(int)
    
    # Save summary table
    summary_df.to_csv("desert_summary.tsv", sep="\t", index=False)

    # --- Lengths Table Generation ---
    lengths_df = pd.concat(all_lengths, ignore_index=True)
    lengths_df.to_csv("all_desert_lengths.tsv", sep="\t", index=False)

    print("✔ desert_summary.tsv generated")
    print("✔ all_desert_lengths.tsv generated")

if __name__ == "__main__":
    main()
