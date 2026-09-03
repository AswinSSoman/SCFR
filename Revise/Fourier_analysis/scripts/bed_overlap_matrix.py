#!/usr/bin/env python3
import sys
import os
import subprocess
import glob
import re

def count_lines(filepath):
    """Counts non-empty, non-header lines in a BED file."""
    count = 0
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith(('#', 'track', 'browser')):
                count += 1
    return count

def compute_total_bp(filepath):
    """Calculates total base pairs in a BED file (sum of end - start)."""
    total_bp = 0
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith(('#', 'track', 'browser')):
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    try:
                        start = int(parts[1])
                        end = int(parts[2])
                        total_bp += (end - start)
                    except ValueError:
                        continue
    return total_bp

def run_bedtools_intersections(bed_a, bed_b):
    """
    Runs bedtools intersect to get:
    1. Feature overlap count (-u)
    2. Sum of overlapping base pairs (-wo)
    """
    # 1. Feature Overlap Count
    cmd_feat = ["bedtools", "intersect", "-u", "-a", bed_a, "-b", bed_b]
    feat_count = 0
    try:
        res_feat = subprocess.run(cmd_feat, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        feat_count = len([line for line in res_feat.stdout.split('\n') if line.strip()])
    except (FileNotFoundError, subprocess.CalledProcessError):
        feat_count = 0

    # 2. Base Pair Overlap Sum (-wo outputs overlap length as last column)
    # Pipeline: bedtools intersect -wo -a A -b B | awk '{s+=$NF} END {print s+0}'
    cmd_bp = f"bedtools intersect -wo -a {bed_a} -b {bed_b} | awk '{{s+=$NF}} END {{print s+0}}'"
    bp_overlap = 0
    try:
        res_bp = subprocess.run(cmd_bp, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        out = res_bp.stdout.strip()
        bp_overlap = int(out) if out else 0
    except (subprocess.CalledProcessError, ValueError):
        bp_overlap = 0

    return feat_count, bp_overlap

def color_cell(pct):
    """Returns ANSI colored percentage string."""
    if pct == 100.0:
        bg_color = "\033[41;30m"    # Red BG
    elif pct >= 75.0:
        bg_color = "\033[43;30m"    # Yellow BG
    elif pct >= 50.0:
        bg_color = "\033[103;30m"   # Light Yellow BG
    elif pct >= 25.0:
        bg_color = "\033[44;37m"    # Blue BG
    elif pct > 0.0:
        bg_color = "\033[100;37m"   # Dark Gray BG
    else:
        bg_color = "\033[90m"       # Dim Gray
    
    reset = "\033[0m"
    return f"{bg_color}{pct:5.1f}%{reset}"

def pad_cell(content, target_width):
    """Pads string to target_width while ignoring invisible ANSI escape characters."""
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    visible_length = len(ansi_escape.sub('', content))
    needed_padding = max(0, target_width - visible_length)
    
    left_pad = needed_padding // 2
    right_pad = needed_padding - left_pad
    
    return " " * left_pad + content + " " * right_pad

def print_matrix(labels, matrix, col_widths, is_percentage=False, title=""):
    """Renders structured tables with dynamic column widths aligned to borders."""
    n = len(labels)
    max_label_len = max(len(l) for l in labels)
    label_margin = max(max_label_len + 6, 40)
    
    print(f"\033[1;36m{title}\033[0m")
    
    # Column Number Headers
    col_hdr = " " * (label_margin + 2)
    for j in range(n):
        hdr_str = f"[{j+1}]"
        col_hdr += pad_cell(hdr_str, col_widths[j]) + " "
    print(col_hdr)

    # Top Border
    top_border = " " * label_margin + "┌" + "┬".join("─" * w for w in col_widths) + "┐"
    print(top_border)

    # Rows
    for i in range(n):
        row_title = f"[{i+1:>2}] {labels[i]}"
        row_str = f"{row_title:<{label_margin}} │"
        
        for j in range(n):
            if is_percentage:
                cell_content = color_cell(matrix[i][j])
            else:
                # Format numbers with commas for readability if integer
                val = matrix[i][j]
                cell_content = f"{val:,}" if isinstance(val, int) else str(val)
                
            row_str += pad_cell(cell_content, col_widths[j]) + "│"
            
        print(row_str)
        
        if i < n - 1:
            mid_border = " " * label_margin + "├" + "┼".join("─" * w for w in col_widths) + "┤"
            print(mid_border)

    # Bottom Border
    bot_border = " " * label_margin + "└" + "┴".join("─" * w for w in col_widths) + "┘\n"
    print(bot_border)

def main():
    if len(sys.argv) > 1:
        bed_files = sys.argv[1:]
    else:
        bed_files = sorted(glob.glob("*.bed"))

    bed_files = [f for f in bed_files if os.path.isfile(f)]
    if len(bed_files) < 2:
        print("Error: Provide at least 2 valid BED files.", file=sys.stderr)
        sys.exit(1)

    labels = [os.path.basename(f).replace(".bed", "") for f in bed_files]
    n = len(bed_files)

    print("\nCounting features and calculating base pair sizes...")
    feat_counts = [count_lines(f) for f in bed_files]
    total_bps = [compute_total_bp(f) for f in bed_files]

    count_matrix = [[0] * n for _ in range(n)]
    pct_matrix = [[0.0] * n for _ in range(n)]
    bp_matrix = [[0] * n for _ in range(n)]

    print("Computing directional pairwise overlaps and base pair intersections...\n")
    for i in range(n):
        for j in range(n):
            if i == j:
                count_matrix[i][j] = feat_counts[i]
                pct_matrix[i][j] = 100.0 if feat_counts[i] > 0 else 0.0
                bp_matrix[i][j] = total_bps[i]
            else:
                feat_overlap, bp_overlap = run_bedtools_intersections(bed_files[i], bed_files[j])
                count_matrix[i][j] = feat_overlap
                pct_matrix[i][j] = (feat_overlap / feat_counts[i] * 100.0) if feat_counts[i] > 0 else 0.0
                bp_matrix[i][j] = bp_overlap

    # Summary Stats
    print("=" * 85)
    print(" SUMMARY STATISTICS")
    print("=" * 85)
    for idx, (label, cnt, bp) in enumerate(zip(labels, feat_counts, total_bps), 1):
        print(f"  [{idx}] {label:<45} : {cnt:>6} features | {bp:>12,} total bp")
    print("=" * 85 + "\n")

    # Column Width Calculations
    pct_col_widths = [8] * n

    raw_col_widths = []
    for j in range(n):
        max_val = max(count_matrix[i][j] for i in range(n))
        raw_col_widths.append(max(7, len(f"{max_val:,}") + 4))

    bp_col_widths = []
    for j in range(n):
        max_val = max(bp_matrix[i][j] for i in range(n))
        bp_col_widths.append(max(9, len(f"{max_val:,}") + 4))

    # 1. Print Percentage Heatmap
    print_matrix(
        labels, 
        pct_matrix, 
        pct_col_widths, 
        is_percentage=True, 
        title="1. DIRECTIONAL % FEATURE OVERLAP MATRIX (% of ROW present in COLUMN)"
    )

    # 2. Print Raw Feature Counts Heatmap
    print_matrix(
        labels, 
        count_matrix, 
        raw_col_widths, 
        is_percentage=False, 
        title="2. RAW FEATURE OVERLAP COUNT MATRIX"
    )

    # 3. Print Total Base Pair Overlap Heatmap
    print_matrix(
        labels, 
        bp_matrix, 
        bp_col_widths, 
        is_percentage=False, 
        title="3. TOTAL BASE PAIR (bp) OVERLAP MATRIX"
    )

if __name__ == "__main__":
    main()
