#!/usr/bin/env python3
"""
Stream a BED file of SCFRs and compute sliding window strand asymmetries.
Writes output incrementally to save memory and improve speed.
"""

import sys
import argparse
from pathlib import Path
from collections import defaultdict

FORWARD = [1, 2, 3]
REVERSE = [-1, -2, -3]


def asym(a, b):
    return (a - b) / (a + b) if (a + b) else 0.0


def get_windows(start, end, window_size, slide_size):
    """
    Return all window start positions that overlap a given interval.
    """
    windows = []
    win_start = (start // slide_size) * slide_size
    while win_start < end:
        win_end = win_start + window_size
        if win_end > start:
            windows.append(win_start)
        win_start += slide_size
    return windows


def flush_completed_windows(chrom, active, current_pos, slide_size, window_size, out_fh):
    """
    Write completed windows (those that ended before current_pos) to file.
    """
    to_delete = []
    for win_start in active[chrom]:
        if win_start + window_size <= current_pos:
            data = active[chrom][win_start]
            f_cnt = sum(data[f][0] for f in FORWARD)
            r_cnt = sum(data[f][0] for f in REVERSE)
            f_len = sum(data[f][1] for f in FORWARD)
            r_len = sum(data[f][1] for f in REVERSE)
            out_fh.write(
                f"{chrom},{win_start},{win_start + window_size},"
                f"{f_cnt},{r_cnt},{asym(f_cnt, r_cnt):.6f},"
                f"{f_len},{r_len},{asym(f_len, r_len):.6f}\n"
            )
            to_delete.append(win_start)

    for win_start in to_delete:
        del active[chrom][win_start]


def main():
    parser = argparse.ArgumentParser(description="Stream SCFR BED file and compute strand asymmetries in sliding windows.")
    parser.add_argument("bed_file", help="Input BED file")
    parser.add_argument("--window-size", type=int, required=True, help="Window size (bp)")
    parser.add_argument("--slide-size", type=int, required=True, help="Slide size (bp)")
    parser.add_argument("--output", type=str, default="scfr_window_metrics.csv", help="Output CSV file")

    args = parser.parse_args()

    active_windows = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0, 0])))

    with open(args.bed_file) as infile, open(args.output, "w") as out_fh:
        out_fh.write("chrom,window_start,window_end,f_count,r_count,strand_count_asym,f_length,r_length,strand_length_asym\n")

        for line in infile:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, frame = line.strip().split("\t")[:4]
            start, end, frame = int(start), int(end), int(frame)
            length = end - start

            # Flush completed windows
            flush_completed_windows(chrom, active_windows, start, args.slide_size, args.window_size, out_fh)

            # Update overlapping windows
            for win_start in get_windows(start, end, args.window_size, args.slide_size):
                active_windows[chrom][win_start][frame][0] += 1
                active_windows[chrom][win_start][frame][1] += length

        # Flush any remaining windows
        for chrom in active_windows:
            flush_completed_windows(chrom, active_windows, float('inf'), args.slide_size, args.window_size, out_fh)

    print(f"\nDone. Output written to: {args.output}\n")


if __name__ == "__main__":
    main()
