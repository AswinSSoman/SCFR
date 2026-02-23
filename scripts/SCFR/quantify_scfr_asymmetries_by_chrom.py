#!/usr/bin/env python3
"""
Parse a BED file of Stop-Codon-Free Regions (SCFRs) whose 4th column
contains the reading frame (1,2,3,-1,-2,-3).

For every chromosome (and for ALL chromosomes together) the script outputs:
    Forward or reverse SCFR counts and total lengths
    Strand-level asymmetry metrics
    Counts per individual reading frame

The same statistics are written to CSV so they can be plotted in R.
"""

import pandas as pd
import sys
from collections import defaultdict
from pathlib import Path

FORWARD = [1, 2, 3]
REVERSE = [-1, -2, -3]

# ----------------------------------------------------------------------
def parse_bed(path):
    """
    Returns
    -------
    chrom_data : dict
        {chrom : {frame : [count, total_len]}}
    """
    chrom_data = defaultdict(lambda: defaultdict(lambda: [0, 0]))

    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, frame = line.rstrip().split("\t")[:4]
            frame = int(frame)
            length = int(end) - int(start)

            chrom_data[chrom][frame][0] += 1
            chrom_data[chrom][frame][1] += length

    return chrom_data


# ----------------------------------------------------------------------
def asym(a, b):
    return (a - b) / (a + b) if (a + b) else 0.0


def summarise(chrom, frame_dict):
    """
    Build a dict with all metrics for one chromosome.
    """
    f_cnt = sum(frame_dict[f][0] for f in FORWARD)
    r_cnt = sum(frame_dict[f][0] for f in REVERSE)
    f_len = sum(frame_dict[f][1] for f in FORWARD)
    r_len = sum(frame_dict[f][1] for f in REVERSE)

    row = {
        "chrom": chrom,
        "f_count": f_cnt,
        "r_count": r_cnt,
        "strand_count_asym": asym(f_cnt, r_cnt),
        "f_length": f_len,
        "r_length": r_len,
        "strand_length_asym": asym(f_len, r_len),
    }

    # individual frame counts (missing frames default to 0)
    for f in FORWARD + REVERSE:
        row[f"frame{f}_count"] = frame_dict.get(f, [0, 0])[0]
        row[f"frame{f}_length"] = frame_dict.get(f, [0, 0])[1]

    return row


# ----------------------------------------------------------------------
def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        sys.exit(
            "Usage: python quantify_scfr_asymmetries_by_chrom.py <SCFR.bed> [output.csv]"
        )

    bed_file = Path(sys.argv[1])
    out_csv = Path(sys.argv[2]) if len(sys.argv) == 3 else Path("scfr_metrics_per_chrom.csv")

    chrom_data = parse_bed(bed_file)

    # per-chromosome rows
    rows = [summarise(chrom, data) for chrom, data in chrom_data.items()]

    # combined ALL row
    total_frame = defaultdict(lambda: [0, 0])
    for data in chrom_data.values():
        for fr, vals in data.items():
            total_frame[fr][0] += vals[0]
            total_frame[fr][1] += vals[1]
    rows.append(summarise("ALL", total_frame))

    df = pd.DataFrame(rows).sort_values("chrom")
    df.to_csv(out_csv, index=False)

    # ------------------- console summary -------------------
    for _, r in df.iterrows():
        print(f"\n=== {r.chrom} ===")
        print(f"Forward SCFR count : {r.f_count}")
        print(f"Reverse SCFR count : {r.r_count}")
        print(f"Strand count asym  : {r.strand_count_asym:.6f}")
        print(f"Forward SCFR length: {r.f_length}")
        print(f"Reverse SCFR length: {r.r_length}")
        print(f"Strand length asym : {r.strand_length_asym:.6f}")
        print(
            "Forward frame counts:",
            {f: int(r[f"frame{f}_count"]) for f in FORWARD},
        )
        print(
            "Reverse frame counts:",
            {f: int(r[f"frame{f}_count"]) for f in REVERSE},
        )

    print(f"\nMetrics written to: {out_csv}\n")


if __name__ == "__main__":
    main()
