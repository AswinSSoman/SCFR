#!/usr/bin/env python3
import sys
import argparse
from collections import defaultdict
import bisect
import os


class Exon:
    __slots__ = ("start", "end", "gene", "frame")

    def __init__(self, start, end, gene, frame):
        self.start = start
        self.end = end
        self.gene = gene
        self.frame = frame


class SCFR:
    __slots__ = ("chrom", "start", "end", "frame")

    def __init__(self, chrom, start, end, frame):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.frame = frame


def parse_exon_file(path):
    """
    Exons file format:
      chrom  start  end  gene  signed_frame
    white space or tab separated, 0 based, end exclusive.
    """
    exons_by_key = defaultdict(list)
    starts_by_key = defaultdict(list)

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            chrom, s, e, gene, frame_str = fields[:5]
            try:
                start = int(s)
                end = int(e)
                frame = int(frame_str)
            except ValueError:
                continue
            key = (chrom, frame)
            ex = Exon(start, end, gene, frame)
            exons_by_key[key].append(ex)

    # sort exons by start coordinate and build starts arrays for bisect
    for key in list(exons_by_key.keys()):
        exons = exons_by_key[key]
        exons.sort(key=lambda x: x.start)
        starts = [ex.start for ex in exons]
        starts_by_key[key] = starts

    return exons_by_key, starts_by_key


def scfr_stream(path):
    """
    Generator over SCFR entries.
    SCFR file format:
      chrom  start  end  signed_frame
    """
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 4:
                continue
            chrom, s, e, frame_str = fields[:4]
            try:
                start = int(s)
                end = int(e)
                frame = int(frame_str)
            except ValueError:
                continue
            yield SCFR(chrom, start, end, frame)


def find_exons_in_scfr(scfr, exons_by_key, starts_by_key):
    """
    Find exons that are fully contained in the SCFR, with same chrom and signed frame.
    """
    key = (scfr.chrom, scfr.frame)
    if key not in exons_by_key:
        return []

    exons = exons_by_key[key]
    starts = starts_by_key[key]

    # first exon whose start >= scfr.start
    idx = bisect.bisect_left(starts, scfr.start)
    inside = []

    # iterate until exon.start >= scfr.end
    for i in range(idx, len(exons)):
        ex = exons[i]
        if ex.start >= scfr.end:
            break
        # require whole exon inside SCFR
        if ex.start >= scfr.start and ex.end <= scfr.end:
            inside.append(ex)

    return inside


def compute_up_down_lengths(scfr, exons_in_scfr):
    """
    For a list of exons inside one SCFR (all same signed frame, sorted by start),
    compute upstream and downstream non exonic lengths inside the SCFR for each exon.

    For plus strand (frame > 0):
      upstream = from SCFR start or previous exon end to exon start
      downstream = from exon end to next exon start or SCFR end

    For minus strand (frame < 0):
      upstream (toward transcription start) is on genomic right side:
        upstream = from exon end to next exon start or SCFR end
      downstream = from SCFR start or previous exon end to exon start
    """
    n = len(exons_in_scfr)
    if n == 0:
        return []

    strand_positive = exons_in_scfr[0].frame > 0

    results = []
    for i, ex in enumerate(exons_in_scfr):
        left_boundary = scfr.start if i == 0 else exons_in_scfr[i - 1].end
        right_boundary = scfr.end if i == n - 1 else exons_in_scfr[i + 1].start

        if strand_positive:
            up_len = ex.start - left_boundary
            down_len = right_boundary - ex.end
        else:
            up_len = right_boundary - ex.end
            down_len = ex.start - left_boundary

        if up_len < 0:
            up_len = 0
        if down_len < 0:
            down_len = 0

        results.append((up_len, down_len))

    return results


def identify_exitron_candidates(scfr, exons_in_scfr):
    """
    For a SCFR with >= 2 exons, identify pairs of exons that can produce exitrons.

    Criteria:
      * same gene
      * same signed frame (true by construction)
      * gap between them is positive and multiple of 3

    For plus strand:
      gap = next.start - current.end

    For minus strand:
      exons are still sorted by genomic start
      transcription order is reversed
      consider adjacent exons in genomic order:
        gap = current.start - next.end
    """
    candidates = []
    n = len(exons_in_scfr)
    if n < 2:
        return candidates

    strand_positive = exons_in_scfr[0].frame > 0

    for i in range(n - 1):
        ex1 = exons_in_scfr[i]
        ex2 = exons_in_scfr[i + 1]

        if ex1.gene != ex2.gene:
            continue

        if strand_positive:
            gap_len = ex2.start - ex1.end
            gap_start = ex1.end
            gap_end = ex2.start
        else:
            gap_len = ex2.start - ex1.end
            gap_start = ex1.end
            gap_end = ex2.start

        if gap_len <= 0:
            continue

        multiple_of_3 = (gap_len % 3 == 0)
        if multiple_of_3:
            candidates.append(
                {
                    "gene": ex1.gene,
                    "ex1_start": ex1.start,
                    "ex1_end": ex1.end,
                    "ex2_start": ex2.start,
                    "ex2_end": ex2.end,
                    "gap_start": gap_start,
                    "gap_end": gap_end,
                    "gap_len": gap_len,
                }
            )

    return candidates


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Match exons to stop codon free regions (SCFRs) and report "
            "single exon SCFRs, multi exon SCFRs and exitron candidates."
        )
    )
    parser.add_argument("exons", help="Exons BED file (chrom start end gene signed_frame)")
    parser.add_argument("scfrs", help="SCFR BED file (chrom start end signed_frame)")
    parser.add_argument(
        "-p",
        "--prefix",
        help="Output prefix (default: derived from SCFR filename)",
    )
    args = parser.parse_args()

    exon_file = args.exons
    scfr_file = args.scfrs

    if args.prefix:
        prefix = args.prefix
    else:
        base = os.path.basename(scfr_file)
        prefix = os.path.splitext(base)[0] or "scfr_exon"

    out_single = prefix + ".single_exon.txt"
    out_multi = prefix + ".multi_exon.txt"
    out_exitron = prefix + ".exitron_candidates.txt"

    # read and index exons
    exons_by_key, starts_by_key = parse_exon_file(exon_file)

    # open outputs
    f_single = open(out_single, "w")
    f_multi = open(out_multi, "w")
    f_exit = open(out_exitron, "w")

    # headers
    f_single.write(
        "#chrom\tSCFR_start\tSCFR_end\tSCFR_frame\t"
        "exon_start\texon_end\tgene\texon_frame\t"
        "upstream_len_in_SCFR\tdownstream_len_in_SCFR\n"
    )

    f_multi.write(
        "#chrom\tSCFR_start\tSCFR_end\tSCFR_frame\t"
        "exon_start\texon_end\tgene\texon_frame\t"
        "upstream_len_in_SCFR\tdownstream_len_in_SCFR\t"
        "num_exons_in_SCFR\tgroup_id\n"
    )

    f_exit.write(
        "#chrom\tSCFR_start\tSCFR_end\tSCFR_frame\tgene\t"
        "exon1_start\texon1_end\texon2_start\texon2_end\t"
        "gap_start\tgap_end\tgap_len\n"
    )

    # summary counters
    total_scfr = 0
    scfr_zero = 0
    scfr_one = 0
    scfr_multi = 0
    exitron_scfr_count = 0
    exitron_pair_count = 0

    # group id for multi exon SCFRs
    multi_group_id = 0

    # process SCFRs streaming
    for scfr in scfr_stream(scfr_file):
        total_scfr += 1
        exons_in = find_exons_in_scfr(scfr, exons_by_key, starts_by_key)
        exons_in.sort(key=lambda x: x.start)
        m = len(exons_in)

        if m == 0:
            scfr_zero += 1
            continue

        # compute upstream and downstream lengths
        up_down = compute_up_down_lengths(scfr, exons_in)

        if m == 1:
            scfr_one += 1
            ex = exons_in[0]
            up_len, down_len = up_down[0]
            f_single.write(
                f"{scfr.chrom}\t{scfr.start}\t{scfr.end}\t{scfr.frame}\t"
                f"{ex.start}\t{ex.end}\t{ex.gene}\t{ex.frame}\t"
                f"{up_len}\t{down_len}\n"
            )
        else:
            scfr_multi += 1
            multi_group_id += 1  # new group id for this multi exon SCFR

            for (ex, (up_len, down_len)) in zip(exons_in, up_down):
                f_multi.write(
                    f"{scfr.chrom}\t{scfr.start}\t{scfr.end}\t{scfr.frame}\t"
                    f"{ex.start}\t{ex.end}\t{ex.gene}\t{ex.frame}\t"
                    f"{up_len}\t{down_len}\t{m}\t{multi_group_id}\n"
                )

            # exitron candidates
            candidates = identify_exitron_candidates(scfr, exons_in)
            if candidates:
                exitron_scfr_count += 1
                for c in candidates:
                    exitron_pair_count += 1
                    f_exit.write(
                        f"{scfr.chrom}\t{scfr.start}\t{scfr.end}\t{scfr.frame}\t"
                        f"{c['gene']}\t"
                            f"{c['ex1_start']}\t{c['ex1_end']}\t"
                            f"{c['ex2_start']}\t{c['ex2_end']}\t"
                            f"{c['gap_start']}\t{c['gap_end']}\t{c['gap_len']}\n"
                    )

    f_single.close()
    f_multi.close()
    f_exit.close()

    # summary
    print("Summary")
    print("=======")
    print(f"Total SCFRs processed: {total_scfr}")
    print(f"SCFRs with 0 whole exons: {scfr_zero}")
    print(f"SCFRs with exactly 1 whole exon: {scfr_one}")
    print(f"SCFRs with >1 whole exon: {scfr_multi}")
    print(f"SCFRs with at least one exitron candidate pair: {exitron_scfr_count}")
    print(f"Total exitron candidate exon pairs: {exitron_pair_count}")
    print()
    print("Output files:")
    print(f"  Single exon SCFRs: {out_single}")
    print(f"  Multi exon SCFRs:  {out_multi}")
    print(f"  Exitrons:          {out_exitron}")


if __name__ == "__main__":
    main()
