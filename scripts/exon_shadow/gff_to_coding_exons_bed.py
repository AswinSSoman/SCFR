#!/usr/bin/env python3
import sys
import argparse

def parse_attributes(attr_str):
    """
    Parse GFF3 or GTF-style attributes into a dict.
    Handles both:
      - key=value;key2=value2
      - key "value"; key2 "value2";
    """
    attrs = {}
    if not attr_str:
        return attrs

    # Remove trailing semicolons and split
    fields = [f.strip() for f in attr_str.strip().strip(";").split(";") if f.strip()]
    for f in fields:
        if "=" in f:
            key, value = f.split("=", 1)
        else:
            # GTF style: key "value"
            parts = f.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
        value = value.strip().strip('"')
        attrs[key] = value
    return attrs


def get_gene_name(attrs):
    """
    Try to get a gene name from common attribute keys.
    """
    for key in ["gene_name", "gene", "Name", "gene_id", "ID", "Parent"]:
        if key in attrs:
            # If multiple values (e.g. Parent=a,b), pick the first
            return attrs[key].split(",")[0]
    return "."


def signed_frame(strand, phase):
    """
    Given strand ('+' or '-') and GFF phase ('0','1','2'),
    return signed reading frame in {1,2,3,-1,-2,-3}.

    GFF phase definitions:
      0 -> first base is codon position 1
      1 -> first base is codon position 2
      2 -> first base is codon position 3
    We map that simply as frame = phase + 1, and apply strand sign.
    """
    if phase not in {"0", "1", "2"}:
        return "0"
    frame = int(phase) + 1
    if strand == "-":
        frame = -frame
    return str(frame)


def gff_to_bed(gff_handle, bed_handle):
    for line in gff_handle:
        if not line.strip() or line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue

        seqid, source, ftype, start, end, score, strand, phase, attributes = parts

        # We only want coding exons, so take CDS entries
        if ftype != "CDS":
            continue

        try:
            start_i = int(start)
            end_i = int(end)
        except ValueError:
            continue

        # BED uses 0-based start, end-exclusive
        bed_start = start_i - 1
        bed_end = end_i

        attrs = parse_attributes(attributes)
        gene_name = get_gene_name(attrs)
        sf = signed_frame(strand, phase)

        # Output: chrom, start, end, gene_name, signed_frame
        bed_handle.write(f"{seqid}\t{bed_start}\t{bed_end}\t{gene_name}\t{sf}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert CDS entries from a GFF file to a BED file of coding exons."
    )
    parser.add_argument(
        "-g", "--gff", required=True, help="Input GFF/GTF file"
    )
    parser.add_argument(
        "-o", "--out", default="-",
        help="Output BED file (default: stdout)"
    )
    args = parser.parse_args()

    if args.gff == "-":
        gff_handle = sys.stdin
    else:
        gff_handle = open(args.gff, "r")

    if args.out == "-":
        bed_handle = sys.stdout
    else:
        bed_handle = open(args.out, "w")

    try:
        gff_to_bed(gff_handle, bed_handle)
    finally:
        if gff_handle is not sys.stdin:
            gff_handle.close()
        if bed_handle is not sys.stdout:
            bed_handle.close()


if __name__ == "__main__":
    main()
