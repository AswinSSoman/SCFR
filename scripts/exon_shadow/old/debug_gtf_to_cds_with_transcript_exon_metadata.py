#!/usr/bin/env python3

import argparse
import glob
from collections import defaultdict

def dprint(debug, *args):
    if debug:
        print("DEBUG:", *args)

def parse_attributes(attr_str, debug=False):
    attrs = {}
    for field in attr_str.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
        else:
            parts = field.split(None, 1)
            if len(parts) != 2:
                continue
            k, v = parts
        attrs[k] = v.strip().strip('"')

    dprint(debug, "Parsed attributes:", attrs)
    return attrs

def get_gene_name(attrs):
    for k in ["gene_name", "gene", "Name", "gene_id", "ID", "Parent"]:
        if k in attrs:
            return attrs[k].split(",")[0]
    return "NA"

def get_chrom_sizes(species_pattern, debug=False):
    chrom_sizes = {}
    base_path = "/media/aswin/SCFR/SCFR-main/genome_reports/*.tsv"

    dprint(debug, "Genome report glob path:", base_path)

    report_files = glob.glob(base_path)
    dprint(debug, "Genome report files found:", report_files)

    target_files = [f for f in report_files if species_pattern.lower() in f.lower()]
    dprint(debug, "Filtered genome report files:", target_files)

    if not target_files:
        print(f"ERROR: No genome report found for species '{species_pattern}'")
        return chrom_sizes

    target_file = target_files[0]
    dprint(debug, "Using genome report:", target_file)

    with open(target_file) as f:
        for line in f:
            if line.startswith("#") or "Seq length" in line:
                continue
            parts = line.rstrip().split("\t")
            if len(parts) >= 11:
                acc = parts[8]
                try:
                    size = int(parts[10])
                    chrom_sizes[acc] = size
                except ValueError:
                    continue

    dprint(debug, "Loaded chromosome sizes (sample):",
           list(chrom_sizes.items())[:5])
    return chrom_sizes

def calculate_genomic_frame(start, end, strand, phase, chrom_size, debug=False):
    try:
        dprint(debug, "Frame input:",
               f"start={start}",
               f"end={end}",
               f"strand={strand}",
               f"phase={phase}",
               f"chrom_size={chrom_size}")

        if phase == ".":
            dprint(debug, "Phase is '.', returning NA")
            return "NA"

        p = int(phase)

        if strand == "+":
            frame = ((start - (start - (start % 3)) + p - 1) % 3) + 1
            dprint(debug, "Plus strand frame =", frame)
            return str(frame)

        elif strand == "-":
            val = chrom_size - end + 1
            mod = val % 3
            frame = abs(val - val - mod - p)
            final_frame = frame if frame != 0 else 3

            dprint(debug, "Minus strand calc:",
                   f"val={val}",
                   f"val%3={mod}",
                   f"p={p}",
                   f"raw_frame={frame}",
                   f"final_frame={final_frame}")

            return str(final_frame)

    except Exception as e:
        dprint(debug, "Frame calculation error:", e)
        return "NA"

    return "NA"

def main():
    parser = argparse.ArgumentParser(description="Convert GTF to BED with Genomic Frame")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-s", "--species", default="human")
    parser.add_argument("--debug", action="store_true",
                        help="Print intermediate values")
    parser.add_argument("--debug-tx",
                        help="Debug only this transcript_id")
    parser.add_argument("--debug-chr",
                        help="Debug only this chromosome")

    args = parser.parse_args()

    chrom_sizes = get_chrom_sizes(args.species, args.debug)

    cds_by_tx = defaultdict(list)
    gene_to_txs = defaultdict(set)
    exon_to_txs = defaultdict(set)

    # PASS 1
    with open(args.input) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            chrom, _, feature, start, end, _, strand, phase, attrs_str = fields
            if feature != "CDS":
                continue

            if args.debug_chr and chrom != args.debug_chr:
                continue

            attrs = parse_attributes(attrs_str, args.debug)
            tx = attrs.get("transcript_id", "NA")

            if args.debug_tx and tx != args.debug_tx:
                continue

            dprint(args.debug, "Processing CDS:",
                   chrom, start, end, strand, phase, tx)

            cs = chrom_sizes.get(chrom, 0)
            dprint(args.debug, "Chrom size lookup:", chrom, cs)

            frame = calculate_genomic_frame(
                int(start), int(end), strand, phase, cs, args.debug
            )

            gene = get_gene_name(attrs)

            cds_by_tx[tx].append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "gene": gene,
                "tx": tx,
                "frame": frame,
                "strand": strand
            })

            gene_to_txs[gene].add(tx)
            exon_key = (chrom, int(start), int(end), gene, frame)
            exon_to_txs[exon_key].add(tx)

            dprint(args.debug, "Stored CDS dict:", cds_by_tx[tx][-1])

    # PASS 2
    with open(args.output, "w") as out:
        for tx, exons in cds_by_tx.items():
            strand = exons[0]["strand"]
            exons.sort(
                key=lambda x: x["start"],
                reverse=(strand == "-")
            )

            total_exons = len(exons)

            for exon_number, e in enumerate(exons, start=1):
                gene = e["gene"]
                exon_key = (e["chrom"], e["start"], e["end"], gene, e["frame"])

                gene_tx_count = len(gene_to_txs[gene])
                exon_tx_count = len(exon_to_txs[exon_key])
                sharing = exon_tx_count / gene_tx_count

                exon_order = (
                    "first" if exon_number == 1 else
                    "last" if exon_number == total_exons else
                    "middle"
                )

                splicing = (
                    "constitutive" if sharing == 1.0 else
                    "unique" if exon_tx_count == 1 else
                    "alternative"
                )

                bed_line = (
                    f"{e['chrom']}\t{e['start']}\t{e['end']}\t"
                    f"{gene}\t{e['tx']}\t{e['strand']}\t{e['frame']}\t"
                    f"{exon_number}\t{gene_tx_count}\t{exon_tx_count}\t"
                    f"{exon_order}\t{sharing:.3f}\t{splicing}\n"
                )

                dprint(args.debug, "Writing BED line:", bed_line.strip())
                out.write(bed_line)

if __name__ == "__main__":
    main()
