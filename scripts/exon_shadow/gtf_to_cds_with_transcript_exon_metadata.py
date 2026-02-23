#!/usr/bin/env python3

import argparse
import glob
import os
from collections import defaultdict

def parse_attributes(attr_str):
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
    return attrs

def get_gene_name(attrs):
    for k in ["gene_name", "gene", "Name", "gene_id", "ID", "Parent"]:
        if k in attrs:
            return attrs[k].split(",")[0]
    return "NA"
                                                              
def get_chrom_sizes(species_pattern):
    chrom_sizes = {}
    base_path = "/media/aswin/SCFR/SCFR-main/genome_reports/*.tsv"
    
    try:
        report_files = glob.glob(base_path)
        target_files = [f for f in report_files if species_pattern.lower() in f.lower()]
        if not target_files:
            print(f"Error: No report file found for {species_pattern}")
            return chrom_sizes
            
        target_file = target_files[0]
        with open(target_file, 'r') as f:
            for line in f:
                if line.startswith("#") or "Seq length" in line: 
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 11:
                    acc = parts[8]
                    try:
                        size = int(parts[10])
                        chrom_sizes[acc] = size
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Warning: Could not load genome report: {e}")
    return chrom_sizes

def calculate_genomic_frame(start, end, strand, phase, chrom_size):
    try:
        if phase == ".": return "NA"
        p = int(phase)
        if strand == "+":
            # awk: (($4-($4-($4%3))+$8-1)%3)+1
            frame = ((start - (start - (start % 3)) + p - 1) % 3) + 1
            return str(frame)
        elif strand == "-":
            # Modified to return positive integer frame
            val = (chrom_size - end + 1)
            # awk logic was: 0 - (val % 3) - p.   if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}
            # To get a positive frame (1,2,3), we take the absolute or normalized value
            #frame = abs((val % 3) + p)
            frame = abs(val - val - (val%3) - p)
            if frame == 4: frame=1
            #return str(frame)
            return str(frame if frame != 0 else 3)
    except:
        return "NA"
    return "NA"

def main():
    parser = argparse.ArgumentParser(description="Convert GTF to BED with Genomic Frame")
    parser.add_argument("-i", "--input", required=True, help="Input GTF file")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    parser.add_argument("-s", "--species", default="human", help="Species name for genome report")
    args = parser.parse_args()

    chrom_sizes = get_chrom_sizes(args.species)
    cds_by_tx = defaultdict(list)
    gene_to_txs = defaultdict(set)
    exon_to_txs = defaultdict(set)

    # PASS 1: Parse GTF
    with open(args.input) as fh:
        for line in fh:
            if line.startswith("#"): continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9: continue

            chrom, _, feature, start, end, _, strand, phase, attrs_str = fields
            if feature != "CDS": continue

            cs = chrom_sizes.get(chrom, 0)
            frame = calculate_genomic_frame(int(start), int(end), strand, phase, cs)

            attrs = parse_attributes(attrs_str)
            gene = get_gene_name(attrs)
            tx = attrs.get("transcript_id", "NA")

            cds_by_tx[tx].append({
                "chrom": chrom, "start": int(start), "end": int(end),
                "gene": gene, "tx": tx, "frame": frame, "strand": strand
            })
            gene_to_txs[gene].add(tx)
            exon_key = (chrom, int(start), int(end), gene, frame)
            exon_to_txs[exon_key].add(tx)

    # PASS 2: Write BED
    with open(args.output, "w") as out:
        for tx, exons in cds_by_tx.items():
            strand = exons[0]["strand"]
            if strand == "+":
                exons.sort(key=lambda x: x["start"])
            else:
                exons.sort(key=lambda x: x["start"], reverse=True)

            total_exons = len(exons)
            for exon_number, e in enumerate(exons, start=1):
                gene = e["gene"]
                exon_key = (e["chrom"], e["start"], e["end"], gene, e["frame"])
                gene_tx_count = len(gene_to_txs[gene])
                exon_tx_count = len(exon_to_txs[exon_key])
                sharing = exon_tx_count / gene_tx_count

                exon_order = "first" if exon_number == 1 else ("last" if exon_number == total_exons else "middle")
                splicing = "constitutive" if sharing == 1.0 else ("unique" if exon_tx_count == 1 else "alternative")

                # Column 6: Strand, Column 7: Frame
                out.write(f"{e['chrom']}\t{e['start']}\t{e['end']}\t{gene}\t{e['tx']}\t"
                          f"{e['strand']}\t{e['frame']}\t{exon_number}\t{gene_tx_count}\t"
                          f"{exon_tx_count}\t{exon_order}\t{sharing:.3f}\t{splicing}\n")

if __name__ == "__main__":
    main()
