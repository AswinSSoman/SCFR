import sys

def parse_attributes(attr_str):
    """Parse the attributes column of a GTF line into a dictionary."""
    attrs = {}
    for attr in attr_str.strip().split(';'):
        if attr.strip() == "":
            continue
        key, value = attr.strip().split(' ', 1)
        attrs[key] = value.strip('"')
    return attrs

def gtf_to_bed(gtf_file, bed_file):
    with open(gtf_file, 'r') as gtf, open(bed_file, 'w') as bed:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = fields

            # Only include CDS features (coding exons)
            if feature != "CDS":
                continue

            attr_dict = parse_attributes(attributes)
            gene_name = attr_dict.get('gene_name', attr_dict.get('gene_id', 'unknown'))

            # BED is 0-based, half-open; GTF is 1-based, closed
            bed_start = int(start) - 1
            bed_end = int(end)

            bed.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene_name}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gtf_to_bed.py input.gtf output.bed")
        sys.exit(1)

    gtf_to_bed(sys.argv[1], sys.argv[2])
