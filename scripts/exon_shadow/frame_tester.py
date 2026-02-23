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
            print(val)
            # awk logic was: 0 - (val % 3) - p.   if($7=="-") print$0,(a-$5+1)-(a-$5+1)-((a-$5+1)%3)-$8}
            # To get a positive frame (1,2,3), we take the absolute or normalized value
            #frame = abs((val % 3) + p)
            frame = abs((val-(val-(val%3))-p))
            #return str(frame)
            return str(frame if frame != 0 else 3)
    except:
        return "NA"
    return "NA"

a = calculate_genomic_frame(39768302, 39769239, "-", 2, 209442601)
print(a)
