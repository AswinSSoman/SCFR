#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input.fasta>"
    exit 1
fi

FASTA=$1

# -----------------------------------
# Step 1: GC / AT content (seqkit)
# -----------------------------------
seqkit fx2tab -n -g "$FASTA" \
| awk 'BEGIN{OFS="\t"} {print $1, $2, 100-$2}' \
> gc_at.tmp

# -----------------------------------
# Step 2: Average GC stretch (Perl)
# -----------------------------------
perl -ne '
BEGIN {
    print "";
}
if (/^>(\S+)/) {
    process() if $name;
    $name=$1;
    $seq="";
    next;
}
$seq.=$_;
END { process() }

sub process {
    my @gc = ($seq =~ /([GC]{3,})/ig);
    if (@gc) {
        my $sum=0;
        $sum+=length($_) for @gc;
        printf "%s\t%.2f\n", $name, $sum/@gc;
    } else {
        printf "%s\t0\n", $name;
    }
}
' "$FASTA" > gc_stretch.tmp

# -----------------------------------
# Step 3: Merge results
# -----------------------------------
awk 'BEGIN{OFS="\t"}
NR==FNR {a[$1]=$2"\t"$3; next}
{print $1, a[$1], $2}
' gc_at.tmp gc_stretch.tmp \
| awk 'BEGIN{print "seqname\tGC\tAT\tGC_avg_stretch"}1'

# -----------------------------------
# Cleanup
# -----------------------------------
rm gc_at.tmp gc_stretch.tmp

