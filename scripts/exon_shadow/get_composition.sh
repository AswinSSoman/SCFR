#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input.fasta>"
    exit 1
fi

FASTA=$1

GC_AT=$(mktemp)
STATS=$(mktemp)

# -----------------------------------
# Step 1: GC / AT content (seqkit)
# -----------------------------------
seqkit fx2tab -n -g "$FASTA" \
| awk 'BEGIN{OFS="\t"} {print $1, $2, 100-$2}' \
> "$GC_AT"

# -----------------------------------
# Step 2: GC skew + GC3 + GC & AT average stretch
# -----------------------------------
perl -ne '
if (/^>(\S+)/) {
    process() if defined $name;
    $name=$1;
    $seq="";
    next;
}
chomp;
$seq.=$_;

END { process() if defined $name }

sub process {
    my $s = uc($seq);
    my $len = length($s);

    # --------------------
    # Base counts
    # --------------------
    my $G = ($s =~ tr/G//);
    my $C = ($s =~ tr/C//);

    my $gc_skew = ($G+$C) ? ($G-$C)/($G+$C) : 0;

    # --------------------
    # GC3 calculation
    # --------------------
    my ($gc3, $codon3) = (0, 0);
    for (my $i = 2; $i < $len; $i += 3) {   # 0-based index → 3rd codon position
        my $base = substr($s, $i, 1);
        if ($base =~ /[GC]/) {
            $gc3++;
            $codon3++;
        } elsif ($base =~ /[AT]/) {
            $codon3++;
        }
    }
    my $gc3_pct = $codon3 ? ($gc3/$codon3)*100 : 0;

    # --------------------
    # GC stretches
    # --------------------
    my @gc = ($s =~ /([GC]{3,})/g);
    my $gc_avg = 0;
    if (@gc) {
        my $sum=0;
        $sum+=length($_) for @gc;
        $gc_avg = $sum/@gc;
    }

    # --------------------
    # AT stretches
    # --------------------
    my @at = ($s =~ /([AT]{3,})/g);
    my $at_avg = 0;
    if (@at) {
        my $sum=0;
        $sum+=length($_) for @at;
        $at_avg = $sum/@at;
    }

    printf "%s\t%.5f\t%.2f\t%.2f\t%.2f\n",
        $name, $gc_skew, $gc3_pct, $gc_avg, $at_avg;
}
' "$FASTA" > "$STATS"

# -----------------------------------
# Step 3: Merge results
# -----------------------------------
awk 'BEGIN{OFS="\t"}
NR==FNR {a[$1]=$2"\t"$3; next}
{print $1, a[$1], $2, $3, $4, $5}
' "$GC_AT" "$STATS" \
| awk 'BEGIN{
    print "seqname\tGC\tAT\tGC_skew\tGC3\tGC_avg_stretch\tAT_avg_stretch"
}1'

# -----------------------------------
# Cleanup
# -----------------------------------
rm "$GC_AT" "$STATS"

