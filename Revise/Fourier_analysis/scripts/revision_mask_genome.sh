######################################################################################################################################################################################################################################################################################################
#
######################################################################################################################################################################################################################################################################################################



#In ceglab25
#cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology
GENOME=/media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna
GTF=/media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf
#scp $GENOME $GTF ceglab27@172.28.65.127:~/aswin/SCFR/Fourier_analysis/human/300/Non_homology/

#In ceglab27
cd ~/aswin/SCFR/Fourier_analysis/human/300/Non_homology
GENOME=GCA_009914755.4_T2T-CHM13v2.0_genomic.fna
GTF=GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf

grep ">" $GENOME | head -3
cut -f1 $GTF | grep -v "^#" | sort -u | head -3

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#1. Mask CDS in the genome (hard mask)

# GTF is 1-based; BED is 0-based half-open
awk 'BEGIN{OFS="\t"} $3=="CDS"{print $1,$4-1,$5}' $GTF | sort -k1,1 -k2,2n | bedtools merge -i - > cds.bed

bedtools maskfasta -fi $GENOME -bed cds.bed -fo genome.cds_masked.fa
samtools faidx genome.cds_masked.fa

#Get cds from ncbi
#datasets download genome accession GCA_054883195.1 --include rna,cds,protein,seq-report

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#2. Prepare CDS as query (spliced, per-transcript)

#1m1.603s
gffread -x cds_query.fa -g $GENOME $GTF

#1. Dedupe your query (biggest win, free)
#Human CDS from a GTF has massive redundancy — every isoform re-lists shared exons. Searching the same sequence 5–20× wastes most of your runtime.
cd-hit-est -i cds_query.fa -o cds_query.nr.fa -c 1.0 -T 28 -M 0

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.Prepare BLAST database & Run BLAST CDS against the masked genome

time makeblastdb -in genome.cds_masked.fa -dbtype nucl -out genome_masked_db -parse_seqids
#For blastn/dc-megablast this speeds up seed lookup substantially, free accuracy-wise:
#3m21.427s
time makembindex -input genome_masked_db -iformat blastdb

#Use -mt_mode 1 and split query across processes, not just threads
#BLAST's internal multithreading doesn't scale linearly much past ~8–16 threads for query-side parallelism. Splitting the query file and running independent processes scales better:
seqkit split2 -p 28 -O split_query cds_query.nr.fa
time ls split_query/*.fa | parallel -j 28 blastn -query {} -db genome_masked_db -task blastn -evalue 1e-3 -use_index true -mt_mode 1 -outfmt 6 -out {}.blast.tsv
time cat split_query/*.blast.tsv > cds_vs_masked.blast.tsv

#Consider dc-megablast as a first pass, not -task blastn alone
#-task blastn (word_size 11, no discontiguous seeding) is the slowest sensitive option. dc-megablast uses discontiguous seeds — much faster, and it's specifically tuned for the 65–95% identity range where most pseudogenes and cross-paralog homology fall:
#450m11.209s
time blastn -query cds_query.nr.fa -db genome_masked_db -task dc-megablast -evalue 1e-3 -use_index true -outfmt 6 -num_threads 28 -out pass1.tsv

#Trade-off: dc-megablast is somewhat less sensitive than plain -task blastn for very short or highly diverged (<65% identity) hits — old, decayed pseudogenes could fall in that zone. 
#Safe strategy: run dc-megablast first (fast, catches the bulk), then run plain -task blastn only on CDS queries that got zero hits, since those are the ones at risk of being missed:
#5m28.969s
comm -23 <(grep ">" cds_query.nr.fa | sed 's/>//' | sort) <(cut -f1 pass1.tsv | sort -u) > no_hit_ids.txt
seqkit grep -n -f no_hit_ids.txt cds_query.nr.fa > cds_query.leftover.fa
#1130m43.787s
time blastn -query cds_query.leftover.fa -db genome_masked_db -task blastn -evalue 1e-3 -use_index true -outfmt 6 -num_threads 28 -out pass2.tsv
#32m19.647s
time cat pass1.tsv pass2.tsv > cds_vs_masked.blast.tsv

#Mask the homology hits
#144m37.402s
time awk 'BEGIN{OFS="\t"}{if($9<$10) print $2,$9-1,$10; else print $2,$10-1,$9}' cds_vs_masked.blast.tsv | sort -k1,1 -k2,2n | bedtools merge -i - > pseudo_hits.bed
time bedtools maskfasta -fi genome.cds_masked.fa -bed pseudo_hits.bed -fo genome.fully_masked.fa -mc n


awk '/^>/ {next}
{   total += length($0)
    N += gsub(/N/, "", $0)
    n += gsub(/n/, "", $0)}
END {printf "Total bases : %d\n", total
    printf "N bases     : %d (%.4f%%)\n", N, 100*N/total
    printf "n bases     : %d (%.4f%%)\n", n, 100*n/total
    printf "N + n       : %d (%.4f%%)\n", N+n, 100*(N+n)/total}' genome.fully_masked.fa

scp pseudo_hits.bed genome.fully_masked.fa ceglab25@172.28.65.125:/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology/


######################################################################################################################################################################################################################################################################################################


cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology
time makeblastdb -in genome.fully_masked.fa -out genome.fully_masked.fa -dbtype nucl


#Visualize the bed files in genome
scp /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna \
	/media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf \
	/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/human_scfr_all_atleast_300bp*.bed \
	/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology/pseudo_hits.bed \
	/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology/genome.fully_masked.fa ceglab8@172.28.65.118:~/Downloads/SCFR

#Sort the GTF file
igvtools sort /home/ceglab8/Downloads/SCFR/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf /home/ceglab8/Downloads/SCFR/GCF_009914755.1_T2T-CHM13v2.0_genomic.sorted.gtf

#Generate the .idx index file
igvtools index /home/ceglab8/Downloads/SCFR/GCF_009914755.1_T2T-CHM13v2.0_genomic.sorted.gtf


#Check overlap between intergenic SCFR & regions with homology
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300
bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique.bed -b Non_homology/pseudo_hits.bed | wc -l
bedtools intersect -a human_scfr_all_atleast_300bp_antisense_intergenic.bed -b Non_homology/pseudo_hits.bed | wc -l

#Save Intergenic SCFRs with no overlap with homolohy regions
bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique.bed -b Non_homology/pseudo_hits.bed -v > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed
bedtools intersect -a human_scfr_all_atleast_300bp_antisense_intergenic.bed -b Non_homology/pseudo_hits.bed -v > human_scfr_all_atleast_300bp_antisense_intergenic_with_no_homology.bed

#scp human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed human_scfr_all_atleast_300bp_antisense_intergenic_with_no_homology.bed ceglab8@172.28.65.118:~/Downloads/SCFR/


#Concatenate DFT results with pure intergenic SCFRs with no homology
time while read s
do
sid=$(echo -e "$s" | awk -F "\t" '{print$1,$2,$3,"frame"$4}' OFS="_" | tr "-" "_")
dtf=$(awk -v sid="$sid" '$1~sid' /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/human_scfr_all_atleast_300bp_only_intergenic_unique/chromosome_wise_summary/summary.tsv | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1')
echo $s $dtf
unset sid dtf
done < human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed	| sed 's/[ ]\+/\t/g' > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_dft_results.bed


time LC_ALL=C awk -F'\t' -v OFS='\t' '
# --- Pass 1: load summary.tsv into memory, keyed by its first column ---
NR==FNR {
    data[$1] = $0
    next
}

# --- Pass 2: stream the bed file, build matching key, look up in hash ---
{
    chrom = $1
    gsub(/\./, "_", chrom)          # NC_060925.1 -> NC_060925_1

    frame = $4
    if (frame ~ /^-/) {
        gsub(/-/, "_", frame)        # -2 -> _2
        key = "SCFR_" chrom "_" $2 "_" $3 "_frame" frame
    } else {
        key = "SCFR_" chrom "_" $2 "_" $3 "_frame" frame   # 1 -> frame1
    }

    lookup = key "_fft"
    out = $0
    if (lookup in data) {
        n = split(data[lookup], f, "\t")
        for (i = 1; i <= n; i++) if (f[i] == "") f[i] = "-"
        rest = f[1]
        for (i = 2; i <= n; i++) rest = rest "\t" f[i]
        out = out "\t" rest
    }
    print out
}
' /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/human_scfr_all_atleast_300bp_only_intergenic_unique/chromosome_wise_summary/summary.tsv \
  human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed \
  > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_dft_results.bed




