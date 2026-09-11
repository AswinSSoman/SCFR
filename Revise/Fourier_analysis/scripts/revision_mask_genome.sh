######################################################################################################################################################################################################################################################################################################
#																																																																																																																																				Mask human genome
######################################################################################################################################################################################################################################################################################################

#1. Mask human CDS in human genome 

#1. Prepare inputs

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
#2. Mask regions based on CDS annotations in the genome (hard mask)

	# GTF is 1-based; BED is 0-based half-open
	awk 'BEGIN{OFS="\t"} $3=="CDS"{print $1,$4-1,$5}' $GTF | sort -k1,1 -k2,2n | bedtools merge -i - > cds.bed

	bedtools maskfasta -fi $GENOME -bed cds.bed -fo genome.cds_masked.fa
	samtools faidx genome.cds_masked.fa

	#Get cds from ncbi to cross-check
	#datasets download genome accession GCA_054883195.1 --include rna,cds,protein,seq-report

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3. Mask regions based on CDS blast
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#3.1. Prepare CDS as query (spliced, per-transcript)

	#1m1.603s
	gffread -x cds_query.fa -g $GENOME $GTF

	#1. Dedupe your query (biggest change in size)
	#Human CDS from a GTF has massive redundancy — every isoform re-lists shared exons. Searching the same sequence 5–20× wastes most of your runtime.
	cd-hit-est -i cds_query.fa -o cds_query.nr.fa -c 1.0 -T 28 -M 0

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.2. Prepare BLAST database & 

	time makeblastdb -in genome.cds_masked.fa -dbtype nucl -out genome_masked_db -parse_seqids
	#For blastn/dc-megablast this speeds up seed lookup substantially, free accuracy-wise:
	#3m21.427s
	time makembindex -input genome_masked_db -iformat blastdb

	#Use -mt_mode 1 and split query across processes, not just threads
	#BLAST's internal multithreading doesn't scale linearly much past ~8–16 threads for query-side parallelism. Splitting the query file and running independent processes scales better:
	seqkit split2 -p 28 -O split_query cds_query.nr.fa
	time ls split_query/*.fa | parallel -j 28 blastn -query {} -db genome_masked_db -task blastn -evalue 1e-3 -use_index true -mt_mode 1 -outfmt 6 -out {}.blast.tsv
	time cat split_query/*.blast.tsv > cds_vs_masked.blast.tsv

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.3. Run BLAST CDS against the masked genome

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

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3.4. Further mask the CDS homology hits 

#144m37.402s
time awk 'BEGIN{OFS="\t"}{if($9<$10) print $2,$9-1,$10; else print $2,$10-1,$9}' cds_vs_masked.blast.tsv | sort -k1,1 -k2,2n | bedtools merge -i - > pseudo_hits.bed
time bedtools maskfasta -fi genome.cds_masked.fa -bed pseudo_hits.bed -fo genome.fully_masked.fa -mc n

#Print proportion of masked bases 
awk '/^>/ {next}
{   total += length($0)
    N += gsub(/N/, "", $0)
    n += gsub(/n/, "", $0)}
END {printf "Total bases : %d\n", total
    printf "N bases     : %d (%.4f%%)\n", N, 100*N/total
    printf "n bases     : %d (%.4f%%)\n", n, 100*n/total
    printf "N + n       : %d (%.4f%%)\n", N+n, 100*(N+n)/total}' genome.fully_masked.fa

scp pseudo_hits.bed genome.fully_masked.fa ceglab25@172.28.65.125:/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology/

cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology
cat cds.bed pseudo_hits.bed | sort -k1,1 -k2,2n | bedtools merge -i - > all_human_homology_regions.bed

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#4. Visualization of masked regins

cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology
time makeblastdb -in genome.fully_masked.fa -out genome.fully_masked.fa -dbtype nucl

#Transfer the relevant files
scp /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna \
	/media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf \
	/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/human_scfr_all_atleast_300bp*.bed \
	/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology/pseudo_hits.bed \
	/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology/all_human_homology_regions.bed
	/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology/genome.fully_masked.fa ceglab8@172.28.65.118:~/Downloads/SCFR

#Sort the GTF file
igvtools sort /home/ceglab8/Downloads/SCFR/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf /home/ceglab8/Downloads/SCFR/GCF_009914755.1_T2T-CHM13v2.0_genomic.sorted.gtf
#Generate the .idx index file
igvtools index /home/ceglab8/Downloads/SCFR/GCF_009914755.1_T2T-CHM13v2.0_genomic.sorted.gtf



