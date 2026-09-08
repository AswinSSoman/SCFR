######################################################################################################################################################################################################################################################################################################
#																																																																																																																																				Idenitfy Proto-genes
######################################################################################################################################################################################################################################################################################################

######################################################################################################################################################################################################################################################################################################
#Prepare inputs

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
#1. Mask CDS based on annotation in the genome (hard mask)

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
#3. Prepare BLAST database & 

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
#4. Run BLAST CDS against the masked genome

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
#5. Further mask the CDS homology hits 

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

######################################################################################################################################################################################################################################################################################################
#Visualization

cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology
time makeblastdb -in genome.fully_masked.fa -out genome.fully_masked.fa -dbtype nucl

#Transfer the relevant files
scp /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna \
	/media/aswin/SCFR/SCFR-main/genes/human/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf \
	/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/human_scfr_all_atleast_300bp*.bed \
	/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology/pseudo_hits.bed \
	/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Non_homology/genome.fully_masked.fa ceglab8@172.28.65.118:~/Downloads/SCFR

#Sort the GTF file
igvtools sort /home/ceglab8/Downloads/SCFR/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf /home/ceglab8/Downloads/SCFR/GCF_009914755.1_T2T-CHM13v2.0_genomic.sorted.gtf

#Generate the .idx index file
igvtools index /home/ceglab8/Downloads/SCFR/GCF_009914755.1_T2T-CHM13v2.0_genomic.sorted.gtf

######################################################################################################################################################################################################################################################################################################
#Intergenic SCFRs

	#Check overlap between intergenic SCFR & regions with homology
	cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300
	bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique.bed -b Non_homology/pseudo_hits.bed | wc -l
	bedtools intersect -a human_scfr_all_atleast_300bp_antisense_intergenic.bed -b Non_homology/pseudo_hits.bed | wc -l

	#Save Intergenic SCFRs with no overlap with homolohy regions
	bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique.bed -b Non_homology/pseudo_hits.bed -v > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed
	bedtools intersect -a human_scfr_all_atleast_300bp_antisense_intergenic.bed -b Non_homology/pseudo_hits.bed -v > human_scfr_all_atleast_300bp_antisense_intergenic_with_no_homology.bed

	#scp human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed human_scfr_all_atleast_300bp_antisense_intergenic_with_no_homology.bed ceglab8@172.28.65.118:~/Downloads/SCFR/

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Concatenate DFT results with pure intergenic SCFRs with no homology

time while read s
do
sid=$(echo -e "$s" | awk -F "\t" '{print$1,$2,$3,"frame"$4}' OFS="_" | tr "-" "_")
dtf=$(awk -v sid="$sid" '$1~sid' /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/human_scfr_all_atleast_300bp_only_intergenic_unique/chromosome_wise_summary/summary.tsv | awk -F"\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1')
echo $s $dtf
unset sid dtf
done < human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed	| sed 's/[ ]\+/\t/g' > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_dft_results.bed

time LC_ALL=C awk -F'\t' -v OFS='\t' 'NR==FNR {
    data[$1] = $0
    next
}{
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
}' /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/human_scfr_all_atleast_300bp_only_intergenic_unique/chromosome_wise_summary/summary.tsv human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_dft_results.bed


#View in IGV
/home/ceglab8/workspace/1.master_thesis/6.source_codes/IGV_Linux_2.5.0/igv.sh

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Filter only SCFRs with 0.3 DFT frequency

awk -F "\t" '$10~/0\.3/' human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_dft_results.bed > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results.bed
awk -F "\t" '$10~/0\.33/' human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_dft_results.bed > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed

#Visualize length
awk '{print$3-$2}' human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results.bed | statplot.R --biplot human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results_length_distribution.png
awk '{print$3-$2}' human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed | statplot.R --biplot human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_length_distribution.png
awk -F "\t" '{print$3-$2}' human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed | sort | uniq -c | sed 's/^[ ]\+//g' | tr " " "\t" | ./term_hist2.py -s 50 -H 30 | less -SR

#sort based on length
awk -F "\t" 'BEGIN{OFS="\t"} {print ($3 - $2), $0}' human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results.bed | sort -k1,1nr | cut -f2- > length_sorted_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results.bed

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Visualize the intergenic SCFR with 0.33 DFT frequency

cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/claude
#Download cytoband track data 
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.cytoBandMapped/chm13v2.0_cytobands_allchrs.bed.gz
#All data related to this genome is available at: https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/
gzip -d chm13v2.0_cytobands_allchrs.bed.gz

#Set up python environment
conda activate scfr
which python3
python3 -c "import numpy, pandas; print(numpy.__version__, pandas.__version__)"

time python3 scfr_report.py \
	../filtering_intergenic_SCFR/human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed \
	-o report_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.html \
	--fai /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.fai \
	--cytoband chm13v2.0_cytobands_allchrs.bed \
	--main-only

time python3 scfr_report.py \
	../filtering_intergenic_SCFR/human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene.bed \
	-o report_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene.html \
	--fai /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.fai \
	--cytoband chm13v2.0_cytobands_allchrs.bed \
	--main-only

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Add header
sed -i '1i Chr Start End Frame_strand Filler Strand SCFR_Name Num_Raw_Peaks Top_Peak_Count Frequencies Magnitudes Periods' length_sorted_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results.bed 
#make tab separated
sed -i 's/[ ]\+/\t/g' length_sorted_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results.bed
#Make it standard bed format: Take 1st 6 columns
awk -F "\t" '{print$1,$2,$3,$4,$5,$6}' OFS="\t" human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results.bed > ucsc.bed

#scp human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results.bed ucsc.bed ceglab8@172.28.65.118:~/Downloads/SCFR/
#scp human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed ceglab8@172.28.65.118:~/Downloads/SCFR/

######################################################################################################################################################################################################################################################################################################
#Download NCBI Refseq Functional elements for human T2T genome from UCSC

wget https://ftp.ncbi.nlm.nih.gov/refseq/FunctionalElements/trackhub/data/human/GCF_009914755.1-RS_2023_10/FEbiolregions_RS_2023_10_T2T-CHM13v2.0.bb
bigBedToBed FEbiolregions_RS_2023_10_T2T-CHM13v2.0.bb FEbiolregions_RS_2023_10_T2T-CHM13v2.0.bed

#Rename chromosome IDs
awk -F "\t" 'FNR==NR {map[$1]=$2; next} $1 in map {$1=map[$1]} 1' OFS="\t" /media/aswin/SCFR/SCFR-main/genomes/human/map.tsv FEbiolregions_RS_2023_10_T2T-CHM13v2.0.bed > tempFE.bed
sed -i 's/\r//' tempFE.bed
awk -F'\t' 'BEGIN{OFS="\t"} {
    # Combine everything from column 10 onwards into a single block of text
    desc = $10
    for(i=11; i<=NF; i++) {
        if($i != "") desc = desc "; " $i
    }
    # Remove hidden Windows carriage returns and internal line breaks from the text
    gsub(/\r/, "", desc)
    gsub(/\n/, " ", desc)
    
    # Print a standard 12-column BED format
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,($3-$2),0,desc
}' tempFE.bed > NCBI_RefSeq_Functional_Elements_GCF_009914755.1-RS_2023_10_biological_regions.bed 

#scp NCBI_RefSeq_Functional_Elements_GCF_009914755.1-RS_2023_10_biological_regions.bed ceglab8@172.28.65.118:~/Downloads/SCFR/
#scp NCBI_RefSeq_Functional_Elements_GCF_009914755.1-RS_2023_10_biological_regions_fomratted.bed ceglab8@172.28.65.118:~/Downloads/SCFR/

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download RepeatMasker data for human from UCSC

cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300
wget https://hgdownload.soe.ucsc.edu/gbdb/hs1/t2tRepeatMasker/chm13v2.0_rmsk.bb
bigBedToBed chm13v2.0_rmsk.bb chm13v2.0_rmsk.bed

awk -F'"' '{
    refseq = ""; ucsc = "";
    for(i=1; i<=NF; i++) {
        if($i == "refseqAccession") refseq = $(i+2)
        if($i == "ucscStyleName") ucsc = $(i+2)
    }
    if(refseq != "" && ucsc != "") print  ucsc "\t" refseq
}' /media/aswin/SCFR/SCFR-main/genomes/human/sequence_report.jsonl > refseq_to_chrnames.tsv

awk -F "\t" 'FNR==NR {map[$1]=$2; next} $1 in map {$1=map[$1]} 1' OFS="\t" refseq_to_chrnames.tsv chm13v2.0_rmsk.bed > Repetitive_Elements.bed

#scp Repetitive_Elements.bed ceglab8@172.28.65.118:~/Downloads/SCFR/

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Download RNA-seq track data from UCSC

wget https://hgdownload.soe.ucsc.edu/gbdb/hs1/bbi/xenoRefGene.bb

######################################################################################################################################################################################################################################################################################################
#Prepare SCFR classes

#Some borderline cases occurs i.e. some SCFRs sometimes overlaps with both intergenic & genic features. But since the overlap length with intergenic reigion is >300bp these would be saved inside "human_scfr_all_atleast_300bp_only_intergenic_unique.bed"
#Partially these SCFRs overlap with genic regions, importantly CDS, which contributes to it's 0.33 DFT frequency, hence these can't be proto-genes & must be removed
#Hence identify all such SCFRs (in the intergenic region which has no homology with human CDS & has DFT frequency 0.33) that has overlap with genic regions

bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed -b final_gtf_features/genes.bed -wo | awk '{print$NF}' | ministat 

cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300
mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/annotation_overlaping_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results

bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed -b final_gtf_features/genes.bed -wo > annotation_overlaping_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results/overlapping_genes.bed
bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed -b final_gtf_features/cds_unique.bed -wo > annotation_overlaping_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results/overlapping_cds_unique.bed
bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed -b final_gtf_features/pseudogene_unique.bed -wo > annotation_overlaping_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results/overlapping_pseudogene_unique.bed
bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed -b final_gtf_features/RNA_genes_unique.bed -wo > annotation_overlaping_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results/overlapping_RNA_genes_unique.bed
bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed -b final_gtf_features/utr_refined_unique.bed -wo > annotation_overlaping_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results/overlapping_utr_refined_unique.bed
bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed -b final_gtf_features/only_introns_refined_unique.bed -wo > annotation_overlaping_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results/overlapping_only_introns_refined_unique.bed
bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed -b final_gtf_features/only_intergenic_unique.bed -wo > annotation_overlaping_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results/overlapping_only_intergenic_unique.bed

#check overlaps between bed files
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/annotation_overlaping_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results
./bed_overlap_terminal.py overlapping_cds_unique.bed overlapping_genes.bed overlapping_only_intergenic_unique.bed overlapping_only_introns_refined_unique.bed overlapping_pseudogene_unique.bed overlapping_RNA_genes_unique.bed overlapping_utr_refined_unique.bed
./bed_overlap_matrix.py overlapping_cds_unique.bed overlapping_genes.bed overlapping_only_intergenic_unique.bed overlapping_only_introns_refined_unique.bed overlapping_pseudogene_unique.bed overlapping_RNA_genes_unique.bed overlapping_utr_refined_unique.bed

#NOTE: SCFRs overlapping genes are overlapping with all SCFRs overlapping other features like UTR, CDS, introns, RNA-genes, etc

bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed -b final_gtf_features/genes.bed -v > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene.bed
awk -F "\t" 'BEGIN{OFS="\t"} {print ($3 - $2), $0}' human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene.bed | sort -k1,1nr | cut -f2- > length_sorted_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results.bed

#scp human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene.bed ceglab8@172.28.65.118:~/Downloads/SCFR/

######################################################################################################################################################################################################################################################################################################
#BLAST AGAINST DATABASES OF EXISTING PROTEINS/GENES

#merge SCFRs as many filtered SCFRs are overlapping with short shift in frame
bedtools merge -i human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene.bed > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged.bed

mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/nr_blast
mv human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged.bed nr_blast/
#Get fasta
cd nr_blast/
bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged.bed -name+ > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged.fa

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Strategy:


#NCBI has mainly 2 databases which serve as  ultimate, all-inclusive catch-all archives [1, 2] of public sequence data:
	#1. nr (non-redundant protein): Stores translated amino acid sequences from GenBank, RefSeq, UniProt, and other archives.
	#2. nt (non-redundant nucleotide): Stores DNA and RNA nucleic acid sequences.

#The core message is: To prove a gene is truly brand new (de novo), you must prove it has absolutely no ancestral relatives. Searching the protein database (nr) is a much harsher, more decisive test for this than searching the nucleotide database (nt).
------------------------------
## 1. The Core Problem: Why Nucleotide Search (nt) Fails for De Novo Genes
	#* DNA changes too fast: DNA sequences mutate rapidly. Due to codon degeneracy (different DNA triplets coding for the exact same amino acid), two genes can look completely different at the DNA level but still produce the exact same protein.
	#* The "False De Novo" Trap: If an ancient gene mutates heavily at the DNA level, a blastn search against nt will find zero matches. You might falsely celebrate thinking you found a brand-new de novo gene. However, it’s actually just a fast-evolving remnant of an old gene.
#2. Why Protein Search (nr) is the "Decisive Test"
	#* Deep evolutionary sight: Protein-level searches use substitution matrices (like BLOSUM62), which know that certain amino acids behave similarly (e.g., swapping Isoleucine for Leucine is a "chemically silent" change).
	#* The True Test: By translating your query and searching against nr (using a tool like blastx), you can detect incredibly faint, ancient evolutionary signals that raw DNA comparisons (blastn) completely miss. In de novo gene literature, if your sequence matches anything in nr, it is immediately excluded from being a true de novo gene.
## 3. Why Use a Lenient E-value (1e-3 instead of 1e-5)?
	#* Discovery vs. Exclusion: When looking for a specific gene, you want a strict E-value (like $10^{-5}$ or $10^{-10}$) so you only get certain, high-quality matches.
	#* Erring on the side of caution: Here, your goal is exclusion—you want to catch any potential hint of an old relative. By loosening the threshold to $10^{-3}$ ($1e-3$), you are telling BLAST: "Show me even the faintest, most borderline matches." A borderline match means your candidate is suspicious and should be discarded, not kept.
## 4. Why Use -seg yes?
	#* Masking "Junk" Repetitions: Many non-coding regions have low-complexity sequences (like ATATATATAT or AAAAAA).
	#* Preventing Fake Matches: If your candidate gene has a repetitive string of text, it might accidentally match a repetitive string in a completely unrelated gene. Turning -seg yes on "blurs out" these low-complexity zones so BLAST only focuses on meaningful sequence data.
## 5. Why use nt only as a Secondary Check?
	#The text notes that nt shouldn't be ignored entirely, but it serves a different purpose. A secondary blastn against nt is used to catch:
	#* RNA genes: Functional genes that never turn into proteins (like tRNAs or lncRNAs), which won't show up in a protein (nr) database.
	#* Pseudogenes: "Dead" genes that have suffered mutations throwing off their reading frames (no clear Open Reading Frame), meaning they can no longer be cleanly translated into a protein sequence but still share DNA sequence identity.
	
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#In mopheus

#Get some space
cd ~/pavo
find . -type f ! \( -name "*.sh" -o -name "*.R" -o -name "*.py" -o -name "*.pl" \) -delete
cd ~/Dros_rna
find . -type f ! \( -name "*.sh" -o -name "*.R" -o -name "*.py" -o -name "*.pl" \) -delete

#Install latest blast+
	cd
	wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
	mkdir -p ~/tools
	tar -xzf ncbi-blast-2.17.0+-x64-linux.tar.gz -C ~/tools/

#Download nr database
	#space required before downloading
	mkir ~/blastdb_new/nr
	cd ~/blastdb_new/nr
	curl -s https://ftp.ncbi.nlm.nih.gov/blast/db/nt-nucl-metadata.json | python3 -m json.tool
	#315m0.730s
	nohup ~/blastdb_new/nr/download_nr.sh > ~/blastdb_new/nr/nr_download.log 2>&1 &

#Download nt database
	mkir ~/blastdb_new/nt
	cd ~/blastdb_new/nt
	curl -s https://ftp.ncbi.nlm.nih.gov/blast/db/nt-nucl-metadata.json | python3 -m json.tool | grep -E "bytes-total"
	nohup ~/blastdb_new/nt/download_nt.sh > ~/blastdb_new/nt/nt_download.log 2>&1 &

#Prepare inputs
	mkdir -p ~/aswin/SCFR/Fourier_analysis/nr_blast
	cd ~/aswin/SCFR/Fourier_analysis/nr_blast
	scp ceglab25@172.28.65.125:/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/nr_blast/human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged.* .

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#RUN BLAST

#Primary test — translated nucleotide vs protein database:
	#tblastn is for protein query vs translated nt — you want the reverse: nucleotide query, 6-frame translated, vs protein db. That's blastx.
	#Use nr (protein), not nt, as the database here — this is the actual test the field uses for "no homology to known coding genes."
	#E-value 1e-3, not the usual 1e-5/1e-10 — you're doing exclusion, not discovery, so you want to be lenient about what counts as a hit. A borderline hit here is a reason to flag/exclude a candidate, not something you want to miss by being too strict. Err toward sensitivity.
	#-seg yes masks low-complexity regions, important since your regions may have compositional bias that could create spurious hits.

#1482m47.609s
nohup bash -c 'time /home/morpheus/tools/ncbi-blast-2.17.0+/bin/blastx \
 -query human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged.fa \
 -db /home/morpheus/blastdb_new/nr/nr \
 -out human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx.tsv \
 -outfmt "6 qseqid sseqid qlen length qstart qend evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand" \
 -evalue 1e-3 \
 -num_threads 56 \
 -max_target_seqs 10 \
 -seg yes' &> stdout_blastx_run.out &

#Inspect the blast hits
cd ~/aswin/SCFR/Fourier_analysis/nr_blast
#Quick view blast results with header names
head human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx.tsv | sed '1i Query Subject Query_Nt_length Alignment_Aa_length Q_Nt_start Q_Nt_end E_value Bit_score Raw_score %_Query_covered_per_sub %_Query_covered_per_hsp %_ident Matches Mismatches Gaps Strand' | column -t
transeq -frame 6 --auto --stdout \
#Cross check individually if SCFR sequence gives same hits using ncbi website blastx & matching exact sequence co-ordinate with local bed file 
transeq --frame 6 --auto --stdout <(bedtools getfasta -fi GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed <(grep XP_055898630.1 human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx.tsv | grep 1878 | awk '{print$1,$5,$6}'| sed 's/^:://g' | tr ":-" "\t" | awk '{print$1,$2+$5, $2+$4}' OFS="\t") -name+) | sed '/^>/! s/.*/\U&/' | ~/aswin/programmes/myfasta -comb | egrep "GLTAALHRVPIASASPHGVPIASASPHGVPIASASPHGVPIASASPHGVPIASATPHGVP|VASASPHGVPVASASPHGVPVASASPHGVPVASASPHGVPVASASPHGVPVASATPHGVPVASATPHGVPVAS" -z
grep XP_055898630.1 human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx.tsv | grep 1809 | ./check_scfr_fasta.sh  | egrep "ASPHGVPVASASPHGVPVASASPHGVPVASATPHGVPVASATPHGVPVASQ|ASPHGVPIASASPHGVPIASASPHGVPIASATPHGVPIASASPHGVPIASATPHGVPIAS" -z

#Convert blastx hits to genomic bed
awk '{
    sub(/^::/, "", $1)
    gsub(/[:-]/, " ", $1)
    split($1, a)
    print a[1] "\t" a[2] "\t" a[3] "\tSub:" $2 ",Len:" $4 ",Qcov:" $11 ",%id:" $12 ",Eval:" $7
}' human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx.tsv > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx_all_hits.bed

#More concise table of blastx hits (one subject hit at same location is shown only once)
awk '
{
    # 1. Parse chromosome, start, and end from $1
    sub(/^::/, "", $1)
    gsub(/[:-]/, " ", $1)
    split($1, a)
    chr = a[1]; start = a[2]; end = a[3]

    # 2. Extract values needed for key, evalue, and output formatting
    sub_id = $2
    len    = $4
    qcov   = $11
    pid    = $12
    evalue = $7 + 0  # Cast to numeric float for scientific notation comparison

    # 3. Create unique group key (cols 1, 2, 3, and 4)
    key = chr "\t" start "\t" end "\tSub:" sub_id

    # 4. Count total occurrences for this group key
    count[key]++

    # 5. Store/update record if key is new OR if a lower E-value is found
    if (!(key in min_eval) || evalue < min_eval[key]) {
        min_eval[key] = evalue
        best_row[key] = chr "\t" start "\t" end "\tSub:" sub_id \
                        ",Len:" len ",Qcov:" qcov ",%id:" pid ",Eval:" $7
    }
}
END {
    # Print the best row for each key along with the total hit count appended
    for (k in best_row) {
        print best_row[k] "\tHits:" count[k]
    }
}' human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx.tsv > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx_unique_hits.bed

sort -u human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx_unique_hits.bed > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx_hits_unique.bed

scp check_scfr_fasta.sh \
human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx.tsv \
human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx_all_hits.bed \
human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx_unique_hits.bed ceglab25@172.28.65.125:/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/nr_blast/

scp human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx.tsv \
human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx_all_hits.bed \
human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx_unique_hits.bed ceglab8@172.28.65.118:~/Downloads/SCFR

#scp human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx.tsv  human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx_hits.bed stdout_blastx_run.out ceglab25@172.28.65.125:/media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/nr_blast/

#In ceglab25
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300
bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene.bed -b nr_blast/human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx_hits.bed | wc -l
bedtools intersect -a human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene.bed -b nr_blast/human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged_results_blastx_hits.bed -v > human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_and_nrblastx_hits.bed

mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/filtering_intergenic_SCFR
cp human_scfr_all_atleast_300bp_only_intergenic_unique.bed \
human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed \
human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_dft_results.bed \
human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results.bed \
human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene.bed \
human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_and_nrblastx_hits.bed \
nr_blast/human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_merged.bed filtering_intergenic_SCFR/

awk -F "\t" 'BEGIN{OFS="\t"} {print ($3 - $2), $0}' human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results.bed | sort -k1,1nr | cut -f2- > length_sorted_human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.3_dft_results.bed

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Secondary test — nucleotide vs nt (catches RNA genes, pseudogenes, recent duplicates that blastx might miss if there's no ORF-preserving frame):
blastn -query all_559.fasta -db nt -out results_blastn.tsv -outfmt 6 -evalue 1e-3 -num_threads 72 -max_target_seqs 10


#scp /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna morpheus@172.30.1.121:~/aswin/SCFR/Fourier_analysis/nr_blast/

blastn -query all_559.fasta -db nt -out results.tsv -outfmt 6 -evalue 1e-10 -num_threads 16 -max_target_seqs 5
blastn_new -query all_559.fasta -out results.tsv -outfmt 6 -evalue 1e-10 -num_threads 72
blastx -query all_559.fasta -db nr -out results_blastx.tsv -outfmt 6 -evalue 1e-3 -num_threads 72 -max_target_seqs 10 -seg yes


