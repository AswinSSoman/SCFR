





cd /media/aswin/SCFR/SCFR-main/genes/gibbon				
awk '!/^#/ {print $2}' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf | sort | uniq -c
awk '!/^#/ {print $3}' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf | sort | uniq -c
awk -F ";" '{for(n=1;n<=NF;n++) if($n~"gbkey") print $n}' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf | sort | uniq -c
awk -F ";" '{for(n=1;n<=NF;n++) if($n~"gene_biotype") print $n}' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf | sort | uniq -c


awk '!/^#/ {print $2}' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf | sort | uniq -c 
   9315 cmsearch
2654096 Gnomon
    126 RefSeq
   1307 tRNAscan-SE

awk '!/^#/ {print $3}' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf | sort | uniq -c 
1091668 CDS
1249900 exon
  38436 gene
  89434 start_codon
  89482 stop_codon
 105924 transcript

awk -F ";" '{for(n=1;n<=NF;n++) if($n~"gbkey") print $n}' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf | sort | uniq -c
1270584  gbkey "CDS"
     10  gbkey "C_region"
  38436  gbkey "Gene"
   4304  gbkey "misc_RNA"
  89208  gbkey "mRNA"
  11709  gbkey "ncRNA"
    155  gbkey "rRNA"
    448  gbkey "tRNA"
     90  gbkey "V_segment"


awk -F ";" '{for(n=1;n<=NF;n++) if($n~"gene_biotype") print $n}' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf | sort | uniq -c
     10  gene_biotype "C_region"
   6590  gene_biotype "lncRNA"
      2  gene_biotype "misc_RNA"
     44  gene_biotype "ncRNA"
  22115  gene_biotype "protein_coding"
   6322  gene_biotype "pseudogene"
    153  gene_biotype "rRNA"
   1118  gene_biotype "snoRNA"
   1535  gene_biotype "snRNA"
     31  gene_biotype "transcribed_pseudogene"
    426  gene_biotype "tRNA"
     90  gene_biotype "V_segment"


#
cd /media/aswin/SCFR/SCFR-main/genes/gibbon/test
wc -l GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.bed

head -3 GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.bed
NC_014047.1	0	70	unassigned_gene_1	.	+	RefSeq	exon	.	gene_id "unassigned_gene_1"; transcript_id "unassigned_transcript_527"; product "tRNA-Phe"; transcript_biotype "tRNA"; exon_number "1"; 
NC_014047.1	0	70	unassigned_gene_1	.	+	RefSeq	transcript	.	gene_id "unassigned_gene_1"; transcript_id "unassigned_transcript_527"; gbkey "tRNA"; product "tRNA-Phe"; transcript_biotype "tRNA"; 
NC_014047.1	70	1022	unassigned_gene_2	.	+	RefSeq	exon	.	gene_id "unassigned_gene_2"; transcript_id "unassigned_transcript_528"; note "12S ribosomal RNA"; product "s-rRNA"; transcript_biotype "rRNA"; exon_number "1"; 

head -3 ../../../exon_shadow/gibbon/gibbon_scfr_all.bed 
NC_072423.2	0	338	-3	1	-
NC_072423.2	0	339	1	1	+
NC_072423.2	0	4	-1	1	-

#Using scripts from claude 
mkdir output
cp classify_scfr.py extract_annotation.py run_classify.sh output/
#0m14.097s
time python3 output/extract_annotation.py GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.bed output
#46m21.848s
time output/run_classify.sh GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.bed gibbon_scfr_all.bed output/
#Took too much ram & system froze (absolute rubbish)
time python3 output/classify_scfr.py output/


#16m35.547s
time ./scfr_overlap.sh GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.bed gibbon_scfr_all.bed output2.bed

If there is absolutely no overlap between SCFR & annotation, then keep the entry with feaure as "intergenic" & biotypes as "intergenic" & 
If there is overlap check feature type, there are total 5 major types I'm interested in: 
	1. protein coding/CDS
	2. RNA genes: many types
	3. Pseudogenes
	4. introns
	5. anyother region within whole gene region except exons & introns.
If SCFR is overlapping with any of this 5 types, keep the output of bedtools intersect -wao only if the strand & frame is same between 2 bed files & save  
If SCFR is overlapping with any of this 5 types, but in different frame 

If there is overlap but in different strand irrespective of any frame, then keep feature as "gene" but biotype as "antisense"
If there is overlap but in different frame then also keep features as "gene" but biotype as "off_frame"
This code don
If an SCFR & feature has same coordinates, strand & frame but are duplicated as different features such as gene, exon, transcript, only use gene, not exon or transcript.
 
If a feature is biotype
group all RNA biotypes together irrespective of features

The output2 line count shouldn't be higher than scfr, unless a an scfr can be classified inside a protein coding/cds & non-coding part of gene & UTR
wc -l output2.bed 
1189017210 output2.bed
wc -l gibbon_scfr_all.bed 
329013580 gibbon_scfr_all.bed


#I want to all SCFRs of atleast 300bp, for running Discrete fourier transform analysis (DFT) to find if SCFRs from non-genic DNA (intergenic, intronic) has any specific patterns of signals, if proto-gene stage as a spectrum, shows any spectrum or classes of signals from complete non-genic to proto-genes. For that I need to classify SCFRs into coding portion of gene (CDS), intergenic, intronic, parts within gene other than exon & introns, pseudogenes, & RNA genes. But complication arise when considering strand & frame. If an SCFR is overlapping a CDS, it cpould be divided into further 6 types, as there are 2 strands & 3 frames each. Conceptually when we change frame not strand, the composition remains same, & DFT simply looks for periodicity, 

#I want to all SCFRs of atleast 300bp, for running Discrete fourier transform analysis (DFT) to find if SCFRs from non-genic DNA (intergenic, intronic) has any specific patterns of signals, if proto-gene stage as a spectrum, shows any spectrum or classes of signals from complete non-genic to proto-genes. For that I need to classify SCFRs into coding portion of gene (CDS), intergenic, intronic, parts within gene other than exon & introns, pseudogenes, & RNA genes. But complication arise when considering strand & frame. If an SCFR is overlapping a CDS, it cpould be divided into further 6 types, as there are 2 strands & 3 frames each. Conceptually when we change frame not strand, the composition remains same, & DFT simply looks for periodicity,

#########################################################################################################################################################################################################################################################################################################


#Filter SCFRs
#3m43.931s
time awk -F "\t" '($3-$2)>=300' gibbon_scfr_all.bed > gibbon_scfr_all_atleast_300bp.bed  

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Filter GTF

#Get relevant info from gtf in bed format
#1m58.695s
time awk -f gtf_to_bed.awk GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf > GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.bed

#Extract features
awk '$4=="CDS"' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.bed > cds.bed
awk '!seen[$1,$2,$3,$4,$5,$6]++' cds.bed > cds_unique.bed

awk '$7=="pseudogene"' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.bed > pseudogene.bed
awk '!seen[$1,$2,$3,$4,$5,$6]++' pseudogene.bed > pseudogene_unique.bed 

awk '$7~"RNA"' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.bed > RNA_genes.bed
awk '!seen[$1,$2,$3,$4,$5,$6]++' RNA_genes.bed > RNA_genes_unique.bed

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get introns & UTRs bed from GTF

#Install tool agat
conda create -n agat_env -c bioconda -c conda-forge agat
unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
conda activate agat_env

#Extract all features including introns & UTRs
#18m36.006s
time agat_sp_add_introns.pl --gff GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf --out introns.bed
#keep only introns
awk -F "\t" '$3=="intron"' introns.bed > only_introns.bed

#Keep relevant columns
awk 'BEGIN{OFS="\t"}
{gene_id=Parent=".";
    n=split($9,a,";");
    for(i=1;i<=n;i++){
        split(a[i],b,"=");
        if(b[1]=="gene_id" && b[2]!="") gene_id=b[2];
        if(b[1]=="Parent" && b[2]!="") Parent=b[2];
    }
    print $1,$4,$5,$3,1,$7,gene_id,Parent
}' only_introns.bed > only_introns_refined.bed

#remove duplicates
awk '!seen[$1,$2,$3,$4,$5,$6]++' only_introns_refined.bed > only_introns_refined_unique.bed

#keep only UTRs
time awk -F "\t" '$3~"_prime_UTR"' introns.bed > utr.bed
#Keep relevant columns
awk 'BEGIN{OFS="\t"}
{gene_id=Parent=".";
    n=split($9,a,";");
    for(i=1;i<=n;i++){
        split(a[i],b,"=");
        if(b[1]=="gene_id" && b[2]!="") gene_id=b[2];
        if(b[1]=="Parent" && b[2]!="") Parent=b[2];
    }
    print $1,$4,$5,$3,1,$7,gene_id,Parent
}' utr.bed > utr_refined.bed
#remove duplicates
awk '!seen[$1,$2,$3,$4,$5,$6]++' utr_refined.bed > utr_refined_unique.bed

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get intergenic regions

awk 'OFS="\t" {print $1, "0", $2}' chromSizes.txt | sort -k1,1 -k2,2n > chromSizes.bed
#8m23.934s
time agat_sp_add_intergenic_regions.pl --gff GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf --out intergenic.bed --cpu 28

time awk -F "\t" '$3~"intergenic_region"' intergenic.bed > only_intergenic.bed
#Keep relevant columns
awk -F "\t" '!/^#/ {print$1,$4,$5,$9}' only_intergenic.bed | awk '!seen[$1,$2,$3,$4]++' | tr " " "\t" > only_intergenic_unique.bed


#Simply check
#total genome coverage by gene
awk -F "\t" '$4=="gene"' GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.bed | awk '!seen[$1,$2,$3,$6]++' | awk '{sum+=($3-$2-1);} END{print sum}'
#total genome coverage by intergenic
awk '{sum+=($3-$2);} END{print sum}' only_intergenic_unique.bed


mkdir final_gtf_features
mv cds_unique.bed pseudogene_unique.bed RNA_genes_unique.bed only_introns_refined_unique.bed utr_refined_unique.bed only_intergenic_unique.bed final_gtf_features/
cd final_gtf_features/

#Find strand aware overlaps (minimum 300bp)

#with cds (last column tells % of total SCFR length covered by overlap)
time bedtools intersect -a gibbon_scfr_all_atleast_300bp.bed -b final_gtf_features/cds_unique.bed -wao -s | awk '$NF>299' | awk '{print$0,($17/($3-$2))*100}' OFS="\t" > gibbon_scfr_all_atleast_300bp_cds_unique.bed
time bedtools intersect -a gibbon_scfr_all_atleast_300bp.bed -b final_gtf_features/pseudogene_unique.bed -wao -s | awk '$NF>299' | awk '{print$0,($17/($3-$2))*100}' OFS="\t" > gibbon_scfr_all_atleast_300bp_pseudogene_unique.bed
time bedtools intersect -a gibbon_scfr_all_atleast_300bp.bed -b final_gtf_features/RNA_genes_unique.bed -wao -s | awk '$NF>299' | awk '{print$0,($17/($3-$2))*100}' OFS="\t" > gibbon_scfr_all_atleast_300bp_RNA_genes_unique.bed
time bedtools intersect -a gibbon_scfr_all_atleast_300bp.bed -b final_gtf_features/only_introns_refined_unique.bed -wao -s | awk '$NF>299' | awk '{print$0,($15/($3-$2))*100}' OFS="\t" > gibbon_scfr_all_atleast_300bp_only_introns_refined_unique.bed
time bedtools intersect -a gibbon_scfr_all_atleast_300bp.bed -b final_gtf_features/utr_refined_unique.bed -wao -s | awk '$NF>299' | awk '{print$0,($15/($3-$2))*100}' OFS="\t" > gibbon_scfr_all_atleast_300bp_utr_refined_unique.bed
time bedtools intersect -a gibbon_scfr_all_atleast_300bp.bed -b final_gtf_features/only_intergenic_unique.bed -wao | awk '$NF>299' | awk '{print$0,($11/($3-$2))*100}' OFS="\t" > gibbon_scfr_all_atleast_300bp_only_intergenic_unique.bed


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run DFT

time bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/gibbon/GCA_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna -bed gibbon_scfr_all_atleast_300bp_cds_unique.bed -name+ -s > gibbon_scfr_all_atleast_300bp_cds_unique.fa

#create conda environment
conda create -n scfr python=3.10 numpy scipy biopython matplotlib tqdm -y
conda activate scfr

#21m35.241s
time python3 /media/aswin/SCFR/SCFR-main/Fourier_analysis/scfr_parallel_fft_motif_report_grouped.py -o gibbon_scfr_all_atleast_300bp_cds_unique -t 32 gibbon_scfr_all_atleast_300bp_cds_unique.fa
#
time python3 /media/aswin/SCFR/SCFR-main/Fourier_analysis/scfr_fourier_chromosome_wise_summary.py gibbon_scfr_all_atleast_300bp_cds_unique --top 3 --cores 32

#Plot density graph
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/plot_fourier_frequencies.R summary.tsv summary.pdf all
cd /media/aswin/SCFR/SCFR-main/genes/gibbon/test/dft_test/chromosome_wise_summary
Rscript plot_fourier_frequencies2.R summary.tsv summary.pdf

#PLot heatmap
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/plot_heatmap_fourier_peaks_genes.R all_species_all_genes_positive_freq_mag.tsv all_species_all_genes_positive_freq_mag.pdf

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#create conda environment
conda create -n scfr python=3.10 numpy scipy biopython matplotlib tqdm pandas -y
conda activate scfr

for scfr in $(ls gibbon_scfr_all_atleast_300bp_*.bed | sed 's/.bed//g')
do
echo ">"$scfr

#Get fasta
time bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/gibbon/GCA_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna -bed $scfr".bed" -name+ -s > $scfr".fa"


#Run DFT
#22m4.830s
time python3 /media/aswin/SCFR/SCFR-main/Fourier_analysis/scfr_parallel_fft_motif_report_grouped.py -o $scfr -t 32 $scfr".fa"
#Get DFT summary
time python3 /media/aswin/SCFR/SCFR-main/Fourier_analysis/scfr_fourier_chromosome_wise_summary.py $scfr --top 3 --cores 32

#Plot density graph
Rscript /media/aswin/SCFR/SCFR-main/genes/gibbon/test/revision_plot_fourier_frequencies.R $scfr/chromosome_wise_summary/summary.tsv $scfr"_frequency_density.pdf"

#Plot heatmap
Rscript /media/aswin/SCFR/SCFR-main/genes/gibbon/test/revision_plot_heatmap_fourier_peaks_genes.R $scfr/chromosome_wise_summary/summary.tsv $scfr"_frequency_histogram.pdf"
done


