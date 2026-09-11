
#Coding potential


#About:
#Coding Potential Calculator distinguishes protein-coding from non-coding RNAs based on the sequence features of the input transcripts.
#Last version of CPC1 is widely used by worldwide researchers. For better serve for scientific community, we update CPC1 to CPC2.
#It can discriminate the coding and non-coding transcripts faster and more accurately. 
#Input:
#CPC2 accepts RNA transcript sequences as input. Both fasta format and GTF/GFF/BED format are supported.


cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300
#4m20.793s
time awk -F "\t" '($3-$2)<300' /media/aswin/SCFR/SCFR-main/exon_shadow/human/human_scfr_all.bed > human_scfr_all_lesser_than_300bp.bed

awk -F "\t" '$3-$2>=300' all_scfrs_overlapping_cllu1.bed > all_scfrs_overlapping_cllu1_atleast_300.bed
awk -F "\t" '$3-$2<300' all_scfrs_overlapping_cllu1.bed > all_scfrs_overlapping_cllu1_less_than_300.bed

cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300
bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed all_scfrs_overlapping_cllu1_atleast_300.bed -name+ -s | less > all_scfrs_overlapping_cllu1_atleast_300.fa
time python3 /media/aswin/programs/CPC2_standalone-1.0.1/bin/CPC2.py -i all_scfrs_overlapping_cllu1_atleast_300.fa -o cpc2_all_scfrs_overlapping_cllu1_atleast_300


NC_060936.1:92,374,600-92,444,426
bedtools intersect -a /media/aswin/SCFR/SCFR-main/exon_shadow/human/human_scfr_all.bed -b <(echo -e "NC_060936.1\t92374600\t92444426\tCLLU1") -wa 

bedtools intersect -a /media/aswin/SCFR/SCFR-main/exon_shadow/human/human_scfr_all.bed -b <(echo -e "NC_060936.1\t92374600\t92444426\tCLLU1") -wa -u > all_scfrs_overlapping_cllu1.bed


cd /media/aswin/SCFR/SCFR-main/genomes/human/test
python3 /media/aswin/programs/CPC2_standalone-1.0.1/bin/CPC2.py -i ACTA1_gene.fa -o ACTA1_gene_cpc2.out



#
#De Novo Origin of Human Protein-Coding Genes (Dong-Dong Wu, November 10, 2011) [https://doi.org/10.1371/journal.pgen.1002379]
#This papers tells:
# In 2009, Knowles and McLysaght identified three putative protein coding genes: CLLU1, c22orf45, and DNAH10OS, which had a de novo origin in the human genome.
# These genes were identified by employing a straightforward, but rigorous, procedure which provided transcriptional and translational evidence, and allowed them to estimate that about 0.075% of the human protein coding genes may have originated de novo from noncoding regions
 
# But databases like ncbi ensembl shows this as long non-coding RNA.
# In NCBI, this gene is annotated in many primates, even in a shrew.
 
 
# Expression of this gene has been shown to be upregulated in some individuals with chronic lymphocytic leukemia (CLL), and has been used for prognostic and diagnostic purposes.
# This gene was originally identified as a human-specific putative protein-coding gene due to the presence of a peptide (PAp00140670, HIIYSTFLSK) that could have supported translation at this locus. 
# This peptide is not present in more recent builds of PeptideAtlas, and the presence of a protein product at this locus has not been independently verified. 
# For this reason, this gene is being represented as non-coding. Sequence comparisons to other primates indicates that no other primate is predicted to contain an open reading frame. [provided by RefSeq, Feb 2017]


#Potential to code for short peptide

#(CLLU1 as an emerging biomarker in chronic lymphoid leukemia (Chunmeng Rong, 20 March 2024)  https://link.springer.com/article/10.1007/s13577-024-01051-4)

#Previous studies considered CLLU1 to be a non-coding RNA; however, recent research has discovered that its coding sequence region possesses the potential to encode a short peptide similar to interleukin-4.
#Remarkably, abnormally elevated expression of CLLU1 has only been detected in chronic lymphoid leukemia among all hematological cancers. 
#High CLLU1 expression often indicates more malignant pathological features and an unfavorable prognosis for patients.

#CLLU1 is located on chromosome 12q22 and comprises three exons, flanked by BTG1 and EEA1.
#It encodes six mRNA transcripts that do not exhibit sequence homology with any known genes.
#Among these transcripts, CLLU1-203 and the coding sequence (CDS) display the highest expression levels.
#The majority of these transcripts cluster on chromosome 12q22, with most being non-coding, while a few, such as cDNA 4 and 5, potentially encode a peptide similar to interleukin-4 (IL-4).
#The CDS likely encodes a short peptide chain consisting of 121 amino acids

Protein Binding to Cis-Motifs in mRNAs Coding Sequence Is Common and Regulates Transcript Stability and the Rate of Translation [https://doi.org/10.3390/cells10112910]


tblastn \
  -db /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna \
  -query CLLU1_protein.fa \
  -num_threads $(nproc) \
  -max_target_seqs 10 \
  -evalue 1e-5 \
  -outfmt "6 sseqid sstart send qseqid bitscore sstrand" | \
awk 'BEGIN {OFS="\t"} {
  # Assign strand (+ or -)
  strand = ($6 == "minus") ? "-" : "+";
  
  # Ensure sstart is smaller than send for start/end coordinates
  if ($2 < $3) {
    start = $2 - 1; # Convert 1-based BLAST to 0-based BED start
    end = $3;
  } else {
    start = $3 - 1;
    end = $2;
  }
  
  print $1, start, end, $4, $5, strand;
}' > tblastn_hits.bed



bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed tblastn_hits.bed -s -name+ > CLLU1_orf.fa



mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Coding_potential
time python3 /media/aswin/programs/CPC2_standalone-1.0.1/bin/CPC2.py -i human_scfr_all_atleast_300bp_cds_unique.fa -o Coding_potential/human_scfr_all_atleast_300bp_cds_unique_cpc2

cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Coding_potential
awk '{if($NF=="noncoding") print$1}' human_scfr_all_atleast_300bp_cds_unique_cpc2.txt | cut -f3- -d ":" | tr ":-" "\t" | grep -v "^#" | tr -d "()+-" > cds_scfrs_with_cpc2_non_coding_status.bed
awk '{if($NF=="coding") print$1}' human_scfr_all_atleast_300bp_cds_unique_cpc2.txt | cut -f3- -d ":" | tr ":-" "\t" | grep -v "^#" | tr -d "()+-" > cds_scfrs_with_cpc2_coding_status.bed

bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed ../cds.bed -s -name+ > cds.fa
time python3 /media/aswin/programs/CPC2_standalone-1.0.1/bin/CPC2.py -i cds.fa -o cpc2_cds

























