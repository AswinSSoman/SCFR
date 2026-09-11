
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



