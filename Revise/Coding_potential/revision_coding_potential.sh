######################################################################################################################################################################################################################################################################################################
#																																																																																																																																				CHECK CODING POTENTIAL OF SCFRS
######################################################################################################################################################################################################################################################################################################


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#CPC2:
		#Coding Potential Calculator distinguishes protein-coding from non-coding RNAs based on the sequence features of the input transcripts.
		#Last version of CPC1 is widely used by worldwide researchers. For better serve for scientific community, we update CPC1 to CPC2.
		#It can discriminate the coding and non-coding transcripts faster and more accurately. 
		#Input:
		#CPC2 accepts RNA transcript sequences as input. Both fasta format and GTF/GFF/BED format are supported.

#NOTE:
#CPC2 works best for spliced RNA/transcript/CDS not on individual exons. Hence if a de novo gene is single exonic it might work but otherwise, CPC2 won't work.

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Coding_potential
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/human/300/Coding_potential

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run CPC2 on all SCFRs

#SCFRs shorter than 300bp
#4m20.793s
time awk -F "\t" '($3-$2)<300' /media/aswin/SCFR/SCFR-main/exon_shadow/human/human_scfr_all.bed > human_scfr_all_lesser_than_300bp.bed
#
time bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed human_scfr_all_lesser_than_300bp.bed -s -name+ > human_scfr_all_lesser_than_300bp.fa
#
python3 /media/aswin/programs/CPC2_standalone-1.0.1/bin/CPC2.py -i human_scfr_all_lesser_than_300bp.fa -o cpc2_human_scfr_all_lesser_than_300bp


#SCFRs of atleast 300bp
#0m38.741s
time bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed ../human_scfr_all_atleast_300bp.bed -s -name+ > human_scfr_all_atleast_300bp.fa
#
time python3 /media/aswin/programs/CPC2_standalone-1.0.1/bin/CPC2.py -i human_scfr_all_atleast_300bp.fa -o cpc2_human_scfr_all_atleast_300bp

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




