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
#44m9.003s
time bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed human_scfr_all_lesser_than_300bp.bed -s -name+ > human_scfr_all_lesser_than_300bp.fa
#
time python3 /media/aswin/programs/CPC2_standalone-1.0.1/bin/CPC2.py -i human_scfr_all_lesser_than_300bp.fa -o cpc2_human_scfr_all_lesser_than_300bp


#SCFRs of atleast 300bp
#0m38.741s
time bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed ../human_scfr_all_atleast_300bp.bed -s -name+ > human_scfr_all_atleast_300bp.fa
#18m54.841s
time python3 /media/aswin/programs/CPC2_standalone-1.0.1/bin/CPC2.py -i human_scfr_all_atleast_300bp.fa -o cpc2_human_scfr_all_atleast_300bp

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Convert CPC2 output text to gff to be compatible to view more details when clicked 
time python3 ./cpc2_to_gff3.py cpc2_human_scfr_all_atleast_300bp.txt 

#If it's slow to load, index it:
grep -v '^#track' cpc2_human_scfr_all_atleast_300bp.gff3 | sort -k1,1 -k4,4n | bgzip > cpc2_human_scfr_all_atleast_300bp.gff3.gz
tabix -p gff cpc2_human_scfr_all_atleast_300bp.gff3.gz

scp cpc2_human_scfr_all_atleast_300bp.gff3.gz cpc2_human_scfr_all_atleast_300bp.gff3.gz.tbi ceglab8@172.28.65.118:~/Downloads/SCFR/
