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
#912m3.824s
time python3 /media/aswin/programs/CPC2_standalone-1.0.1/bin/CPC2.py -i human_scfr_all_lesser_than_300bp.fa -o cpc2_human_scfr_all_lesser_than_300bp


#SCFRs of atleast 300bp
#0m38.741s
time bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed ../human_scfr_all_atleast_300bp.bed -s -name+ > human_scfr_all_atleast_300bp.fa
#18m54.841s
time python3 /media/aswin/programs/CPC2_standalone-1.0.1/bin/CPC2.py -i human_scfr_all_atleast_300bp.fa -o cpc2_human_scfr_all_atleast_300bp

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Convert CPC2 output text to gff to be compatible to view more details when clicked 
time python3 ./cpc2_to_gff3.py cpc2_human_scfr_all_atleast_300bp.txt
#87m7.050s
time python3 ./cpc2_to_gff3.py cpc2_human_scfr_all_lesser_than_300bp.txt

#If it's slow to load, index it:
grep -v '^#track' cpc2_human_scfr_all_atleast_300bp.gff3 | sort -k1,1 -k4,4n | bgzip > cpc2_human_scfr_all_atleast_300bp.gff3.gz
tabix -p gff cpc2_human_scfr_all_atleast_300bp.gff3.gz

#99m34.365s
time grep -v '^#track' cpc2_human_scfr_all_lesser_than_300bp.gff3 | sort -k1,1 -k4,4n | bgzip > cpc2_human_scfr_all_lesser_than_300bp.gff3.gz
#6m7.863s
time tabix -p gff cpc2_human_scfr_all_lesser_than_300bp.gff3.gz

scp  cpc2_human_scfr_all_atleast_300bp.gff3.gz cpc2_human_scfr_all_atleast_300bp.gff3.gz.tbi cpc2_human_scfr_all_lesser_than_300bp.gff3.gz cpc2_human_scfr_all_lesser_than_300bp.gff3.gz.tbi ceglab8@172.28.65.118:~/Downloads/SCFR/

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Annotate an SCFR BED file with CPC2 coding/noncoding status.
time python3 annotate_scfr_cpc2.py \
  human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology_with_0.33_dft_results_no_overlap_with_gene_and_ntblastn_hits.bed \
  Coding_potential/cpc2_human_scfr_all_atleast_300bp.txt \
  -o scfr_with_cpc2_status.tsv \
  --unmatched scfr_no_cpc2_match.bed

#SCFRs with coding status
grep -w coding scfr_with_cpc2_status.tsv | awk '{print$1,$2,$3,$7,1,$6}' OFS="\t"

#Convert the 84 SCFRs with DFT & CPC2 results into ucsc compatible bed 
python3 scfr_with_cpc2_status_to_ucsc_custom_track.py scfr_with_cpc2_status.tsv -o scfr_with_cpc2_status_ucsc.bed

#
bedtools getfasta -fi /media/aswin/SCFR/SCFR-main/genomes/human/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna -bed <(grep -w coding scfr_with_cpc2_status.tsv | awk '{print$1,$2,$3,$7,1,$6}' OFS="\t") -s -name+ 

bedtools intersect -a <(grep -w coding scfr_with_cpc2_status.tsv | awk '{print$1,$2,$3,$7,1,$6}' OFS="\t") -b <(grep -v "^chrM" Repetitive_Elements.bed)


#Observations of 3 SCFRs with coding potential:
>NC_060926.1	68462816	68463848	3	1	+	SCFR_NC_060926_1_68462816_68463848_frame3_fft	5;3;0.444767;0.148256;0.333333;483.194;167.722;147.049;2;7;3	frame+coords	3::NC_060926.1:68462816-68463848(+)	1032	281	0.37781000000000003	9.460429191589355	1	0.99932	coding
	-For ucsc visualization: chr2:68,456,185-68,482,903
	-This region is shared only with chimpanzee & pygmy chim chain, absent in other primates
	-Contain simple repeats, a LINE elment is also inserted 
	-These features are absent in nearby genes

>NC_060927.1	109352916	109353879	-1	1	-	SCFR_NC_060927_1_109352916_109353879_frame_1_fft	9;3;0.444444;0.333333;0.389408;313.356;206.221;188.980;2;3;3	frame+coords	-1::NC_060927.1:109352916-109353879(-)	963	180	0.36944000000000005	10.11052837371826	1	0.888501	coding	
	-Only chimp shared this region that too very fragmented, not even pygmy chimp share this region
	-Structural variants present tis region
	-Simple repeats present in this region 
	-These features are absent in nearby genes

>NC_060929.1	1761824	1762922	3	1	+	SCFR_NC_060929_1_1761824_1762922_frame3_fft	6;3;0.493625;0.495446;0.333333;405.225;353.389;201.383;2;2;3	frame+coords	3::NC_060929.1:1761824-1762922(+)	1098	325	0.34168	7.750473976135253	10.9999	coding	
	-Shared with chimp, pygmy chimp & even gorilla
	-This region is unique to the T2T-CHM13 v2.0 assembly compared to the GRCh38/hg38 and GRCh37/hg19 reference assemblies.

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


awk 'NR > 1 { count[$NF]++; total++ } END { for (val in count) printf "%s: %d (%.2f%%)\n", val, count[val], (count[val]/total)*100 }' cpc2_human_scfr_all_atleast_300bp.txt



#Merge coding & intergenic status
awk -F'\t' -v OFS='\t' -v tag="intergenic" -v other="genic" '
  NR==FNR { key[$4 "::" $1 ":" $2 "-" $3 "(" $6 ")"] = 1; next }
  /^#/    { print $0, "Region_type"; next }
  { print $0, ($1 in key ? tag : other) }
' ../filtering_intergenic_SCFR/human_scfr_all_atleast_300bp_only_intergenic_unique.bed cpc2_human_scfr_all_atleast_300bp.txt > cpc2_human_scfr_all_atleast_300bp_with_intergenic_info.txt

awk -F'\t' -v OFS='\t' -v tag="No_homology" -v other="homology" '
  NR==FNR { key[$4 "::" $1 ":" $2 "-" $3 "(" $6 ")"] = 1; next }
  /^#/    { print $0, "homology_status"; next }
  { print $0, ($1 in key ? tag : other) }
' ../filtering_intergenic_SCFR/human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed cpc2_human_scfr_all_atleast_300bp.txt > cpc2_human_scfr_all_atleast_300bp_with_intergenic_homology_info.txt

#Merge coding & homology status
awk -F'\t' -v OFS='\t' -v tag="No_homology" -v other="homology" '
  NR==FNR { key[$4 "::" $1 ":" $2 "-" $3 "(" $6 ")"] = 1; next }
  /^#/    { print $0, "Homology_status"; next }
  { print $0, ($1 in key ? tag : other) }
' ../filtering_intergenic_SCFR/human_scfr_all_atleast_300bp_only_intergenic_unique_with_no_homology.bed cpc2_human_scfr_all_atleast_300bp_with_intergenic_info.txt > cpc2_human_scfr_all_atleast_300bp_with_intergenic_homology_info.txt


#Count amount of homology
./crosstab.sh -f cpc2_human_scfr_all_atleast_300bp_with_intergenic_info.txt
./crosstab.sh -f cpc2_human_scfr_all_atleast_300bp_with_intergenic_homology_info.txt 
awk -F "\t" '{print$(NF-2),$(NF-1), $NF}' cpc2_human_scfr_all_atleast_300bp_with_intergenic_homology_info.txt | sort | uniq -c

awk 'NR > 1 { count[$NF]++; total++ } END { for (val in count) printf "%s: %d (%.2f%%)\n", val, count[val], (count[val]/total)*100 }' cpc2_human_scfr_all_atleast_300bp_with_homology.txt
./crosstab.sh -f cpc2_human_scfr_all_atleast_300bp_with_homology.txt -r -2 -c -1




