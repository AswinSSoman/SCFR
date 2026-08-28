######################################################################################################################################################################################################################################################################################################
#CODE TO RUN DFT ON ALL SCFRs FROM ALL SPECIES 
######################################################################################################################################################################################################################################################################################################


######################################################################################################################################################################################################################################################################################################
#211m55.810s for gibbon without intergenic region + 154m44.505s

#Time taken: 354m32.204s; gibbon: 26m46.632s, bonobo: 31m36.784s, chimpanzee: 29m40.073s ,gorilla: 37m36.156s,borangutan: 34m54.705s,sorangutan: 34m56.297s
cd /media/aswin/SCFR/SCFR-main
#time for species in human gibbon bonobo chimpanzee gorilla borangutan sorangutan 
time for species in gibbon bonobo chimpanzee gorilla borangutan sorangutan 
	do
	echo ">"$species

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
#Create folder
mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/$species/300
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/$species/300

#Set varaiables
pscfr=$(find /media/aswin/SCFR/SCFR-main/exon_shadow/ -mindepth 2 -maxdepth 2 -name "*_scfr_all.bed" | grep "$species")
pgtf=$(find /media/aswin/SCFR/SCFR-main/genes/ -mindepth 2 -maxdepth 2 -name "*.gtf" | grep "$species")
xscfr=$(basename "$pscfr" .bed)
xgtf=$(basename "$pgtf" .gtf)

######################################################################################################################################################################################################################################################################################################
#Filter SCFRs

#3m43.931s
time awk -F "\t" '($3-$2)>=300' $pscfr > $xscfr"_atleast_300bp.bed"

######################################################################################################################################################################################################################################################################################################
#Filter GTF

#Get relevant info from gtf in bed format
#1m58.695s
time awk -f /media/aswin/SCFR/SCFR-main/Fourier_analysis/revision_gtf_to_bed.awk $pgtf > $xgtf".bed"

#Extract features
awk -F "\t" '$4=="CDS"' $xgtf".bed" > cds.bed
awk '!seen[$1,$2,$3,$4,$5,$6]++' cds.bed > cds_unique.bed

awk -F "\t" '$7=="pseudogene"' $xgtf".bed" > pseudogene.bed
awk '!seen[$1,$2,$3,$4,$5,$6]++' pseudogene.bed > pseudogene_unique.bed 

awk -F "\t" '$7~"RNA"' $xgtf".bed" > RNA_genes.bed
awk '!seen[$1,$2,$3,$4,$5,$6]++' RNA_genes.bed > RNA_genes_unique.bed

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get introns & UTRs bed from GTF

#Install tool agat
#conda create -n agat_env -c bioconda -c conda-forge agat
unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
conda activate agat_env

#Extract all features including introns & UTRs
#18m36.006s
time agat_sp_add_introns.pl --gff $pgtf --out introns.bed --cpu 28
#keep only introns
awk -F "\t" '$3=="intron"' introns.bed > only_introns.bed

#Keep relevant columns
awk 'BEGIN{FS="\t"; OFS="\t"}
{
    gene_id="."; Parent=".";
    n=split($9, a, ";");
    for(i=1; i<=n; i++){
        # Clean up potential spaces in GFF attributes
        gsub(/^ /, "", a[i]); 
        split(a[i], b, "=");
        if(b[1]=="gene_id") gene_id=b[2];
        if(b[1]=="Parent") Parent=b[2];
    }
    # Standard BED6 + extra columns
    # $1=chrom, $4=start, $5=end, $3=name, $6=score(placeholder), $7=strand
    if($4 ~ /^[0-9]+$/) print $1, $4, $5, $3"_"gene_id, "1", $7, gene_id, Parent
}' only_introns.bed > only_introns_refined.bed

#remove duplicates
awk '!seen[$1,$2,$3,$4,$5,$6]++' only_introns_refined.bed | awk '{print$1,$2,$3,$4"_"$8,$5,$6}' | tr " " "\t" > only_introns_refined_unique.bed

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
awk '!seen[$1,$2,$3,$4,$5,$6]++' utr_refined.bed | awk '{print$1,$2,$3,$4,$5,$6}' | tr " " "\t" > utr_refined_unique.bed

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get intergenic regions

#8m23.934s, 16m6.086s
time agat_sp_add_intergenic_regions.pl --gff $pgtf --out intergenic.bed --cpu 28
awk -F "\t" '$3~"intergenic_region"' intergenic.bed > only_intergenic.bed
#Keep relevant columns
awk -F "\t" '!/^#/ {print$1,$4,$5,$9}' only_intergenic.bed | awk '!seen[$1,$2,$3,$4]++' | tr " " "\t" > only_intergenic_unique.bed

#Simply check
#total genome coverage by gene
awk -F "\t" '$4=="gene"' $xgtf".bed" | awk '!seen[$1,$2,$3,$6]++' | awk '{sum+=($3-$2-1);} END{print sum}'
#total genome coverage by intergenic
awk '{sum+=($3-$2);} END{print sum}' only_intergenic_unique.bed

#Get intergenic regions in antisense strand
awk -F "\t" '$4=="gene"' $xgtf".bed" > genes.bed

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Collect all features

mkdir final_gtf_features
mv cds_unique.bed pseudogene_unique.bed RNA_genes_unique.bed only_introns_refined_unique.bed utr_refined_unique.bed only_intergenic_unique.bed genes.bed final_gtf_features/

#Find strand aware overlaps (minimum 300bp)

#with cds (last column tells % of total SCFR length covered by overlap)
time bedtools intersect -a $xscfr"_atleast_300bp.bed" -b final_gtf_features/cds_unique.bed -wao -s | awk '$NF>299' | awk '{print$0,($17/($3-$2))*100}' OFS="\t" > $xscfr"_atleast_300bp_cds_unique.bed"
time bedtools intersect -a $xscfr"_atleast_300bp.bed" -b final_gtf_features/pseudogene_unique.bed -wao -s | awk '$NF>299' | awk '{print$0,($17/($3-$2))*100}' OFS="\t" > $xscfr"_atleast_300bp_pseudogene_unique.bed"
time bedtools intersect -a $xscfr"_atleast_300bp.bed" -b final_gtf_features/RNA_genes_unique.bed -wao -s | awk '$NF>299' | awk '{print$0,($17/($3-$2))*100}' OFS="\t" > $xscfr"_atleast_300bp_RNA_genes_unique.bed"
time bedtools intersect -a $xscfr"_atleast_300bp.bed" -b final_gtf_features/only_introns_refined_unique.bed -wao -s | awk '$NF>299' | awk '{print$0,($15/($3-$2))*100}' OFS="\t" > $xscfr"_atleast_300bp_only_introns_refined_unique.bed"
time bedtools intersect -a $xscfr"_atleast_300bp.bed" -b final_gtf_features/utr_refined_unique.bed -wao -s | awk '$NF>299' | awk '{print$0,($15/($3-$2))*100}' OFS="\t" > $xscfr"_atleast_300bp_utr_refined_unique.bed"
#time bedtools intersect -a $xscfr"_atleast_300bp.bed" -b final_gtf_features/only_intergenic_unique.bed -wao | awk '$NF>299' | awk '{print$0,($11/($3-$2))*100}' OFS="\t" > $xscfr"_atleast_300bp_only_intergenic_unique.bed"
time bedtools intersect -a $xscfr"_atleast_300bp.bed" -b final_gtf_features/only_intergenic_unique.bed -wao | awk '$NF>299' | awk '{print$0,($11/($3-$2))*100}' OFS="\t" \
| awk 'BEGIN{OFS="\t"} 
     $1 != "" && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2 < $3 {
         print $1,$2,$3,$4,$5,$6
     }' > $xscfr"_atleast_300bp_only_intergenic_unique.bed"
time bedtools intersect -a $xscfr"_atleast_300bp.bed" -b final_gtf_features/genes.bed -wao | awk '$NF>299' | awk -F "\t" '$6!=$12' | awk '{print$0,($17/($3-$2))*100}' OFS="\t" >  $xscfr"_atleast_300bp_antisense_intergenic.bed"

######################################################################################################################################################################################################################################################################################################

#create conda environment
#conda create -n scfr python=3.10 numpy scipy biopython matplotlib tqdm pandas -y
conda activate scfr

time for scfr in $(ls $xscfr"_atleast_300bp_"*.bed | sed 's/.bed//g')
do
echo ">"$scfr

#Get fasta
genome=$(find /media/aswin/SCFR/SCFR-main/genomes/$species/ -mindepth 1 -maxdepth 1 -name "GC*.fna")
#0m12.317s
time bedtools getfasta -fi $genome -bed $scfr".bed" -name+ -s > $scfr".fa"

#Run DFT
#22m4.830s
time python3 /media/aswin/SCFR/SCFR-main/Fourier_analysis/scfr_parallel_fft_motif_report_grouped_novisuals.py -o $scfr -t 32 $scfr".fa"
#Get DFT summary
time python3 /media/aswin/SCFR/SCFR-main/Fourier_analysis/scfr_fourier_chromosome_wise_summary.py $scfr --top 3 --cores 32

#Plot density graph
gt=$(ls $scfr/chromosome_wise_summary/summary.tsv | cut -f1 -d "/" | cut -f6- -d "_" | sed 's/unique//g' | sed 's/only//g' | sed 's/_$//g' | sed 's/^_//g')
Rscript /media/aswin/SCFR/SCFR-main/Fourier_analysis/revision_plot_fourier_frequencies.R $scfr/chromosome_wise_summary/summary.tsv $scfr"_frequency_density.pdf" $gt

#Plot heatmap
Rscript /media/aswin/SCFR/SCFR-main/Fourier_analysis/revision_plot_heatmap_fourier_peaks_genes.R $scfr/chromosome_wise_summary/summary.tsv $scfr"_frequency_histogram.pdf" $gt
unset genome gt
done
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Unset variables
unset pscfr pgtf xscfr xgtf

cd /media/aswin/SCFR/SCFR-main

done

######################################################################################################################################################################################################################################################################################################
#Collect all figures

mkdir /media/aswin/SCFR/SCFR-main/Fourier_analysis/collect_figures
cd /media/aswin/SCFR/SCFR-main/Fourier_analysis/collect_figures

find . -mindepth 3 -maxdepth 3 -name "*frequency_density.pdf" -type f | xargs -n1 sh -c 'cp $0 /media/aswin/SCFR/SCFR-main/Fourier_analysis/collect_figures/'
find . -mindepth 3 -maxdepth 3 -name "*frequency_histogram.pdf" -type f | xargs -n1 sh -c 'cp $0 /media/aswin/SCFR/SCFR-main/Fourier_analysis/collect_figures/'

#Display all results in an html file
python3 make_scfr_gallery.py --input /media/aswin/SCFR/SCFR-main/Fourier_analysis/collect_figures --output /media/aswin/SCFR/SCFR-main/Fourier_analysis/scfr_gallery --dpi 100



