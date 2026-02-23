#######################################################################################################################################################################################################################################################################################################
#PCA
#######################################################################################################################################################################################################################################################################################################


#######################################################################################################################################################################################################################################################################################################
#Length bin:

#Run PCA & k-mean clustering for different length bins & see if clustering gets better with SCFR length bins (189m8.519s)
process_species() {
  species="$1"

  echo ">$species"

  SCFR="SCFR_all/${species}_SCFR_all.out"
  GENOME_FASTA="genomes/${species}/GC*.fna"
  BASEOUT="/media/aswin/SCFR/SCFR-main/Length_bin_PCA_kmeans/${species}"

  mkdir -p "$BASEOUT"

  awk '
    BEGIN {
      OFS = "\t"

      lo["1000_2500"]=1000;   hi["1000_2500"]=2500
      lo["2500_5000"]=2500;   hi["2500_5000"]=5000
      lo["5000_7500"]=5000;   hi["5000_7500"]=7500
      lo["7500_10000"]=7500;  hi["7500_10000"]=10000

      # NEW BIN (>10000)
      lo["gt10000"]=10001;    hi["gt10000"]=1e12
    }
    {
      len = $3 - $2
      strand = ($4 ~ "-") ? "-" : "+"

      for (b in lo) {
        if (len >= lo[b] && len <= hi[b]) {
          print $1, $2, $3, $4, 1, strand > "'"$BASEOUT"'/" b ".bed"
          break
        }
      }
    }
  ' "$SCFR"

  for lenbin in 1000_2500 2500_5000 5000_7500 7500_10000 gt10000
  do
    OUTDIR="$BASEOUT/$lenbin"
    mkdir -p "$OUTDIR"

    BED="$BASEOUT/$lenbin.bed"
    FASTA="$OUTDIR/${species}_${lenbin}.fasta"

    [[ ! -s "$BED" ]] && continue

    sort -k1,1 -k2,2n "$BED" \
    | bedtools getfasta \
        -fi $GENOME_FASTA \
        -bed stdin \
        -name+ -s \
    > "$FASTA"
python3 /media/aswin/SCFR/SCFR-main/my_scripts/PCA/corrected_rscu_calc.py "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/pca_script.R "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/kmeans_script.R "$OUTDIR" "$OUTDIR"
  done
}
export -f process_species

#Run for 7 species
#time parallel -j 7 process_species ::: human bonobo chimpanzee gorilla borangutan sorangutan gibbon
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon; do
    process_species "$species"
done

#Check run
cd /media/aswin/SCFR/SCFR-main/Length_bin_PCA_kmeans
for i in $(ls -d */ | tr -d "/")
do
for f in 1000_2500 2500_5000 5000_7500 7500_10000 gt10000
do
pdf=$(find $i/$f -mindepth 1 -maxdepth 1 -name "pca_cluster_PC1_PC2.pdf" -type f)
if [[ -z $pdf ]]; then p="not-run"; else p="run"; fi
ns=$(grep ">" -c $i/$f/$i"_"$f".fasta")
ts=$(grep -v ">" $i/$f/$i"_"$f".fasta" | wc | awk '{print$3-$1}')
echo $i $f $p $ns $ts
unset pdf p ns ts
done
unset f
done | column -t

#######################################################################################################################################################################################################################################################################################################
#Length threshold: (243m58.707s)

process_species_threshold() {
  species="$1"

  echo ">$species"

  SCFR="SCFR_all/${species}_SCFR_all.out"
  GENOME_FASTA="genomes/${species}/GC*.fna"
  BASEOUT="/media/aswin/SCFR/SCFR-main/Length_threshold_PCA_kmeans/${species}"

  mkdir -p "$BASEOUT"

  awk '
    BEGIN {
      OFS = "\t"
      thr["gt1000"]=1000
      thr["gt2500"]=2500
      thr["gt5000"]=5000
      thr["gt7500"]=7500
      thr["gt10000"]=10000
    }
    {
      len = $3 - $2
      strand = ($4 ~ "-") ? "-" : "+"

      for (t in thr) {
        if (len > thr[t]) {
          print $1, $2, $3, $4, 1, strand > "'"$BASEOUT"'/" t ".bed"
        }
      }
    }
  ' "$SCFR"

  for threshold in gt1000 gt2500 gt5000 gt7500 gt10000
  do
    OUTDIR="$BASEOUT/$threshold"
    mkdir -p "$OUTDIR"

    BED="$BASEOUT/$threshold.bed"
    FASTA="$OUTDIR/${species}_${threshold}.fasta"

    [[ ! -s "$BED" ]] && continue

    sort -k1,1 -k2,2n "$BED" \
    | bedtools getfasta \
        -fi $GENOME_FASTA \
        -bed stdin \
        -name+ -s \
    > "$FASTA"
python3 /media/aswin/SCFR/SCFR-main/my_scripts/PCA/corrected_rscu_calc.py "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/pca_script.R "$OUTDIR" "$OUTDIR"
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/kmeans_script.R "$OUTDIR" "$OUTDIR"
  done
}
export -f process_species_threshold

#Run for 7 species
#time parallel -j 7 process_species_threshold ::: human bonobo chimpanzee gorilla borangutan sorangutan gibbon
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon; do
    process_species_threshold "$species"
done

#Check run
cd /media/aswin/SCFR/SCFR-main/Length_threshold_PCA_kmeans
for i in $(ls -d */ | tr -d "/")
do
for f in gt1000 gt2500 gt5000 gt7500 gt10000
do
pdf=$(find $i/$f -mindepth 1 -maxdepth 1 -name "pca_cluster_PC1_PC2.pdf" -type f)
if [[ -z $pdf ]]; then p="not-run"; else p="run"; fi
ns=$(grep ">" -c $i/$f/$i"_"$f".fasta")
ts=$(grep -v ">" $i/$f/$i"_"$f".fasta" | wc | awk '{print$3-$1}')
echo $i $f $p $ns $ts
unset pdf p ns ts
done
unset f
done | column -t

#######################################################################################################################################################################################################################################################################################################
#label coding vs non-coding SCFRs in PCA

cd /media/aswin/SCFR/SCFR-main
for bed in $(find Length_threshold_PCA_kmeans/ -mindepth 2 -maxdepth 2 -name "gt1000.bed" -type f | xargs readlink -f)
do
species=$(echo $bed | cut -f7 -d "/")
sp=$(grep -i "$species" /media/aswin/SCFR/SCFR-main/genome_reports/species_names | awk '{print$2}' | tr "_" " ")
path=$(echo $bed | sed 's!/gt1000.bed!!g')
cd $path
echo ">"$species
#scfr overlapping exons in same frame
bedtools intersect -a $bed -b /media/aswin/SCFR/SCFR-main/exon_shadow/$species/"$species"_coding_exons.bed -wo | awk '$12$13==$4' | awk '{$NF = $NF + 1; print}' | awk '{print$0,$9-$8+1,$3-$2+1,($20/($9-$8+1))*100}' \
| awk '{print$4"::"$1":"$2"-"$3"("$6")","in_frame_coding",$20,$21,$22,$23}' | sort -k1,1 -k6,6rn | awk '!seen[$1]++' | sed '1i SCFR coding_status overlap_len cds_len scfr_len percent_coding_in_scfr' | sed 's/[ ]\+/\t/g' > scfr_inframe_coding_gt1000.tsv
#scfr overlapping exons in diff frame
bedtools intersect -a $bed -b /media/aswin/SCFR/SCFR-main/exon_shadow/$species/"$species"_coding_exons.bed -wo | awk '$12$13!=$4' | awk '{$NF = $NF + 1; print}' | awk '{print$0,$9-$8+1,$3-$2+1,($20/($9-$8+1))*100}' \
| awk '{print$4"::"$1":"$2"-"$3"("$6")","out_frame_coding",$20,$21,$22,$23}' | sort -k1,1 -k6,6rn | awk '!seen[$1]++' | sed '1i SCFR coding_status overlap_len cds_len scfr_len percent_coding_in_scfr' | sed 's/[ ]\+/\t/g' > scfr_outframe_coding_gt1000.tsv
#scfr_inframe_coding_gt1000.tsv & scfr_outframe_coding_gt1000.tsv have some common scfrs, remove these common ones from scfr_outframe_coding_gt1000.tsv
awk -F'\t' 'NR==FNR {seen[$1]; next} !($1 in seen)' scfr_inframe_coding_gt1000.tsv scfr_outframe_coding_gt1000.tsv | sed '1i SCFR coding_status overlap_len cds_len scfr_len percent_coding_in_scfr' | sed 's/[ ]\+/\t/g'  > scfr_outframe_coding_gt1000_filtered.tsv
#scfr not overlapping exons
awk '{print $1}' scfr_inframe_coding_gt1000.tsv scfr_outframe_coding_gt1000_filtered.tsv > exclude_list.txt
seqkit grep -v -f exclude_list.txt gt1000/"$species"_gt1000.fasta | grep ">" | sed 's/$/ non_coding/g' | sed '1i SCFR coding_status' | sed 's/[ ]\+/\t/g' > scfr_non_coding_gt1000.tsv
#final list of scfrs & coding status
awk '{print$1,$2}' scfr_inframe_coding_gt1000.tsv scfr_outframe_coding_gt1000_filtered.tsv scfr_non_coding_gt1000.tsv | grep -v "coding_status" | sed '1i SCFR coding_status' | sed 's/[ ]\+/\t/g' > scfr_coding_status.tsv
unset species sp path 
cd /media/aswin/SCFR/SCFR-main
done

#######################################################################################################################################################################################################################################################################################################
#Identify repeats in SCFRs of atleast 1000bp length

#Run repeat masker (ran for 3 species (chimpanzee, gorilla, human) using species library, for rest hominidiae is used as species)
cd /media/aswin/SCFR/SCFR-main
time for seq in $(find Length_threshold_PCA_kmeans/ -mindepth 3 -name "*_gt1000.fasta" -type f | xargs readlink -f)
do
mkdir -p repeat_masker/"$species"
species=$(echo $seq | cut -f7 -d "/")
sp=$(grep -i "$species" /media/aswin/SCFR/SCFR-main/genome_reports/species_names | awk '{print$2}' | tr "_" " ")
echo ">"$species
cd repeat_masker
time /media/aswin/programs/RepeatMasker-4.1.7/RepeatMasker/RepeatMasker -pa 8 -a -s -u -gff -html -species "Hominidae" -dir $species -xsmall $seq
unset species sp
cd /media/aswin/SCFR/SCFR-main
done

#######################################################################################################################################################################################################################################################################################################
#Get composition

cd /media/aswin/SCFR/SCFR-main
for seq in $(find Length_threshold_PCA_kmeans/ -mindepth 3 -maxdepth 3 -name "*gt1000.fasta" -type f | xargs readlink -f)
do
species=$(echo $seq | cut -f7 -d "/")
sp=$(grep -i "$species" /media/aswin/SCFR/SCFR-main/genome_reports/species_names | awk '{print$2}' | tr "_" " ")
path=$(echo $seq | awk -F "/" '!($NF="")'  OFS="/")
cd $path
echo ">"$species
/media/aswin/SCFR/SCFR-main/my_scripts/exon_shadow/get_composition.sh $seq > "$species"_gt1000_composition.tsv
cd /media/aswin/SCFR/SCFR-main
done

#######################################################################################################################################################################################################################################################################################################
#Plot PCA clusters & color by different factors

#For length bin
cd /media/aswin/SCFR/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
echo ">"$species
cd Length_bin_PCA_kmeans/"$species"/
coding_status=$(readlink -f /media/aswin/SCFR/SCFR-main/Length_threshold_PCA_kmeans/"$species"/scfr_coding_status.tsv)
repeat_class=$(readlink -f /media/aswin/SCFR/SCFR-main/repeat_masker/"$species"/"$species"_gt1000.fasta.out)
python3 /media/aswin/SCFR/SCFR-main/my_scripts/PCA/fastaout_to_tsv.py $repeat_class /media/aswin/SCFR/SCFR-main/repeat_masker/"$species"/"$species"_gt1000_pca_input.out
repeat_input=$(readlink -f /media/aswin/SCFR/SCFR-main/repeat_masker/"$species"/"$species"_gt1000_pca_input.out)
gc=$(readlink -f /media/aswin/SCFR/SCFR-main/Length_threshold_PCA_kmeans/"$species"/gt1000/"$species"_gt1000_composition.tsv)
for lenbin in 2500_5000 5000_7500 7500_10000 gt10000
do
echo " -"$lenbin
cd $lenbin
#plot PCA by color
if [[ "$lenbin" == "2500_5000" ]]; then point="3"; elif [[ "$lenbin" == "5000_7500" ]]; then point="6"; elif [[ "$lenbin" == "7500_10000" ]] || [[ "$lenbin" == "gt10000" ]]; then point="7"; else :; fi
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_color_by_clade.R . . $species 14 10 $point $lenbin
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_cluster.R . . $species 14 10 $point $lenbin
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_coding_status.R . . $coding_status $species 14 10 $point $lenbin
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_repeat_family.R . . $repeat_input $species 14 10 $point $lenbin
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_GC_content.R . . $gc $species 14 10 $point $lenbin
unset point
cd ../
done
unset lenbin coding_status repeat_class repeat_input
cd /media/aswin/SCFR/SCFR-main/
done

#For length threshold
cd /media/aswin/SCFR/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
echo ">"$species
cd Length_threshold_PCA_kmeans/"$species"/
coding_status=$(readlink -f /media/aswin/SCFR/SCFR-main/Length_threshold_PCA_kmeans/"$species"/scfr_coding_status.tsv)
repeat_class=$(readlink -f /media/aswin/SCFR/SCFR-main/repeat_masker/"$species"/"$species"_gt1000.fasta.out)
python3 /media/aswin/SCFR/SCFR-main/my_scripts/PCA/fastaout_to_tsv.py $repeat_class /media/aswin/SCFR/SCFR-main/repeat_masker/"$species"/"$species"_gt1000_pca_input.out
repeat_input=$(readlink -f /media/aswin/SCFR/SCFR-main/repeat_masker/"$species"/"$species"_gt1000_pca_input.out)
gc=$(readlink -f /media/aswin/SCFR/SCFR-main/Length_threshold_PCA_kmeans/"$species"/gt1000/"$species"_gt1000_composition.tsv)
for lenbin in gt2500 gt5000 gt7500 gt10000
do
echo " -"$lenbin
cd $lenbin
#plot PCA by color
if [[ "$lenbin" == "gt2500" ]]; then point="3"; elif [[ "$lenbin" == "gt5000" ]]; then point="6"; elif [[ "$lenbin" == "gt7500" ]] || [[ "$lenbin" == "gt10000" ]]; then point="7"; else :; fi
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_color_by_clade.R . . $species 14 10 $point $lenbin
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_cluster.R . . $species 14 10 $point $lenbin
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_coding_status.R . . $coding_status $species 14 10 $point $lenbin
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_repeat_family.R . . $repeat_input $species 14 10 $point $lenbin
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_GC_content.R . . $gc $species 14 10 $point $lenbin
unset point
cd ../
done
unset lenbin coding_status repeat_class repeat_input
cd /media/aswin/SCFR/SCFR-main/
done

#######################################################################################################################################################################################################################################################################################################
#Plot PCA loadings

#Visuzalize codons contribution & amino acids
cd /media/aswin/SCFR/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
echo ">"$species
cd Length_threshold_PCA_kmeans/"$species"/
for lenbin in gt2500 gt5000 gt7500 gt10000
do
echo " -"$lenbin
cd $lenbin
len=$(echo $lenbin | sed 's/gt//g')
#Plor heatmap of codon loadings
#Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_codon_heatmap.R pca_loadings.tsv "$species"_"$lenbin"_codon_loadings_heatmap.pdf "$species" "$len"
#Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_codon_heatmap_updated.R pca_loadings.tsv "$species"_"$lenbin"_top_codon_loadings_heatmap.pdf "$species" "$len" 2 5 explained_variance.tsv
#Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_codon_heatmap_updated.R pca_loadings.tsv "$species"_"$lenbin"_top_codon_loadings_heatmap.pdf "$species" "$len" 5 10  explained_variance.tsv
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_codon_contribution_heatmap.R pca_loadings.tsv "$species"_"$lenbin"_top_codon_loadings_heatmap.pdf "$species" "$len" 2 5 explained_variance.tsv
#Plot variance explained by top 10 PCs
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_variance_explained_stat.R explained_variance.tsv "$species"_"$lenbin"_top_10_PC_variance.pdf "$species" "$len"
unset len
cd ../
done
unset lenbin
cd /media/aswin/SCFR/SCFR-main/
done

cd /media/aswin/SCFR/SCFR-main
mkdir shreya/Length_threshold_PCA_kmeans/top_codon_loadings
mkdir shreya/Length_threshold_PCA_kmeans/top_PC_variance
find Length_threshold_PCA_kmeans/ -name "*_top_codon_loadings_heatmap.pdf" -type f | grep -v "test" | xargs -n1 sh -c 'cp $0 /media/aswin/SCFR/SCFR-main/shreya/Length_threshold_PCA_kmeans/top_codon_loadings/'
find Length_threshold_PCA_kmeans/ -name "*_top_10_PC_variance.pdf" -type f | xargs -n1 sh -c 'cp $0 /media/aswin/SCFR/SCFR-main/shreya/Length_threshold_PCA_kmeans/top_PC_variance/'

#######################################################################################################################################################################################################################################################################################################
#Summary of PCA at different lengths

#Summary of Length bin
cd /media/aswin/SCFR/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
cd Length_bin_PCA_kmeans/"$species"/
for lenbin in 2500_5000 5000_7500 7500_10000 gt10000
do
cd $lenbin
i1=$(awk '$1 == "PC1" || $1 == "PC2" {sum += $3} END {print sum}' explained_variance.tsv)
i2=$(awk 'NR>1 && $5 != "NA" {if($5 > max) {max=$5; k=$0}} END {print k}' k_optimization_scores.tsv)
echo $species $lenbin $i1 $i2
unset i1 i2
cd ../
done
unset lenbin
cd /media/aswin/SCFR/SCFR-main/
done | sed '1i Species Length_bin PC1_PC2 k Silhouette DBI WCSS Curvature' | sed 's/[ ]\+/\t/g' > /media/aswin/SCFR/SCFR-main/Length_bin_PCA_kmeans/all_species_pca_clustering_summary.tsv

#Summary of Length threshold
cd /media/aswin/SCFR/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
cd Length_threshold_PCA_kmeans/"$species"/
for lenthr in gt2500 gt5000 gt7500 gt10000
do
cd $lenthr
i1=$(awk '$1 == "PC1" || $1 == "PC2" {sum += $3} END {print sum}' explained_variance.tsv)
i2=$(awk 'NR>1 && $5 != "NA" {if($5 > max) {max=$5; k=$0}} END {print k}' k_optimization_scores.tsv)
echo $species $lenthr $i1 $i2
unset i1 i2
cd ../
done
unset lenthr
cd /media/aswin/SCFR/SCFR-main/
done | sed '1i Species Length_threshold PC1_PC2 k Silhouette DBI WCSS Curvature' | sed 's/[ ]\+/\t/g' > /media/aswin/SCFR/SCFR-main/Length_threshold_PCA_kmeans/all_species_pca_clustering_summary.tsv

#Plot PCA comparison
cd /media/aswin/SCFR/SCFR-main
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_pca_clustering_summary.R Length_bin_PCA_kmeans/all_species_pca_clustering_summary.tsv Length_bin_PCA_kmeans/all_species_pca_clustering_summary.png
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_pca_clustering_summary.R Length_threshold_PCA_kmeans/all_species_pca_clustering_summary.tsv Length_threshold_PCA_kmeans/all_species_pca_clustering_summary.png


#######################################################################################################################################################################################################################################################################################################
#Run PCA for genes
#######################################################################################################################################################################################################################################################################################################

cd /media/aswin/SCFR/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon
do
gene=$(readlink -f /media/aswin/SCFR/SCFR-main/Fourier_analysis/genes/"$species"/GCF_*_cds.fa)
echo $species $gene
mkdir -p PCA_genes/$species/
cd PCA_genes/$species/
#Filter genes (min length: 2500)
time python3 /media/aswin/SCFR/SCFR-main/my_scripts/PCA/filter_genes_for_PCA.py $gene "$species"_cds_longest_per_gene.fa 2500
#Calculate composition
/media/aswin/SCFR/SCFR-main/my_scripts/exon_shadow/get_composition.sh "$species"_cds_longest_per_gene.fa > "$species"_cds_longest_per_gene_composition.tsv
#Run PCA
time python3 /media/aswin/SCFR/SCFR-main/my_scripts/PCA/corrected_rscu_calc.py . .
time Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/pca_script.R . .
time Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/kmeans_script.R . .
#plot
time Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_cluster.R . . $species 14 10 2 genes
time Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_by_GC_content.R . . "$species"_cds_longest_per_gene_composition.tsv $species 14 10 2 genes
cd /media/aswin/SCFR/SCFR-main/
done

#Plot PCA loadings of genes

#Visuzalize codons contribution & amino acids
cd /media/aswin/SCFR/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
echo ">"$species
cd PCA_genes/$species/
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/PCA/plot_codon_contribution_heatmap.R pca_loadings.tsv "$species"_genes_codon_loadings_heatmap.pdf "$species"
cd /media/aswin/SCFR/SCFR-main/
done



#######################################################################################################################################################################################################################################################################################################
#DRAFT
#######################################################################################################################################################################################################################################################################################################

#copy & collect pdf in one location for slide making
cd /media/aswin/SCFR/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
cd Length_threshold_PCA_kmeans/"$species"
for lenthr in gt2500 gt5000 gt7500 gt10000
do
cd $lenthr
mkdir -p /media/aswin/SCFR/SCFR-main/shreya/Length_threshold_PCA_kmeans/$species/
for pdf in $(ls pca_*_PC1_PC2.pdf)
do
pc=$(echo $pdf | cut -f3,4 -d "_" | sed 's/\.pdf//g')
col=$(echo $pdf | cut -f2 -d "_" | sed 's/\.pdf//g')
cp $pdf /media/aswin/SCFR/SCFR-main/shreya/Length_threshold_PCA_kmeans/$species/"$species"_"$lenthr"_"$pc"_"$col".pdf
unset pc col
done
cp pca_loadings.tsv /media/aswin/SCFR/SCFR-main/shreya/Length_threshold_PCA_kmeans/$species/"$species"_"$lenthr"_pca_loadings.tsv
cd ../
done
unset lenthr
cd /media/aswin/SCFR/SCFR-main/
done



#Code to view the pca codon raw loadings & contribution
awk '{print $1, $2, $3}' pca_loadings.tsv | awk 'NR==1 {print $0, "Sum_Abs"} NR>1 {abs2=($2<0?-$2:$2); abs3=($3<0?-$3:$3); print $0, abs2+abs3}' | awk 'NR==1 {print $0, "PC1_sq", "PC2_sq"} NR>1 {print $0, $2*$2, $3*$3}' | awk '
  NR == 1 {header = $0 " Norm_PC1_sq Norm_PC2_sq"; next}
  {
    # Store each line in an array and accumulate sums
    rows[NR] = $0
    val5[NR] = $5
    val6[NR] = $6
    sum5 += $5
    sum6 += $6}
  END {print header
    for (i = 2; i <= NR; i++) {
      # Divide each stored row by the final sums
      print rows[i], (val5[i]/sum5)*100, (val6[i]/sum6)*100}}' | column -t

find . -maxdepth 3 -mindepth 3  -name "*.fasta" -type f | xargs -n1 bash -c 'paste <(echo $0 | cut -f4 -d "/") <(grep -v ">" $0 | wc | awk "{print\$3-\$1}")' > total_fasta_length

