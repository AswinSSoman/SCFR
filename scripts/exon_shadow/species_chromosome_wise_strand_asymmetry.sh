#!/bin/bash

#3 positional arguments: species window-size output-folder

species=$1

mkdir -p $3

#🔹 Step 1: Make fixed genomic windows
#bedtools makewindows -g /media/aswin/SCFR/SCFR-main/genome_sizes/borangutan.genome -w 100000 > borangutan.windows.bed
bedtools makewindows -g /media/aswin/SCFR/SCFR-main/genome_sizes/"$species".genome -w $2 > $3/"$species"_windows.bed

#🔹 Step 2: Intersect BED with windows
awk 'NR>1{print$1,$2,$3,$4,$5,$6}' OFS="\t" "$species"_single_exon.tsv > "$species"_single_exon_scfr.bed
bedtools intersect -a $3/"$species"_windows.bed -b "$species"_single_exon_scfr.bed -wa -wb > $3/window_hits.bed
rm "$species"_single_exon_scfr.bed

#🔹 Step 3: Count + and − per window
awk '{key = $1 FS $2 FS $3
  if ($NF == "+") plus[key]++
  else if ($NF == "-") minus[key]++}
END { for (k in plus) {
    p = plus[k]
    m = minus[k] + 0
    if (p + m > 0)
      print k, p, m, (p-m)/(p+m)}
  for (k in minus) if (!(k in plus)) {
    p = 0; m = minus[k]
    print k, p, m, (p-m)/(p+m)}}' $3/window_hits.bed > $3/window_SA.tsv

#🔹 Step 4: Extract highly asymmetric windows
#awk '$6 > 0.5 {print}' window_asymmetry.tsv > neg_strand_hotspots.bed
awk '$6>=0.6 && ($4+$5)>=10' $3/window_SA.tsv > $3/pos_hotspots.bed
awk '$6<=-0.6 && ($4+$5)>=10' $3/window_SA.tsv > $3/neg_hotspots.bed

#Plot
cd $3/
Rscript /media/aswin/SCFR/SCFR-main/my_scripts/exon_shadow/plot_SA_by_chromosome.R window_SA.tsv "$species"_
cd ../
