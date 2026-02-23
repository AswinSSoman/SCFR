#!/bin/bash

cd /media/aswin/SCFR/SCFR-main
species=$1

echo ">"$species
mkdir -p exon_shadow/"$species"
cd exon_shadow/"$species"

#Prepare inputs:
  #get zero based scfr bed file
#  awk '{if($4~"-") print$0,"1","-"; else print$0,"1","+"}' OFS="\t" /media/aswin/SCFR/SCFR-main/SCFR_all/"$species"_SCFR_all.out > "$species"_scfr_all.bed
  #Convert gtf to zero based cds bed (headers: chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
#  gtf=$(find /media/aswin/SCFR/SCFR-main/genes/"$species" -name "GCF*.gtf")
#  rm "$species"_coding_exons.bed
#  python3 /media/aswin/SCFR/SCFR-main/my_scripts/exon_shadow/gtf_to_cds_with_transcript_exon_metadata.py -i $gtf -s "$species" -o "$species"_coding_exons.bed
  #get SCFR exon overlaps (headers: chrom start end frame filler strand chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
#  bedtools intersect -a "$species"_scfr_all.bed -b "$species"_coding_exons.bed -wo -s > "$species"_scfr_cds_all_overlaps.bed
  #get scfrs containing exons with same frame as exons & calculate their upstream & downstream shadow (scfr is 0-based & cds is 1-based, hence add 1 when cheking overlap & subtract 1 when calculating shadow length)
#  awk '{c4=$4; gsub("-","",c4); if(c4==$13 && $2<=($8+1) && $3>=$9) {if($6=="+") print$0,$8-$2-1,$3-$9; else if($6=="-") print$0,$3-$9,$8-$2-1}}' "$species"_scfr_cds_all_overlaps.bed > "$species"_scfr_containing_cds.bed

#filter & classiffy exon-scfr overlaps & calculate exon shadow
#create output files
#mv "$species"_single_exon.tsv old_13jan_"$species"_single_exon.tsv
#mv "$species"_multi_exon.tsv old_13jan_"$species"_multi_exon.tsv
#mv "$species"_composite_exon.tsv old_13jan_"$species"_composite_exon.tsv
#mv "$species"_exitron_candidates.tsv old_13jan_"$species"_exitron_candidates.tsv
#  echo "chr start end frame filler strand chr exon_start exon_end gene transcript exon_strand exon_frame exon_number gene_tx_count exon_tx_count exon_order exon_sharing exon_splicing overlap_len upstream_len_in_scfr downstream_len_in_scfr merged_txs" | tr " " "\t" > "$species"_single_exon.tsv
#  echo "chr start end frame filler strand chr first_exon_start last_exon_end gene transcript exon_strand exon_frame exon_count upstream_len_in_scfr downstream_len_in_scfr first_exon_number last_exon_number gene_tx_count first_exon_tx_count last_exon_tx_count first_exon_order last_exon_order first_exon_sharing last_exon_sharing first_exon_splicing last_exon_splicing first_exon_overlap last_exon_overlap merged_txs" | tr " " "\t" > "$species"_multi_exon.tsv
#  echo "chr start end frame filler strand chr first_exon_start last_exon_end gene transcript exon_strand exon_frame exon_count upstream_len_in_scfr downstream_len_in_scfr first_exon_number last_exon_number gene_tx_count first_exon_tx_count last_exon_tx_count first_exon_order last_exon_order first_exon_sharing last_exon_sharing first_exon_splicing last_exon_splicing first_exon_overlap last_exon_overlap merged_txs" | tr " " "\t" > "$species"_composite_exon.tsv
  echo "chr start end frame filler frame chr exon_1_start exon_1_end exon_2_start exon_2_end gene transcript exon_strand exon_frame intron_start intron_end intron_length merged_txs" | tr " " "\t" > "$species"_exitron_candidates.tsv

#Loop over each unique SCFRs to calculate exon shadow
time while read scfr
do
grep "$scfr" "$species"_scfr_containing_cds.bed > scfr_temp.bed
#strand
strand=$(awk 'NR==1{print $6; exit}' scfr_temp.bed)
#Total number of exons
tnoe=$(wc -l < scfr_temp.bed)
#Number of overlapping exons
mc=$(awk '{print$7,$8,$9,$10,$11,$12}' OFS="\t" scfr_temp.bed | bedtools sort -i - | bedtools merge -i - -s -c 5 -o count | wc -l)

#Classify exon shadow
#  if [[ $tnoe == 1 ]] || [[ $mc == 1 ]]
#  then
  #save unique single exon containing SCFRs (single exon per SCFR)
  #filter identical exons, keep only one
#  awk '{ k=$7 FS $8 FS $9 FS $10 FS $12 FS $13 FS $17 FS $19 FS $20 FS $21 FS $22; c[k]++; if(!(k in f)){f[k]=$0;o[++n]=k} } END{for(i=1;i<=n;i++) print f[o[i]],c[o[i]]}' scfr_temp.bed | tr " " "\t" >> "$species"_single_exon.tsv
 if [[ $mc > 1 ]]
  then

#Check total exon count for all transcripts (count single exon transcripts: exon count less than 2)

#exitron calculation
  unset ut
  for ut in $(awk '{print$11}' scfr_temp.bed | sort -u)
  do
  if [[ $strand == "+" ]]
  then
  awk -v ut="$ut" '$11==ut' scfr_temp.bed | sort -k8,8n -k9,9n > temp_multi_sorted
 cex=$(awk '{print$14}' temp_multi_sorted | sort -n | paste -s -d " " | awk '{print$2-$1}')
 if [[ $cex == 1 ]]
 then
  awk '{print$1,$2,$3,$4,$5,$6,$7,a,b,$8,$9,$10,$11,$12,$13,b,$8,$8-b-1; a=$8;b=$9}' OFS="\t" temp_multi_sorted | sed '1d'
 elif [[ $cex -gt 1 ]]
 then :
 fi
  rm temp_multi_sorted
  elif [[ $strand == "-" ]]
  then
  awk -v ut="$ut" '$11==ut' scfr_temp.bed | sort -k8,8nr -k9,9nr > temp_multi_sorted
 cex=$(awk '{print$14}' temp_multi_sorted | sort -n | paste -s -d " " | awk '{print$2-$1}')
 if [[ $cex == 1 ]]
 then
  awk '{print$1,$2,$3,$4,$5,$6,$7,a,b,$8,$9,$10,$11,$12,$13,$9,a,a-$9-1; a=$8;b=$9}' OFS="\t" temp_multi_sorted | sed '1d'
 elif [[ $cex -gt 1 ]]
 then :
 fi
  rm temp_multi_sorted
  else :
  fi
unset cex
  done | awk '{ key=$7 FS $8 FS $9 FS $10 FS $11 FS $12 FS $14 FS $15 FS $16 FS $17 FS $18; count[key]++; if (!(key in first)) { first[key]=$0; order[++n]=key } } END { for (i=1; i<=n; i++) print first[order[i]], count[order[i]] }' | tr " " "\t" >> "$species"_exitron_candidates.tsv
unset ut
else :
fi

#delete temp files & unset variables
rm scfr_temp.bed
unset tnoe mc
done < <(awk '{print$1,$2,$3,$4,$5,$6}' OFS="\t" "$species"_scfr_containing_cds.bed | sort -u)

unset scfr
cd /media/aswin/SCFR/SCFR-main
