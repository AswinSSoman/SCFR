#!/bin/bash

cd /media/aswin/SCFR/SCFR-main
species=$1

echo ">"$species
mkdir -p exon_shadow/"$species"
cd exon_shadow/"$species"

#Prepare inputs:

#filter & classiffy exon-scfr overlaps & calculate exon shadow
#create output files
echo "chr start end frame filler strand chr first_exon_start last_exon_end gene transcript exon_strand exon_frame exon_count upstream_len_in_scfr downstream_len_in_scfr first_exon_number last_exon_number gene_tx_count first_exon_tx_count last_exon_tx_count first_exon_order last_exon_order first_exon_sharing last_exon_sharing first_exon_splicing last_exon_splicing first_exon_overlap last_exon_overlap merged_txs" | tr " " "\t" > "$species"_multi_exon.tsv
echo "chr start end frame filler strand chr first_exon_start last_exon_end gene transcript exon_strand exon_frame exon_count upstream_len_in_scfr downstream_len_in_scfr first_exon_number last_exon_number gene_tx_count first_exon_tx_count last_exon_tx_count first_exon_order last_exon_order first_exon_sharing last_exon_sharing first_exon_splicing last_exon_splicing first_exon_overlap last_exon_overlap merged_txs" | tr " " "\t" > "$species"_composite_exon.tsv

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

if [[ $mc > 1 ]]
then

#save multi-exon containing SCFRs (multiple exons per SCFR or single transcript per SCFR)
  for ut in $(awk '{print$11}' scfr_temp.bed | sort -u)
  do
  awk -v ut="$ut" '$11==ut' scfr_temp.bed | sort -k8,8n -k9,9n > temp_multi_sorted
  #plus strand
  if [[ $strand == "+" ]]
  then
  e1=$(head -1 temp_multi_sorted | awk '{print$1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}')
  e2=$(tail -1 temp_multi_sorted| awk '{print$9,$14,$16,$17,$18,$19,$20,$22}')
  ec=$(wc -l < temp_multi_sorted)
  paste <(echo "$e1") <(echo "$e2") | awk -v ec="$ec" '{print$1,$2,$3,$4,$5,$6,$7,$8,$21,$9,$10,$11,$12,ec,$20,$28,$13,$22,$14,$15,$23,$16,$24,$17,$25,$18,$26,$19,$27}'
  #minus strand
  elif [[ $strand == "-" ]]
  then
  e1=$(head -1 temp_multi_sorted | awk '{print$1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$22}')
  e2=$(tail -1 temp_multi_sorted | awk '{print$9,$14,$16,$17,$18,$19,$20,$21}')
  ec=$(wc -l < temp_multi_sorted)
  paste <(echo "$e1") <(echo "$e2") | awk -v ec="$ec" '{print$1,$2,$3,$4,$5,$6,$7,$8,$21,$9,$10,$11,$12,ec,$28,$20,$22,$13,$14,$23,$15,$24,$16,$25,$17,$26,$18,$27,$19}'
  else :
  fi
  rm temp_multi_sorted
  unset e1 e2 ec
  #Filter transcripts with identical exons, keep only one
  done > temp_multi_or_composite

#Check total exon count for all transcripts (count single exon transcripts: exon count less than 2)
  tec=$(awk '{print$11}' scfr_temp.bed | sort -V | uniq -c | awk '$1<2' | wc -l)

#classify multi Vs composite
  if [[ $tec -lt 1 ]]
  then
  cat temp_multi_or_composite | awk '{ key=$7 FS $8 FS $9 FS $10 FS $12 FS $13 FS $15 FS $16 FS $22 FS $23 FS $26 FS $27 FS $28 FS $29; c[key]++; if(!(key in f)){f[key]=$0;o[++n]=key} } END{for(i=1;i<=n;i++) print f[o[i]],c[o[i]]}' | tr " " "\t" >> "$species"_multi_exon.tsv
  else
  cat temp_multi_or_composite | awk '{ key=$7 FS $8 FS $9 FS $10 FS $12 FS $13 FS $15 FS $16 FS $22 FS $23 FS $26 FS $27 FS $28 FS $29; c[key]++; if(!(key in f)){f[key]=$0;o[++n]=key} } END{for(i=1;i<=n;i++) print f[o[i]],c[o[i]]}' | tr " " "\t" >> "$species"_composite_exon.tsv
  fi
  rm temp_multi_or_composite
  unset tec

else :
fi

#delete temp files & unset variables
rm scfr_temp.bed
unset tnoe mc ut
done < <(awk '{print$1,$2,$3,$4,$5,$6}' OFS="\t" "$species"_scfr_containing_cds.bed | sort -u)

unset scfr
cd /media/aswin/SCFR/SCFR-main

