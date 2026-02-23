
cd /media/aswin/SCFR/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
mkdir exon_shadow/"$species"
cd exon_shadow/"$species"

#Prepare inputs:
  #get zero based scfr bed file
  awk '{if($4~"-") print$0,"1","-"; else print$0,"1","+"}' OFS="\t" /media/aswin/SCFR/SCFR-main/SCFR_all/"$species"_SCFR_all.out > "$species"_scfr_all.bed
  #Convert gtf to zero based cds bed (headers: chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
  gtf=$(find /media/aswin/SCFR/SCFR-main/genes/"$species" -name "GCF*.gtf")
  python3 /media/aswin/SCFR/SCFR-main/my_scripts/exon_shadow/gtf_to_cds_with_transcript_exon_metadata.py -i $gtf -s "$species" -o "$species"_coding_exons.bed
  #get SCFR exon overlaps (headers: chrom start end frame filler strand chrom start end gene tx strand frame exon_number gene_tx_count exon_tx_count exon_order sharing splicing)
  bedtools intersect -a "$species"_scfr_all.bed -b "$species"_coding_exons.bed -wo -s > "$species"_scfr_cds_all_overlaps.bed
  #get scfrs containing exons with same frame as exons & calculate their upstream & downstream shadow
  awk '{c4=$4; gsub("-","",c4); if(c4==$13 && $2<=$8 && $3>=$9) {if($6=="+") print$0,$8-$2,$3-$9; else if($6=="-") print$0,$3-$9,$8-$2}}' "$species"_scfr_cds_all_overlaps.bed > "$species"_scfr_containing_cds.bed

#filter & classiffy exon-scfr overlaps & calculate exon shadow
#create output files
  echo "chr start end frame filler strand chr exon_start exon_end gene transcript exon_strand exon_frame exon_number gene_tx_count exon_tx_count exon_order exon_sharing exon_splicing overlap_len upstream_len_in_scfr downstream_len_in_scfr merged_txs" | tr " " "\t" > "$species"_single_exon.tsv
  echo "chr start end frame filler strand chr first_exon_start last_exon_end gene transcript exon_strand exon_frame exon_count upstream_len_in_scfr downstream_len_in_scfr first_exon_number last_exon_number gene_tx_count first_exon_tx_count last_exon_tx_count first_exon_order last_exon_order first_exon_sharing last_exon_sharing first_exon_splicing last_exon_splicing first_exon_overlap last_exon_overlap merged_txs" | tr " " "\t" > "$species"_multi_exon.tsv
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
  if [[ $tnoe == 1 ]] || [[ $mc == 1 ]]
  then
  #save unique single exon containing SCFRs (single exon per SCFR)
  #filter identical exons, keep only one
  awk '{ k=$7 FS $8 FS $9 FS $10 FS $12 FS $13 FS $17 FS $19 FS $20 FS $21 FS $22; c[k]++; if(!(k in f)){f[k]=$0;o[++n]=k} } END{for(i=1;i<=n;i++) print f[o[i]],c[o[i]]}' scfr_temp.bed | tr " " "\t" >> "$species"_single_exon.tsv
  elif [[ $mc > 1 ]]
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
  done | awk '{ key=$7 FS $8 FS $9 FS $10 FS $12 FS $13 FS $15 FS $16 FS $22 FS $23 FS $26 FS $27 FS $28 FS $29; c[key]++; if(!(key in f)){f[key]=$0;o[++n]=key} } END{for(i=1;i<=n;i++) print f[o[i]],c[o[i]]}' | tr " " "\t" >> "$species"_multi_exon.tsv

#exitron calculation
  unset ut
  for ut in $(awk '{print$11}' scfr_temp.bed | sort -u)
  do
  if [[ $strand == "+" ]]
  then
  awk -v ut="$ut" '$11==ut' scfr_temp.bed | sort -k8,8n -k9,9n > temp_multi_sorted
  awk '{print$1,$2,$3,$4,$5,$6,$7,a,b,$8,$9,$10,$11,$12,$13,b,$8,$8-b; a=$8;b=$9}' OFS="\t" temp_multi_sorted | sed '1d'
  rm temp_multi_sorted
  elif [[ $strand == "-" ]]
  then
  awk -v ut="$ut" '$11==ut' scfr_temp.bed | sort -k8,8nr -k9,9nr > temp_multi_sorted
  awk '{print$1,$2,$3,$4,$5,$6,$7,a,b,$8,$9,$10,$11,$12,$13,$9,a,a-$9; a=$8;b=$9}' OFS="\t" temp_multi_sorted | sed '1d'
  rm temp_multi_sorted
  else :
  fi
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
done

