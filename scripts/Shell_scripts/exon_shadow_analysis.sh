####################################################################################################################################################################################################################################################################################################################
#Exon shadow

#set base directory
#base_directory=""

#Calculate exon shadow & exitron
mkdir $base_directory/SCFR-main/exon_shadow
cd $base_directory/SCFR-main
start_time=$(date +%s)
for sp in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
(
echo ">"$sp
#$base_directory/SCFR-main/my_scripts/exon_shadow/species_wise_exon_shadow_exitron_shell_script.sh $sp
$base_directory/SCFR-main/my_scripts/exon_shadow/per_species_run_exon_shadow_exitron.sh $sp
) &
done
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#################################################################################################################################################################################################################################################################################################
#QC

#genes to exclude:
	cd $base_directory/SCFR-main/
	time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	gtf=$(readlink -f genes/"$species"/*.gtf)
	awk '$3="CDS"' $gtf | grep -v "^#" | awk -F ";" '{for (n=1;n<=NF;n++) if($n~/note/) print $n}' | sort | uniq -c > $base_directory/SCFR-main/genes/"$species"_note_gtf
	#grep 'unassigned_transcript' $gtf | awk '$3="CDS"' > $base_directory/SCFR-main/genes/"$species"_unassigned_transcripts_gtf
	awk '$3="CDS"' $gtf | grep -v "^#" | awk -F ";" '{for (n=1;n<=NF;n++) if($n~/unassigned_transcript/) print $0}' | sort | uniq -c > $base_directory/SCFR-main/genes/"$species"_unassigned_transcripts_gtf
	grep 'partial "true"' $gtf | awk '$3="CDS"' > $base_directory/SCFR-main/genes/"$species"_partial_annotation_gtf
	unset gtf
	done

#Filter shadow entries
	cd $base_directory/SCFR-main
	time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo ">"$species
	gtf=$(readlink -f genes/"$species"/*.gtf)
	#exons to exclude
	awk '$3=="CDS"' $gtf | egrep 'partial "true"|insert|delet|stop codon' | awk -F'\t|;' -v OFS='\t' '{tid = exon = "" 
	  for (i = 9; i <= NF; i++) {if ($i ~ /transcript_id/) tid = $i
	    if ($i ~ /exon_number/)  exon = $i
	    if ($i ~ /gene_id/) gid =$i}
	  if (tid != "" || exon != "" || gid != "")
	    print $1, $4, $5, tid, exon, gid, $7}' | sed 's/transcript_id //g' | sed 's/exon_number//g' | sed 's/gene_id//g' | tr -d ' "' > genes/"$species"/"$species"_exons_to_exclude.tsv
	#Since excluding individual exons is more complicated especially in multi-exon SCFRs, removing exons eventually remove whole gene, hence rather than remove entries by exons, remove by genes
	awk '{print$6}' genes/"$species"/"$species"_exons_to_exclude.tsv | sort -u > genes/"$species"/"$species"_unique_genes_to_exclude
	#Filter shadow files
	awk 'NR==FNR { exclude[$1]; next }
	  !($10 in exclude)' genes/"$species"/"$species"_unique_genes_to_exclude exon_shadow/"$species"/"$species"_single_exon.tsv | awk -F "\t" '!($17=="last" && $22>0)' > exon_shadow/"$species"/"$species"_single_exon_filtered.tsv
	awk 'NR==FNR { exclude[$1]; next }
	  !($10 in exclude)' genes/"$species"/"$species"_unique_genes_to_exclude exon_shadow/"$species"/"$species"_multi_exon.tsv | awk -F "\t" '!($23=="last" && $16>0)' > exon_shadow/"$species"/"$species"_multi_exon_filtered.tsv
	awk 'NR==FNR { exclude[$1]; next }
	  !($10 in exclude)' genes/"$species"/"$species"_unique_genes_to_exclude exon_shadow/"$species"/"$species"_composite_exon.tsv | awk -F "\t" '!($23=="last" && $16>0)' > exon_shadow/"$species"/"$species"_composite_exon_filtered.tsv
	awk 'NR==FNR { exclude[$1]; next }
	  !($12 in exclude)' genes/"$species"/"$species"_unique_genes_to_exclude exon_shadow/"$species"/"$species"_exitron_candidates.tsv > exon_shadow/"$species"/"$species"_exitron_candidates_filtered.tsv
	unset gtf
	done

#Summary of how many non-zero downstream shadow values are present for every species with every exon type
	cd $base_directory/SCFR-main
	time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do 
	a=$(awk '$17=="last" && $22>0' exon_shadow/"$species"/"$species"_single_exon_filtered.tsv | wc -l)
	b=$(awk '$23=="last" && $16>0' exon_shadow/"$species"/"$species"_multi_exon_filtered.tsv | wc -l)
	c=$(awk '$23=="last" && $16>0' exon_shadow/"$species"/"$species"_multi_exon_filtered.tsv | wc -l)
	echo $a $b $c | sed "s/^/$species /g"
	unset a b c
	done | sed '1i species single multi composite' | column -t > exon_shadow/non_zero_last_exon_ds_shadow

#Summary of filtered
	cd $base_directory/SCFR-main
	time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do 
	a=$(wc -l < exon_shadow/"$species"/"$species"_single_exon.tsv)
	b=$(wc -l < exon_shadow/"$species"/"$species"_single_exon_filtered.tsv)
	c=$(wc -l < exon_shadow/"$species"/"$species"_multi_exon.tsv)
	d=$(wc -l < exon_shadow/"$species"/"$species"_multi_exon_filtered.tsv)
	e=$(wc -l < exon_shadow/"$species"/"$species"_composite_exon.tsv)
	f=$(wc -l < exon_shadow/"$species"/"$species"_composite_exon_filtered.tsv)
	echo $a $b $c $d $e $f | awk '{print$1, $2, $1-$2, $3, $4, $3-$4, $5, $6, $5-$6}' | sed "s/^/$species /g"
	unset a b c d e f
	done | sed '1i species single single_filtered difference multi multi_filtered difference composite composite_filtered difference' | column -t > exon_shadow/exon_shadow_filtered_summary

#Combine all shadows (single, multi & composite)
	cd $base_directory/SCFR-main
	time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	head -1 exon_shadow/"$species"/"$species"_multi_exon_filtered.tsv > exon_shadow/"$species"/"$species"_filtered_combined.tsv
	awk '{print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,1,$21,$22,$14,$14,$15,$16,$16,$17,$17,$18,$18,$19,$19,$20,$20,$23}' exon_shadow/"$species"/"$species"_single_exon_filtered.tsv | grep -v "upstream_len_in_scfr" >> exon_shadow/"$species"/"$species"_filtered_combined.tsv
	grep -v "upstream_len_in_scfr" exon_shadow/"$species"/"$species"_multi_exon_filtered.tsv >> exon_shadow/"$species"/"$species"_filtered_combined.tsv
	grep -v "upstream_len_in_scfr" exon_shadow/"$species"/"$species"_composite_exon_filtered.tsv >> exon_shadow/"$species"/"$species"_filtered_combined.tsv
	sed 's/[\t ]\+/\t/g' exon_shadow/"$species"/"$species"_filtered_combined.tsv -i
	awk '$15>0 || $16>0' exon_shadow/"$species"/"$species"_filtered_combined.tsv > exon_shadow/"$species"/"$species"_filtered_combined_non_zero.tsv
	done

#################################################################################################################################################################################################################################################################################################
#Plot shadow distriubtion

#Make summary tables
cd $base_directory/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
cd exon_shadow/"$species"
up=$(awk 'NR>1{print$15}' "$species"_filtered_combined_non_zero.tsv | python3 $base_directory/SCFR-main/my_scripts/get_stats.py | egrep "Count|^Minimum|^Maximum|^Mean|^Median|^Q1|^Q3" | awk -F ":" '{print$NF}' | tr -d " ," | paste -s -d " " | sed "s/^/up /g")
dn=$(awk 'NR>1{print$16}' "$species"_filtered_combined_non_zero.tsv | python3 $base_directory/SCFR-main/my_scripts/get_stats.py | egrep "Count|^Minimum|^Maximum|^Mean|^Median|^Q1|^Q3" | awk -F ":" '{print$NF}' | tr -d " ," | paste -s -d " " | sed "s/^/down /g")
echo -e $species $up"\n"$species $dn
unset up dn
cd $base_directory/SCFR-main
done | sed '1i Species direction N min max mean Q1 median Q3' | sed 's/[ ]\+/\t/g' > exon_shadow/all_species_total_shadow_length_distribution.tsv

cd $base_directory/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
cd exon_shadow/"$species"
up=$(awk 'NR>1{print$15}' "$species"_filtered_combined_non_zero.tsv | awk '$1>0' | python3 $base_directory/SCFR-main/my_scripts/get_stats.py | egrep "Count|^Minimum|^Maximum|^Mean|^Median|^Q1|^Q3" | awk -F ":" '{print$NF}' | tr -d " ," | paste -s -d " " | sed "s/^/up /g")
dn=$(awk 'NR>1{print$16}' "$species"_filtered_combined_non_zero.tsv | awk '$1>0' | python3 $base_directory/SCFR-main/my_scripts/get_stats.py | egrep "Count|^Minimum|^Maximum|^Mean|^Median|^Q1|^Q3" | awk -F ":" '{print$NF}' | tr -d " ," | paste -s -d " " | sed "s/^/down /g")
echo -e $species $up"\n"$species $dn
unset up dn
cd $base_directory/SCFR-main
done | sed '1i Species direction N min max mean Q1 median Q3' | sed 's/[ ]\+/\t/g' > exon_shadow/all_species_positive_shadow_length_distribution.tsv

#Plot distribution
cd $base_directory/SCFR-main
Rscript $base_directory/SCFR-main/Aishwarya_dwivedi/all_species_overall_data.R exon_shadow/all_species_total_shadow_length_distribution.tsv exon_shadow/all_species_total_shadow_length_distribution.pdf
Rscript $base_directory/SCFR-main/scripts/plot_all_species_exon_shadow_distribution.R exon_shadow/all_species_positive_shadow_length_distribution.tsv exon_shadow/all_species_positive_shadow_length_distribution.pdf

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Plot proportions

cd $base_directory/SCFR-main/
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
cd exon_shadow/"$species"
i1=$(awk 'NR > 1 {c4[$4]; c5[$5]; count++} END {print length(c4),length(c5), count}' "$species"_coding_exons.bed)
i2=$(awk 'NR > 1 {u10[$10]; u11[$11]; sum14 += $14} END {print length(u10),length(u11),sum14}' "$species"_filtered_combined_non_zero.tsv)
echo $species $i1 $i2 | awk '{print$0,($5/$2)*100, ($6/$3)*100, ($7/$4)*100}' | awk '{print$0,100-$8,100-$9,100-$10}'
unset i1 i2
cd $base_directory/SCFR-main/
done | sed '1i Species total_genes total_tx total_exons shadow_genes shadow_tx shadow_exons pc_genes_with_shadow pc_txs_with_shadow pc_exons_with_shadow pc_genes_wo_shadow pc_txs_wo_shadow pc_exons_wo_shadow' | sed 's/[ ]\+/\t/g' > $base_directory/SCFR-main/exon_shadow/plot/shadow_proportion_in_genes.tsv

cp $base_directory/SCFR-main/SCFR_all/all_scfr_length_stats.tsv exon_shadow/plot/
sed 's/Sumatran orangutan/sorangutan/g' exon_shadow/plot/all_scfr_length_stats.tsv -i
sed 's/Bornean orangutan/borangutan/g' exon_shadow/plot/all_scfr_length_stats.tsv -i
sed 's/./\L&/' exon_shadow/plot/all_scfr_length_stats.tsv -i

#Percent of SCFR with & without shadow: Don't plot since the percent os 99% and 1%
cd $base_directory/SCFR-main/
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
i1=$(grep -i $species $base_directory/SCFR-main/exon_shadow/plot/all_scfr_length_stats.tsv | awk '{print$2}')
i2=$(awk 'END {print NR-1}' exon_shadow/"$species"/"$species"_filtered_combined_non_zero.tsv)
i3=$(awk 'END {print NR-1}' exon_shadow/"$species"/"$species"_exitron_candidates_filtered.tsv)
echo $species $i1 $i2 $i3 | awk '{print$0,($3/$2)*100, ($4/$2)*100}' | awk '{print$0,100-$5,100-$6}'
unset i1 i2 i3
done | sed '1i Species total_scfr shadoe_scfr exitron_scfr pc_scfr_with_shadow pc_scfr_with_exitron pc_scfr_wo_shadow pc_scfr_wo_exitron' | sed 's/[ ]\+/\t/g' > $base_directory/SCFR-main/exon_shadow/plot/shadow_proportion_in_scfrs.tsv

#Exitron proportion in genes & transcripts
cd $base_directory/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
cd exon_shadow/"$species"
i1=$(awk 'NR > 1 {c4[$4]; c5[$5]} END {print length(c4),length(c5)}' "$species"_coding_exons.bed)
i2=$(awk 'NR > 1 {u12[$12]; u13[$13]} END {print length(u12), length(u13)}' "$species"_exitron_candidates_filtered.tsv)
echo $species $i1 $i2 | awk '{pg=($4/$2)*100; pt=($5/$3)*100; print $0, pg, pt, 100-pg, 100-pt}'
unset i1 i2
cd $base_directory/SCFR-main/
done | sed '1i Species total_genes total_tx genes_with_exitron tx_with_exitron pc_genes_with_exitron pc_txs_with_exitron pc_genes_wo_exitron pc_txs_wo_exitron' | sed 's/[ ]\+/\t/g' > $base_directory/SCFR-main/exon_shadow/plot/exitron_percentage_in_genes_txs.tsv

#Proportion of shadow types
cd $base_directory/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
cd exon_shadow/"$species"
t=$(awk 'END {print (NR-1)*2}' "$species"_filtered_combined_non_zero.tsv)
s=$(awk 'NR > 1 {
    if ($15 > 0) pos++; else if ($15 < 0) neg++; else zero++;
    if ($16 > 0) pos++; else if ($16 < 0) neg++; else zero++;} 
END {print pos+0, neg+0, zero+0}' "$species"_filtered_combined_non_zero.tsv)
d=$(awk 'NR > 1 { if ($15 > 0) p15++; if ($16 > 0) p16++ } END { print p15+0, p16+0 }' "$species"_filtered_combined_non_zero.tsv)
echo $t $s $d | awk '{print$0,($2/$1)*100, ($3/$1)*100, ($4/$1)*100, ($5/($5+$6))*100, ($6/($5+$6))*100}' | sed "s/^/$species /g"
unset t s d
cd $base_directory/SCFR-main
done | sed '1i Species total_shadows positive_shadows negative_shadows zero_shadows up_shadows down_shadows pc_positive_shadows pc_negative_shadows pc_zero_shadows pc_up_shadows pc_down_shadows' | sed 's/[ ]\+/\t/g' > $base_directory/SCFR-main/exon_shadow/plot/shadow_types.tsv

#Shadow pair
cd $base_directory/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
cd exon_shadow/"$species"
i1=$(awk 'END{print NR-1}' "$species"_filtered_combined_non_zero.tsv)
i2=$(awk 'END{print NR-1}' "$species"_exitron_candidates_filtered.tsv)
i3=$(awk 'NR > 1 && $16 <= 0 && $15 > 0 {count++} END {print count+0}' "$species"_filtered_combined_non_zero.tsv)
i4=$(awk 'NR > 1 && $15 <= 0 && $16 > 0 {count++} END {print count+0}' "$species"_filtered_combined_non_zero.tsv)
i5=$(awk 'NR > 1 && $15 > 0 && $16 > 0 {count++} END {print count+0}' "$species"_filtered_combined_non_zero.tsv)
echo $i1 $i2 $i3 $i4 $i5 | awk '{print$0,($2/$1)*100,($3/$1)*100,($4/$1)*100,($5/$1)*100}' | awk '{print$1,$2,$3,$4,$5,$6,100-$6,$7,$8,$9}' | sed "s/^/$species /g"
unset i1 i2 i3 i4 i5
cd $base_directory/SCFR-main
done | sed '1i Species total_shadows exitron_shadows only_up_shadows only_down_shadows both_shadow pc_shadow_with_exitron pc_shadow_without_exitron pc_only_up pc_only_down pc_both_shadow' | sed 's/[ ]\+/\t/g' > $base_directory/SCFR-main/exon_shadow/plot/paired_shadow_proportions.tsv

#################################################################################################################################################################################################################################################################################################
#Strand asymmetry in total SCFR, CDS & shadow

mkdir $base_directory/SCFR-main/exon_shadow/asymmetry
cd $base_directory/SCFR-main/exon_shadow/asymmetry
#rm all_species_cds_strand_count_asymmetry.tsv all_species_composite_exon_shadow_strand_count_asymmetry.tsv all_species_multi_exon_shadow_strand_count_asymmetry.tsv all_species_scfr_strand_count_asymmetry.tsv all_species_single_exon_shadow_strand_count_asymmetry.tsv all_species_shadow_strand_count_asymmetry.tsv

cd $base_directory/SCFR-main/exon_shadow/
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
#Strand asymmetry of shadow count
awk '{c[$6]++} END {for (k in c) print c[k], k}' "$species"/"$species"_coding_exons.bed | grep -E '\+|\-' | sort -t' ' -k2,2r | paste -s -d " " | awk '{print$0,$1-$3,($1-$3)/($1+$3)}' | sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_strand_cds_count_asymmetry.tsv
awk '{c[$6]++} END {for (k in c) print c[k], k}' "$species"/"$species"_scfr_all.bed | grep -E '\+|\-' | sort -t' ' -k2,2r | paste -s -d " " | awk '{print$0,$1-$3,($1-$3)/($1+$3)}' |  sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_strand_scfr_count_asymmetry.tsv
awk '{c[$6]++} END {for (k in c) print c[k], k}' "$species"/"$species"_filtered_combined_non_zero.tsv | grep -E '\+|\-' | sort -t' ' -k2,2r | paste -s -d " " | awk '{print$0,$1-$3,($1-$3)/($1+$3)}' |  sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_strand_total_shadow_count_asymmetry.tsv
awk '$15>0' "$species"/"$species"_filtered_combined_non_zero.tsv | awk '{c[$6]++} END {for (k in c) print c[k], k}' | grep -E '\+|\-' | sort -t' ' -k2,2r | paste -s -d " " | awk '{print$0,$1-$3,($1-$3)/($1+$3)}' |  sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_strand_up_shadow_positive_count_asymmetry.tsv
awk '$16>0' "$species"/"$species"_filtered_combined_non_zero.tsv | awk '{c[$6]++} END {for (k in c) print c[k], k}' | grep -E '\+|\-' | sort -t' ' -k2,2r | paste -s -d " " | awk '{print$0,$1-$3,($1-$3)/($1+$3)}' |  sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_strand_down_shadow_positive_count_asymmetry.tsv
#Strand asymmetry of shadow length
awk '{sum[$6] += ($15 - $14 + 1)} END {for (i in sum) print i, sum[i]}' "$species"/"$species"_coding_exons.bed | grep -E '\+|\-' | sort -t' ' -k1,1r | paste -s -d " " | awk '{print$2,$1,$4,$3}' | awk '{print$0,$1-$3,($1-$3)/($1+$3)}' | sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_strand_cds_length_asymmetry.tsv
awk '{sum[$6] += ($3 - $2 + 1)} END {for (i in sum) print i, sum[i]}' "$species"/"$species"_scfr_all.bed | grep -E '\+|\-' | sort -t' ' -k1,1r | paste -s -d " " | awk '{print$2,$1,$4,$3}' | awk '{print$0,$1-$3,($1-$3)/($1+$3)}' | sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_strand_scfr_length_asymmetry.tsv
awk '{s15[$6]+=$15; s16[$6]+=$16} END {for (i in s15) print i, s15[i], s16[i]}'  "$species"/"$species"_filtered_combined_non_zero.tsv | grep -v "strand" | awk '{print$1,$2+$3}' | grep -E '\+|\-' | sort -t' ' -k2,2r | paste -s -d " " | awk '{print$2,$1,$4,$3}' | awk '{print$0,$1-$3,($1-$3)/($1+$3)}' | sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_strand_total_shadow_length_asymmetry.tsv
awk '{s15[$6]+=$15; s16[$6]+=$16} END {for (i in s15) print i, s15[i], s16[i]}'  "$species"/"$species"_filtered_combined_non_zero.tsv | grep -v "strand" | awk '{print$1,$2}' | grep -E '\+|\-' | sort -t' ' -k2,2r | paste -s -d " " | awk '{print$2,$1,$4,$3}' | awk '{print$0,$1-$3,($1-$3)/($1+$3)}' | sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_strand_up_shadow_length_asymmetry.tsv
awk '{s15[$6]+=$15; s16[$6]+=$16} END {for (i in s15) print i, s15[i], s16[i]}'  "$species"/"$species"_filtered_combined_non_zero.tsv | grep -v "strand" | awk '{print$1,$3}' | grep -E '\+|\-' | sort -t' ' -k2,2r | paste -s -d " " | awk '{print$2,$1,$4,$3}' | awk '{print$0,$1-$3,($1-$3)/($1+$3)}' | sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_strand_down_shadow_length_asymmetry.tsv
#Directional asymmetry of shadow count
awk '$15>0 {c15++} $16>0 {c16++} END {print"up", c15; print"down", c16}' "$species"/"$species"_filtered_combined_non_zero.tsv | paste -s -d " " | awk '{print$2,$1,$4,$3}' | awk '{print$0,$1-$3,($1-$3)/($1+$3)}' | sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_direction_total_shadow_count_asymmetry.tsv
#Directional asymmetry of shadow length
awk '{a+=$15; b+=$16} END{print a,b}' "$species"/"$species"_filtered_combined_non_zero.tsv | awk '{print$1,"up",$2,"down",$1-$2,($1-$2)/($1+$2)}' | sed "s/^/$species /g" | sed 's/[ ]\+/\t/g' >> asymmetry/all_species_direction_total_shadow_length_asymmetry.tsv
done

#Exitrons
cd $base_directory/SCFR-main/exon_shadow/
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
awk '{c[$14]++} END {for (k in c) print c[k], k}' "$species"/"$species"_exitron_candidates_filtered.tsv | grep -E '\+|\-' | sort -t' ' -k2,2r | paste -s -d " " | awk '{print$2,$1,$4,$3}' | awk '{print$0,$2-$4, ($2-$4)/($2+$4)}' >> asymmetry/all_species_strand_exitron_count_asymmetry
awk '{s18[$14]+=$18} END {for (i in s18) print i, s18[i]}' "$species"/"$species"_exitron_candidates_filtered_exon_features.tsv | grep -E '\+|\-' | sort -t' ' -k1,1r | paste -s -d " " | awk '{print$0,$2-$4, ($2-$4)/($2+$4)}' >> asymmetry/all_species_strand_exitron_length_asymmetry
done

#Plot Strand asymmetry

#################################################################################################################################################################################################################################################################################################
#Shadow length enrichment

cd $base_directory/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
echo ">"$species
rm -r exon_shadow/"$species"/enrichment_shadow_length
mkdir exon_shadow/"$species"/enrichment_shadow_length
cd exon_shadow/"$species"/enrichment_shadow_length
#plot
Rscript $base_directory/SCFR-main/my_scripts/exon_shadow/enrichment_shadow_length.R ../"$species"_filtered_combined_non_zero.tsv enrichment_"$species"
cd $base_directory/SCFR-main/
done

cd $base_directory/SCFR-main
mkdir -p $base_directory/SCFR-main/exon_shadow/plot/enrichment_shadow_length
find exon_shadow/ -mindepth 3 -maxdepth 3 -name "*top_enriched_shadow_lengths.pdf" | xargs -n1 sh -c 'cp $0 $base_directory/SCFR-main/exon_shadow/plot/enrichment_shadow_length/'
find exon_shadow/ -mindepth 3 -maxdepth 3 -name "*top_enriched_shadow_length_bins.pdf" | xargs -n1 sh -c 'cp $0 $base_directory/SCFR-main/exon_shadow/plot/enrichment_shadow_length/'
find exon_shadow/ -mindepth 3 -maxdepth 3 -name "enrichment_*_shadow_length_1_500.pdf" | xargs -n1 sh -c 'cp $0 $base_directory/SCFR-main/exon_shadow/plot/enrichment_shadow_length/'

#################################################################################################################################################################################################################################################################################################
#Composition of exon-shadow

#rename chromosomes fro extracting shadow fasta
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do 
gr=$(ls $base_directory/SCFR-main/genome_reports/GC*.tsv | grep "$species")
genome=$(readlink -f $base_directory/SCFR-main/genomes/"$species"/GC*.fna)
cd $base_directory/SCFR-main/genomes/$species
awk -F "\t" -v OFS="\t" '{ for(N=1; N<=NF; N++) if($N=="") $N="-" } 1' $gr | sed 's/[ ]\+/_/g' | awk 'NR>1{print$6,$9}' | tr " " "\t" > map.tsv
gawk -i inplace 'BEGIN{while ((getline < "map.tsv") > 0)
    map[$1] = $2}
/^>/{for(k in map) if(index($0, k)) sub(/^>/, ">" map[k] " ")} {print}' $genome
unset gr genome
done

#Extract fasta of exon shadow:
cd $base_directory/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do 
genome=$(readlink -f $base_directory/SCFR-main/genomes/"$species"/GC*.fna)
echo ">"$species ":" $genome
cd $base_directory/SCFR-main/exon_shadow/"$species"
mkdir composition
cd composition
#Get bed coordinates for upstream & downstream shadow (head = upstream, gene, transcript ID, first exon number, last exon number, first exon start, last exon start, frame)
awk '$15>0' ../"$species"_filtered_combined_non_zero.tsv | awk 'NR>1{if($6=="+") print$1,$2,$8,"upstream_"$10"_"$11"_"$17"_"$18"_"$8"_"$9"_"$4,1,$6; else print$1,$9,$3,"upstream_"$10"_"$11"_"$17"_"$18"_"$8"_"$9"_"$4,1,$6}' OFS="\t" > upstream_"$species"_filtered_combined_non_zero.bed
awk '$16>0' ../"$species"_filtered_combined_non_zero.tsv | awk 'NR>1{if($6=="+") print$1,$9,$3,"downstream_"$10"_"$11"_"$17"_"$18"_"$8"_"$9"_"$4,1,$6; else print$1,$2,$8,"downstream_"$10"_"$11"_"$17"_"$18"_"$8"_"$9"_"$4,1,$6}' OFS="\t" > downstream_"$species"_filtered_combined_non_zero.bed
#Extract fasta sequence
bedtools getfasta -fi $genome -bed upstream_"$species"_filtered_combined_non_zero.bed -s -name+ > upstream_"$species"_filtered_combined_non_zero.fa
bedtools getfasta -fi $genome -bed downstream_"$species"_filtered_combined_non_zero.bed -s -name+ > downstream_"$species"_filtered_combined_non_zero.fa
#Get Percentage GC content
#infoseq -sequence upstream_"$species"_filtered_combined_non_zero.fa -auto -only -name -length -pgc > upstream_"$species"_filtered_combined_non_zero_pgc.out
#infoseq -sequence downstream_"$species"_filtered_combined_non_zero.fa -auto -only -name -length -pgc > downstream_"$species"_filtered_combined_non_zero_pgc.out
unset genome
cd $base_directory/SCFR-main
done

#Get intron sequences (12m53.760s)
cd $base_directory/SCFR-main/
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
cd exon_shadow/"$species"
genome=$(readlink -f $base_directory/SCFR-main/genomes/"$species"/GC*.fna)
#Extract intron fasta
time $base_directory/SCFR-main/my_scripts/exon_shadow/extract_introns.sh "$species"_coding_exons.bed $genome "$species"_introns.fa
#filter duplicate introns
time seqkit rmdup -s < "$species"_introns.fa > "$species"_introns_filtered.fa
unset genome
cd $base_directory/SCFR-main/
done 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Add exon related data to exitron results (40m48.259s)
cd $base_directory/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
cd exon_shadow/"$species"
echo "chr start end frame filler frame chr exon_1_start exon_1_end exon_2_start exon_2_end gene transcript exon_strand exon_frame intron_start intron_end intron_length merged_txs first_exon_number last_exon_number gene_tx_count first_exon_tx_count last_exon_tx_count first_exon_order last_exon_order first_exon_sharing last_exon_sharing first_exon_splicing last_exon_splicing" > "$species"_exitron_candidates_filtered_exon_features.tsv
time while read ee
do
e1=$(echo $ee | awk '{print$7,$8,$9,$12,$13,$14,$15}' OFS="\t")
e2=$(echo $ee | awk '{print$7,$10,$11,$12,$13,$14,$15}' OFS="\t")
em1=$(grep "$e1" "$species"_coding_exons.bed | awk '{print$8,$9,$10,$11,$12,$13}')
em2=$(grep "$e2" "$species"_coding_exons.bed | awk '{print$8,$10,$11,$12,$13}')
em=$(echo $em1 $em2 | awk '{print$1,$7,$2,$3,$8,$4,$9,$5,$10,$6,$11,$12}')
echo $ee $em
unset e1 e2 em1 em2 em
done < <(sed '1d' "$species"_exitron_candidates_filtered.tsv) >> "$species"_exitron_candidates_filtered_exon_features.tsv
sed 's/[ \t]\+/\t/g' "$species"_exitron_candidates_filtered_exon_features.tsv -i
unset ee
cd $base_directory/SCFR-main
done

#Get exitron fasta
cd $base_directory/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
cd exon_shadow/"$species"
awk 'NR>1{print$7,$16,$17,$12"_"$13"_"$20"_"$21,1,$14}' OFS="\t" "$species"_exitron_candidates_filtered_exon_features.tsv > composition/"$species"_exitron_candidates_filtered_exon_features.bed
bedtools getfasta -fi $base_directory/SCFR-main/genomes/"$species"/GC*.fna -bed composition/"$species"_exitron_candidates_filtered_exon_features.bed -name+ -s > composition/"$species"_exitron_candidates_filtered_exon_features.fa
cd $base_directory/SCFR-main
done

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calculate composition:

#26m6.043s
cd $base_directory/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
cd exon_shadow/"$species"/composition
up=$(readlink -f $base_directory/SCFR-main/exon_shadow/"$species"/composition/upstream_"$species"_filtered_combined_non_zero.fa)
dn=$(readlink -f $base_directory/SCFR-main/exon_shadow/"$species"/composition/downstream_"$species"_filtered_combined_non_zero.fa)
in=$(readlink -f $base_directory/SCFR-main/exon_shadow/"$species"/"$species"_introns_filtered.fa)
ex=$(readlink -f $base_directory/SCFR-main/Fourier_analysis/genes/"$species"/GC*_cds.fa)
seqkit rmdup -s < $ex | sed -E 's/^>lcl\|([^[:space:]]+).*/>\1/' > "$species"_exon.fa
xt=$(readlink -f "$species"_exitron_candidates_filtered_exon_features.fa)
#calculate GC, AT content, GC skew, average G/C stretch, average A/T stretch
$base_directory/SCFR-main/my_scripts/exon_shadow/get_composition.sh $up > "$species"_upstream_composition.out
$base_directory/SCFR-main/my_scripts/exon_shadow/get_composition.sh $dn > "$species"_downstream_composition.out
$base_directory/SCFR-main/my_scripts/exon_shadow/get_composition.sh $in > "$species"_intron_composition.out
$base_directory/SCFR-main/my_scripts/exon_shadow/get_composition.sh "$species"_exon.fa > "$species"_exon_composition.out
$base_directory/SCFR-main/my_scripts/exon_shadow/get_composition.sh $xt > "$species"_exitron_composition.out
unset up dn in ex xt
cd $base_directory/SCFR-main
done

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calculate & plot distribution

set -euo pipefail

BASE=$base_directory/SCFR-main
OUTDIR=$BASE/exon_shadow/plot/composition
STATS_PY=$BASE/my_scripts/get_stats.py
species_list=(human bonobo chimpanzee gorilla borangutan sorangutan gibbon)

mkdir -p "$OUTDIR"
header="species\tseq_type\tN\tmin\tmax\tmean\tq1\tmedian\tq3"
echo -e "$header" > "$OUTDIR/gc_distribution.tsv"
echo -e "$header" > "$OUTDIR/at_distribution.tsv"
echo -e "$header" > "$OUTDIR/gc_skew_distribution.tsv"
echo -e "$header" > "$OUTDIR/gc3_distribution.tsv"
echo -e "$header" > "$OUTDIR/avg_gc_stretch_distribution.tsv"
echo -e "$header" > "$OUTDIR/avg_at_stretch_distribution.tsv"
cd "$BASE"

time for species in "${species_list[@]}"; do
    echo ">$species"
    COMP_DIR=$BASE/exon_shadow/$species/composition
    cd "$COMP_DIR"
    # Cache file list once per species
    mapfile -t files < <(find . -name "*_composition.out")
    # column → output file mapping
    declare -A cols=(
        [2]=gc_distribution.tsv
        [3]=at_distribution.tsv
        [4]=gc_skew_distribution.tsv
        [5]=gc3_distribution.tsv
        [6]=avg_gc_stretch_distribution.tsv
        [7]=avg_at_stretch_distribution.tsv
    )
    for col in "${!cols[@]}"; do
        for seq in "${files[@]}"; do
            id=$(basename "$seq" | cut -f2 -d "_")
            awk -v c="$col" 'NR>1{print $c}' "$seq" \
            | python3 "$STATS_PY" \
            | egrep "Count|^Minimum|^Maximum|^Mean|^Median|^Q1|^Q3" \
            | awk -F ":" '{print $NF}' \
            | tr -d " ," \
            | paste -s -d " " \
            | sed "s/^/$species $id /"
        done \
        | sed 's/[ ]\+/\t/g' \
        >> "$OUTDIR/${cols[$col]}"
    done
    cd "$BASE"
done

#Plot distribution
for i in $(ls *distribution.tsv)
do
p=$(echo $i | sed 's/_distribution.tsv//g' | sed 's/gc/GC/g' | sed 's/at/AT/g' | sed 's/AT\b/& content/g' | sed 's/GC\b/& content/g' | sed 's/avg/Average /g' | tr "_" " ")
o=$(echo $i | sed 's/\.tsv$/\.pdf/g')
echo $p $o
Rscript $base_directory/SCFR-main/my_scripts/exon_shadow/plot_composition_distribution.R $i $o "$p distribution" "$p"
unset p o
done

#################################################################################################################################################################################################################################################################################################
#Gene set Enrichment

#Get genes with exitrons
cd $base_directory/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
cd exon_shadow/"$species"
awk 'NR>1{print$12}' "$species"_exitron_candidates_filtered_exon_features.tsv | sort | uniq > "$species"_genes_with_exitron_candidates.txt
cd $base_directory/SCFR-main
done

#################################################################################################################################################################################################################################################################################################
#Fourier analysis

cd $base_directory/SCFR-main/
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
cd exon_shadow/"$species"/composition
echo "- Running: Upstream"
seqtk seq -L 30 upstream_"$species"_filtered_combined_non_zero.fa > upstream_"$species"_filtered_combined_non_zero_filtered.fa
python3 $base_directory/SCFR-main/my_scripts/orf_parallel_fft_motif_report_grouped.py upstream_"$species"_filtered_combined_non_zero_filtered.fa -t 32 -o upstream_fourier
python3 $base_directory/SCFR-main/my_scripts/scfr_fourier_chromosome_wise_summary.py upstream_fourier --top 3 --cores 32
echo "- Running: Downstream"
seqtk seq -L 30 downstream_"$species"_filtered_combined_non_zero.fa > downstream_"$species"_filtered_combined_non_zero_filtered.fa
python3 $base_directory/SCFR-main/my_scripts/orf_parallel_fft_motif_report_grouped.py downstream_"$species"_filtered_combined_non_zero_filtered.fa -t 32 -o downstream_fourier
python3 $base_directory/SCFR-main/my_scripts/scfr_fourier_chromosome_wise_summary.py downstream_fourier --top 3 --cores 32
echo "- Running: Extron"
seqtk seq -L 30 "$species"_exitron_candidates_filtered_exon_features.fa > "$species"_exitron_candidates_filtered_exon_features_filtered.fa
python3 $base_directory/SCFR-main/my_scripts/orf_parallel_fft_motif_report_grouped.py "$species"_exitron_candidates_filtered_exon_features_filtered.fa -t 32 -o exitron_fourier
python3 $base_directory/SCFR-main/my_scripts/scfr_fourier_chromosome_wise_summary.py exitron_fourier --top 3 --cores 32
cd $base_directory/SCFR-main/
done

#Get fourier peaks (11m50.294s)
cd $base_directory/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
echo ">"$species
awk '$3>0' exon_shadow/"$species"/composition/upstream_fourier/chromosome_wise_summary/summary.tsv > exon_shadow/"$species"/composition/upstream_fourier/"$species"_upstream_with_peaks.tsv
awk '$3 > 0 && $4 ~ /^0\.3/' exon_shadow/"$species"/composition/upstream_fourier/chromosome_wise_summary/summary.tsv > exon_shadow/"$species"/composition/upstream_fourier/"$species"_3_periodicity_upstream.tsv
awk '$3>0' exon_shadow/"$species"/composition/downstream_fourier/chromosome_wise_summary/summary.tsv > exon_shadow/"$species"/composition/downstream_fourier/"$species"_downstream_with_peaks.tsv
awk '$3 > 0 && $4 ~ /^0\.3/' exon_shadow/"$species"/composition/downstream_fourier/chromosome_wise_summary/summary.tsv > exon_shadow/"$species"/composition/downstream_fourier/"$species"_3_periodicity_downstream.tsv
awk '$3>0' exon_shadow/"$species"/composition/exitron_fourier/chromosome_wise_summary/summary.tsv > exon_shadow/"$species"/composition/exitron_fourier/"$species"_exitron_with_peaks.tsv
awk '$3 > 0 && $4 ~ /^0\.3/' exon_shadow/"$species"/composition/exitron_fourier/chromosome_wise_summary/summary.tsv > exon_shadow/"$species"/composition/exitron_fourier/"$species"_3_periodicity_exitron.tsv
cd $base_directory/SCFR-main/
done

#Plot fourier frequence density
cd $base_directory/SCFR-main
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
cd $base_directory/SCFR-main/exon_shadow/"$species"/composition/
Rscript $base_directory/SCFR-main/my_scripts/plot_fourier_frequencies.R upstream_fourier/"$species"_upstream_with_peaks.tsv "$species"_upstream_fourier_fourier_freq.pdf
Rscript $base_directory/SCFR-main/my_scripts/plot_fourier_frequencies.R downstream_fourier/"$species"_downstream_with_peaks.tsv "$species"_downstream_fourier_fourier_freq.pdf
Rscript $base_directory/SCFR-main/my_scripts/plot_fourier_frequencies.R exitron_fourier/"$species"_exitron_with_peaks.tsv "$species"_exitron_fourier_fourier_freq.pdf
cd $base_directory/SCFR-main
done

#################################################################################################################################################################################################################################################################################################
#Summary 

#Print summary of exitron results
cd $base_directory/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
cd exon_shadow/"$species"/composition
i1=$(grep ">" "$species"_exitron_candidates_filtered_exon_features_filtered.fa -c)
i2=$(awk 'END{print NR-1}' exitron_fourier/chromosome_wise_summary/summary.tsv)
i3=$(awk 'END{print NR-1}' exitron_fourier/exitrons_with_peaks.tsv)
i4=$(wc -l < exitron_fourier/3_periodicity_exitrons.tsv)
echo $species $i1 $i2 $i3 $i4
unset i1 i2 i3 i4
cd $base_directory/SCFR-main/
done | sed '1i Species input_count output_count output_with_peaks output_with_3_periodicity' | sed 's/[ ]\+/\t/g' > $base_directory/SCFR-main/exon_shadow/plot/exitron_summary.tsv

#Print summary of upstream results
cd $base_directory/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
cd exon_shadow/"$species"/composition
i1=$(grep ">" upstream_"$species"_filtered_combined_non_zero_filtered.fa -c)
i2=$(awk 'END{print NR-1}' upstream_fourier/chromosome_wise_summary/summary.tsv)
i3=$(awk 'END{print NR-1}' upstream_fourier/upstream_with_peaks.tsv)
i4=$(wc -l < upstream_fourier/3_periodicity_upstream.tsv)
echo $species $i1 $i2 $i3 $i4
unset i1 i2 i3 i4
cd $base_directory/SCFR-main/
done | sed '1i Species input_count output_count output_with_peaks output_with_3_periodicity' | sed 's/[ ]\+/\t/g' > $base_directory/SCFR-main/exon_shadow/plot/upstream_summary.tsv

#Print summary of downstream results
cd $base_directory/SCFR-main/
time for species in human bonobo borangutan sorangutan chimpanzee gorilla gibbon 
do
cd exon_shadow/"$species"/composition
i1=$(grep ">" downstream_"$species"_filtered_combined_non_zero_filtered.fa -c)
i2=$(awk 'END{print NR-1}' downstream_fourier/chromosome_wise_summary/summary.tsv)
i3=$(awk 'END{print NR-1}' downstream_fourier/downstream_with_peaks.tsv)
i4=$(wc -l < downstream_fourier/3_periodicity_downstream.tsv)
echo $species $i1 $i2 $i3 $i4
unset i1 i2 i3 i4
cd $base_directory/SCFR-main/
done | sed '1i Species input_count output_count output_with_peaks output_with_3_periodicity' | sed 's/[ ]\+/\t/g' > $base_directory/SCFR-main/exon_shadow/plot/downstream_summary.tsv

#################################################################################################################################################################################################################################################################################################
