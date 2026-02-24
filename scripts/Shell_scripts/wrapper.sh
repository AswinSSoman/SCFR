####################################################################################################################################################################################################################################################################################################################
#                                                                                                                              MAIN SCRIPT FOR SCFR PROJECT
####################################################################################################################################################################################################################################################################################################################

####################################################################################################################################################################################################################################################################################################################
#1. Prepare and collect data & script

#set base directory
#base_directory=""

#Download github repository
	cd $base_directory
	wget https://github.com/ceglab/SCFR/archive/refs/heads/main.zip
	unzip main.zip

#Create folders
	cd $base_directory/SCFR-main
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo $species
	mkdir -p chrs/"$species"
	mkdir -p SCFR/"$species"
	mkdir -p SCFR/gaps/"$species"
	mkdir -p genes/"$species"
	mkdir -p GC/"$species"
	mkdir -p SCFR_all
	mkdir -p scripts
	mkdir -p QC/genome_metadata
	mkdir -p SCFR_lists
	mkdir -p Fourier_analysis
	mkdir -p gene_deserts
	mkdir -p SCFR_summaries
	done

#mkdir github
	#mv SCFR_lists PCA mv SCFR_summaries github/

#List of genomic accesions for analysis
	echo -e "GCA_009914755.4\thuman\nGCA_028858775.3\tchimpanzee\nGCA_029289425.3\tbonobo\nGCA_028878055.3\tgibbon\nGCA_029281585.3\tgorilla\nGCA_028885625.3\tborangutan\nGCA_028885655.3\tsorangutan" > QC/genome_accessions

####################################################################################################################################################################################################################################################################################################################
#2. Download all the seven primate genomes

	cd $base_directory/SCFR-main
	time while read species
	do
	genome=$(echo $species | awk '{print$1}')
	name=$(echo $species | awk '{print$2}')
	echo ">"$genome $name
	#Download genome & annotation files
	#time datasets download genome accession $genome --include genome,gtf,seq-report --filename $genome".zip"
	unzip -o $genome".zip"
	cd ncbi_dataset/data/$genome
	#If genome size is too big or the downloads gets failed, use dehydrate
	#datasets download genome accession $genome --include genome,gtf,seq-report --dehydrated --filename $genome.zip
	#unzip $genome.zip -od $genome
	#datasets rehydrate --directory $genome
	#Get accessions
	cat sequence_report.jsonl | json2xml | xtract -pattern opt -def "-"  -element genbankAccession refseqAccession | awk '$2!="-"' > $genome"_genbank_refseq_accession"
	#Add refseq accesion in genome
	time awk 'NR==FNR {map[$1]=$2; next}  /^>/ {hdr=$0; for(k in map) {if(hdr ~ ">"k"($|[[:space:]])") {sub(">"k, ">"map[k]" "k); break}}} {print}' $genome"_genbank_refseq_accession" $genome*".fna" > $genome".fna"
	#Get chromosome fasta files
	outdir="$base_directory/SCFR-main/chrs/$name"
	awk -v outdir="$outdir" '/^>/ { if (out) close(out) ; out = outdir "/" substr($1,2) ".fasta" } { print >> out }' "${genome}.fna"
	unset genome name
	cd $base_directory/SCFR-main
	done < QC/genome_accessions

#Delete mitochondrial genomes
	cd $base_directory/SCFR-main
	grep -vf <(awk '{print$1}' QC/genome_metadata/all_species_chr_length | grep -v "refseq_accession") <(find chrs -name "*.fasta") | xargs rm

####################################################################################################################################################################################################################################################################################################################
#3. Quality check

cd $base_directory/SCFR-main
mkdir -p QC/genome_metadata

#3.1. get genome metadata
	cd $base_directory/SCFR-main/QC
	awk '{print$1}' genome_accessions \
	| xargs -n1 sh -c 'datasets summary genome accession $0 --as-json-lines | json2xml | xtract -pattern root -def "-" -element organism_name assembly_name assembly_level accession contig_l50 contig_n50 scaffold_l50 scaffold_n50 total_number_of_chromosomes total_sequence_length total_ungapped_length genome_coverage' \
	| sed 's/ /_/g' | sed '1i Species Assembly_name Assembly_level Genbank_accession Refseq_accession Contig_l50 Contig_n50 Scaffold_l50 Scaffold_n50 Total_number_of_chromosomes Total_sequence_length Total_ungapped_length Genome_coverage' \
	| sed 's/PRJNA[^ \t]\+//g' | sed 's/SAMN[^ \t]\+//g' | awk '!($1=="Homo_sapiens" || $1=="Species") { $1=""; sub(/^ +/, ""); }1' > genome_metadata/genome_metadata
#Get chromosome metadata using datasets
	cd $base_directory/SCFR-main/QC
	while read i
	  do
	  echo $i
	  i1=$(echo $i | awk '{print$1}')
	  i2=$(echo $i | awk '{print$2}')
	  datasets summary genome accession "$i1" --report sequence --as-json-lines > genome_metadata/$i2"_"$i1"_genome_metadata.json"
	  cat genome_metadata/$i2"_"$i1"_genome_metadata.json" | json2xml | xtract -pattern opt -def "-" -element assembly_accession assembly_unit assigned_molecule_location_type chr_name gc_count gc_percent genbank_accession length refseq_accession role sequence_name | tr " " "_" | sed '1i assembly_accession assembly_unit assigned_molecule_location_type chr_name gc_count gc_percent genbank_accession length refseq_accession role sequence_name' | column -t > genome_metadata/$i2"_"$i1"_genome_metadata.txt"
	  unset i1 i2
	done < genome_accessions
#Get length of chromosomes 
	cd $base_directory/SCFR-main/QC/genome_metadata
	find . -name "*_genome_metadata.txt" | xargs -n1 sh -c 'awk "{print\$9,\$8}" $0' | grep -v "refseq_accession" | sed '1i refseq_accession Length' | grep -v "^\-" > all_species_chr_length

#3.2. Calculate length of chromosome from fetched genome sequences
	cd $base_directory/SCFR-main/chrs
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	 do
	 cd $species
	 for chr in $(ls *.fasta)
	 do
	 length=$(awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' $chr | awk '{print$1,$NF}')
	 echo $species $length
	 unset length
	 done
	 unset chr
	 cd $base_directory/SCFR-main/chrs
	done | sed '1i Species chromosome length' | column -t > $base_directory/SCFR-main/QC/genome_metadata/fetched_length

#3.3. Report errors in input: Compare lengths between fetched and database chromosomes
	cd $base_directory/SCFR-main/QC/genome_metadata
	for chr in $(awk '{print$1}' all_species_chr_length | egrep -v "refseq_accession|-")
	do
	m1=$(awk -v a="$chr" '$1==a' all_species_chr_length)
	m2=$(awk -v a="$chr" '$2==a' fetched_length)
	echo $m1 $m2
	unset m1 m2
	done | awk '{if($2==$5) print$0,"same"; else print$0,"different"}' |  sed '1i ncbi_accession ncbi_length species fetched_accession fetched_length Compare_length' | column -t > compare_chromosome_length
	
####################################################################################################################################################################################################################################################################################################################
#4.Identify and summarize SCFRs

#4.1. Identify SCFRs (51.6833 mins)
	cd $base_directory/SCFR-main
	start_time=$(date +%s)
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo $species
	for chr in `ls -1 chrs/$species/*.fasta|cut -f 3 -d '/'|sed 's/\.fasta//g'`
	do
	(
	echo $chr
	python3 scripts/find_stop_codon_free_regions_with_reverse_gap_report.py chrs/$species/"$chr".fasta --gaps SCFR/gaps/"$species"/gap_regions_"$chr".bed > SCFR/"$species"/"$chr".fasta.SCFRs.out
	) &
	done
	wait
	done
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#4.2. Collate all gaps
cd $base_directory/SCFR-main
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo $species
	cat SCFR/gaps/"$species"/gap*.bed|sort -k1,1 -k2n,2 > SCFR_all/gaps_"$species".bed
	done

#4.3.1. Get summary of SCFRs (35.4167 mins)
	cd $base_directory/SCFR-main
	start_time=$(date +%s)
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	(
	echo $species
	cat SCFR/"$species"/*.fasta.SCFRs.out | sort -k1,1 -k2n,2 | bedtools intersect -v -a stdin -b SCFR_all/gaps_"$species".bed > SCFR_all/"$species"_SCFR_all.out
	) &
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e
#4.3.2. Since running 6 species need huge memory run 3 at a time (1.67528 hours ~100.517 mins)
	start_time=$(date +%s)
	max_jobs=3
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	(
	echo $species
	Rscript scripts/summarize_SCFR_bed_frames_all.R SCFR_all/"$species"_SCFR_all.out
	) &
	# If max_jobs reached, wait for ONE job to finish (not all)
	if (( $(jobs -r | wc -l) >= max_jobs )); then
	wait -n
	fi
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e
	unset max_jobs
#4.3.3. Combine the summaries
	ls | egrep "_SCFR_all_frame_summary_by_chromosome|_SCFR_all_frame_summary_all" | xargs -n1 sh -c 'mv $0 SCFR_summaries/'
	cd $base_directory/SCFR-main/SCFR_summaries
	ls | grep "_SCFR_all_frame_summary_all.csv" | xargs -n1 sh -c 'grep -v Chromosome $0 | sed "s/^/$0 /g"' | sed 's/_SCFR_all_frame_summary_all.csv//g' | tr "," " " | tr -d '"' | sed '1i Species Chromosome Frame N Min Q1 Median Mean Q3 Max SD P95 P99 Q_1Kb Q_5Kb Q_10Kb' > all_species_chromsome_frame_wise_summary
	ls | grep "_SCFR_all_frame_summary_all.csv" | xargs -n1 sh -c 'grep ALL $0 | sed "s/^/$0 /g"' | sed 's/_SCFR_all_frame_summary_all.csv//g' | tr "," " " | tr -d '"' | sed '1i Species Chromosome Frame N Min Q1 Median Mean Q3 Max SD P95 P99 Q_1Kb Q_5Kb Q_10Kb' > all_species_frame_wise_summary

#4.4. Plot data
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	plot_all_species_SCFR_length_stats.R 
	window_wise_all_species_scfr_coding_stats.csv
#Make SCFR length stats plot
#NOTE: Rename species common names before plotting
	cd $base_directory/SCFR-main/SCFR_all
	Rscript $base_directory/SCFR-main/my_scripts/plot_all_species_SCFR_length_stats.R all_scfr_length_stats.tsv all_scfr_length_stats.pdf

####################################################################################################################################################################################################################################################################################################################
#5. Strand asymmetry

#5.1. Get strand asymmetry of SCFRs per chromosome (11.3 mins)
	cd $base_directory/SCFR-main
	start_time=$(date +%s)
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	(
	echo $species
	python3 scripts/quantify_scfr_asymmetries_by_chrom.py SCFR_all/"$species"_SCFR_all.out SCFR_all/"$species"_SCFR_asymmetries_out.csv
	) &
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#5.2. Get strand assymetry of SCFRs in sliding windows across the genomes (9.45 mins)
	cd $base_directory/SCFR-main
	start_time=$(date +%s)
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	(
	echo $species
	python3 scripts/quantify_scfr_asymmetries_by_chrom_window.py SCFR_all/"$species"_SCFR_all.out --window-size 100000 --slide-size 50000 --output SCFR_all/"$species"_SCFR_asymmetries_out_win100000_slide50000.csv
	) &
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#5.3. Get window based strand asymmetry (6 secs)
	cd $base_directory/SCFR-main
	start_time=$(date +%s)
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	(
	echo $species
	Rscript scripts/version_3_plot_strand_asymmetry_sliding_extremes.R --input SCFR_all/"$species"_SCFR_asymmetries_out_win100000_slide50000.csv --pdf SCFR_all/"$species"_sliding_outliers.pdf --bed SCFR_all/"$species"_sliding_outlier_regions.bed --min_region_size 100000 --window_len 5 --min_hits 3
	) &
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#5.4. combine the summaries
	cd $base_directory/SCFR-main/SCFR_all
	ls | grep "_SCFR_asymmetries_out.csv" | xargs -n1 sh -c 'grep -v chrom $0 | sed "s/^/$0 /g"' | sed 's/_SCFR_asymmetries_out.csv//g' | tr "," " " | tr -d '"' \
		| sed '1i Species chrom f_count r_count strand_count_asym f_length r_length strand_length_asym frame1_count frame1_length frame2_count frame2_length frame3_count frame3_length frame-1_count frame-1_length frame-2_count frame-2_length frame-3_count frame-3_length' > all_species_chromosome_wise_SCFR_asymmetries
	ls | grep "_SCFR_asymmetries_out.csv" | xargs -n1 sh -c 'grep ALL $0 | sed "s/^/$0 /g"' | sed 's/_SCFR_asymmetries_out.csv//g' | tr "," " " | tr -d '"' \
		| sed '1i Species chrom f_count r_count strand_count_asym f_length r_length strand_length_asym frame1_count frame1_length frame2_count frame2_length frame3_count frame3_length frame-1_count frame-1_length frame-2_count frame-2_length frame-3_count frame-3_length' > all_species_SCFR_asymmetries

##################################################################################################################################################################################################################################################
#6. GC content & coding gene enrichment

#6.1. get GC content of the SCFRs (48.1667)
	cd $base_directory/SCFR-main
	start_time=$(date +%s)
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	(
	echo $species
	cat chrs/$species/*.fasta > tmp."$species".fasta
	bedtools nuc -fi tmp."$species".fasta -bed SCFR_all/"$species"_SCFR_all.out > SCFR_all/"$species"_SCFR_GC_all.out
	rm tmp."$species".fasta
	cat SCFR_all/"$species"_SCFR_GC_all.out|awk '$13>10000{print $0}'|grep -v "^#"|sed 's/:/\t/g'|cut -f 1-13|sort -k1,1 -k2n,2 > SCFR_all/"$species"_long_SCFRs.bed
	) &
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.2. Download all the seven primate genome annotation files (6.18333 mins)
	cd $base_directory/SCFR-main/genes
	start_time=$(date +%s)
	declare -A urls=(
	  [human]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz"
	  [chimpanzee]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.gtf.gz"
	  [bonobo]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.gtf.gz"
	  [gibbon]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gtf.gz"
	  [gorilla]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.gtf.gz"
	  [borangutan]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/625/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.gtf.gz"
	  [sorangutan]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.gtf.gz"
	)
	for sp in "${!urls[@]}"; do
	    mkdir -p "$sp"
	    (
	      cd "$sp"
	      wget -q "${urls[$sp]}"
	    ) &
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.3. Extract coding exons annotated in the GTF
	cd $base_directory/SCFR-main
	start_time=$(date +%s)
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	(
	echo $species
	gunzip genes/"$species"/*.gtf.gz
	gtf=`ls -1  genes/"$species"/*.gtf|cut -f 3 -d '/'|sed 's/\.gtf//g'`
	#Extract coding exon coordinates from GTF
	python3 scripts/gtf_to_bed.py genes/"$species"/"$gtf".gtf genes/"$species"/"$gtf".bed
	) &
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.4. Window-wise GC content of SCFR

#6.4.1. Bin the SCFR counts by length and GC & AT content
#For AT (19.7833 mins)
	cd $base_directory/SCFR-main/SCFR_all
	start_time=$(date +%s)
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	(
	echo $species
	#For AT
	cat "$species"_SCFR_GC_all.out | awk '{printf "%s %.2f\n", $13, $5}' | sed 's/ /\t/g' |awk ' { if ($1 < 1000) bin1 = int($1 / 100) * 100; else bin1 = int($1 / 1000) * 1000; bin2 = sprintf("%.1f", int($2 * 10) / 10); count[bin1, bin2]++; range[bin1] = ($1 < 1000) ? 100 : 1000; } END { for (key in count) { split(key, bins, SUBSEP); r = range[bins[1]]; printf "%d-%d\t%s-%s\t%d\n", bins[1], bins[1]+r-1, bins[2], sprintf("%.1f", bins[2]+0.1), count[key]; } } ' > "$species"_AT_bins.out
	) &
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e
#For GC (17.9833 mins)
	cd $base_directory/SCFR-main/SCFR_all
	start_time=$(date +%s)
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	(
	echo $species
	#For GC
	cat "$species"_SCFR_GC_all.out | awk '{printf "%s %.2f\n", $13, $6}' | sed 's/ /\t/g' |awk ' { if ($1 < 1000) bin1 = int($1 / 100) * 100; else bin1 = int($1 / 1000) * 1000; bin2 = sprintf("%.1f", int($2 * 10) / 10); count[bin1, bin2]++; range[bin1] = ($1 < 1000) ? 100 : 1000; } END { for (key in count) { split(key, bins, SUBSEP); r = range[bins[1]]; printf "%d-%d\t%s-%s\t%d\n", bins[1], bins[1]+r-1, bins[2], sprintf("%.1f", bins[2]+0.1), count[key]; } } ' > "$species"_GC_bins.out
	) &
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.5. Find genes intersecting with long SCFRs

	cd $base_directory/SCFR-main
	start_time=$(date +%s)
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	(
	echo $species
	gtf=`ls -1  genes/"$species"/*.gtf|cut -f 3 -d '/'|sed 's/\.gtf//g'`
	bedtools intersect -a SCFR_all/"$species"_long_SCFRs.bed -b genes/"$species"/"$gtf".bed -wa -wb|cut -f 5,13,17|sort -u > SCFR_all/"$species"_genes_of_interest.txt
	) &
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#6.6. Plot a 2 dimensional histogram of length vs AT content and label the SCFR longer than 10 Kb that overlap coding exons
	Rscript scripts/combined_2dhist.r

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#6.7. Summaries

#Get the sum of all SCFR counts across all length bins for each GC bins (12.5333 mins)
	cd $base_directory/SCFR-main/SCFR_all
	for species in bonobo borangutan chimpanzee gibbon gorilla sorangutan human; do awk -v a="$species" '{counts[$2]+=$3} END{for(gc in counts)print a,gc, counts[gc]}' ${species}_GC_bins.out | sort -k1,1 -k2,2n; done | column -t > GC_all_species_total_SCFR_count_across_all_length_bins
	for species in bonobo borangutan chimpanzee gibbon gorilla sorangutan human; do awk -v a="$species" '{counts[$2]+=$3} END{for(gc in counts)print a,gc, counts[gc]}' ${species}_AT_bins.out | sort -k1,1 -k2,2n; done | column -t > AT_all_species_total_SCFR_count_across_all_length_bins
	for species in bonobo borangutan chimpanzee gibbon gorilla sorangutan human; do awk -v a="$species" '{counts[$1]+=$3} END{for(len in counts)print a,len, counts[len]}' ${species}_GC_bins.out | sort -k1,1 -k2,2n; done | column -t > GC_all_species_total_SCFR_count_across_all_GC_bins
	for species in bonobo borangutan chimpanzee gibbon gorilla sorangutan human; do awk -v a="$species" '{counts[$1]+=$3} END{for(len in counts)print a,len, counts[len]}' ${species}_AT_bins.out | sort -k1,1 -k2,2n; done | column -t > AT_all_species_total_SCFR_count_across_all_AT_bins

#Print %GC or %AT & SCFR length bins where peak (highest SCFR count) is found
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon; do max=$(grep $species GC_all_species_total_SCFR_count_across_all_length_bins | sort -k3,3nr | head -1); echo $species $max; done | column -t
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon; do max=$(grep $species AT_all_species_total_SCFR_count_across_all_length_bins | sort -k3,3nr | head -1); echo $species $max; done | column -t
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon; do max=$(grep $species GC_all_species_total_SCFR_count_across_all_GC_bins | sort -k3,3nr | head -1); echo $species $max; done | column -t
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon; do max=$(grep $species AT_all_species_total_SCFR_count_across_all_AT_bins | sort -k3,3nr | head -1); echo $species $max; done | column -t

#Print number of genes overlapping long SCFRs in each species
	wc -l *_genes_of_interest.txt  | sort -k1n
#Print number of MUC gene members per species
	for i in $(ls | grep "_genes_of_interest.txt"); do j=$(grep MUC -i $i -c); echo $i $j; unset j; done | column -t

#######################################################################################################################################################################################################################################################################################################
#7. Gene deserts

#7.1. Identify gene deserts (1 min)
	mkdir $base_directory/SCFR-main/gene_deserts
	cd $base_directory/SCFR-main
	start_time=$(date +%s)
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	(
	echo $species
	bed=$(find genes/$species/ -name "GCF_*.bed")
	#SCFRs from both intronic & intergenic region
	python3 scripts/gene_desert_finder.py $bed --z 2 --out gene_deserts/$species"_intronic_intergenic"
	#SCFRs from only intergenic region
	python3 my_scripts/gene_desert_finder_intergenic.py $bed --z 2 --out gene_deserts/$species"_only_intergenic"
	) &
	done
	wait
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#7.2. Plot length stats of gene deserts

#Get genome sizes
	mkdir -p genome_sizes
	for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
	do
	    out=genome_sizes/${species}.genome
	    > "$out"
	    for fa in chrs/$species/*.fasta
	    do
	        chr=$(basename "$fa" .fasta)
	        len=$(awk '/^>/ {next} {L+=length($0)} END {print L}' "$fa")
	        echo -e "${chr}\t${len}" >> "$out"
	    done
	    echo "Created $out"
	done
#Get gene desert bed file
	for f in gene_deserts/*_only_intergenic_gene_deserts.tsv; do
	species=$(basename "$f" | cut -d'_' -f1)
	awk 'BEGIN{OFS="\t"} NR>1 {print $1,$2,$3}' "$f" > gene_deserts/${species}_only_intergenic_gene_deserts.bed
	done
#Get SCFR in bed file
	for f in SCFR_all/*SCFR_all.out; do
	species=$(basename "$f" | cut -d'_' -f1)
	awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$f" > SCFR_all/${species}_SCFR_all.bed
	done

#7.3. Plot length stats
	#Gene deserts
	cd $base_directory/SCFR-main/
	python3 my_scripts/compute_desert_stats.py gene_deserts/*_only_intergenic_gene_deserts.bed
	sed -e 's/_only_intergenic_gene_deserts.bed//g' -e 's/sorangutan/Sumatran orangutan/g' -e 's/borangutan/Bornean orangutan/g' desert_summary.tsv -i
	sed '/species/! s/^\(.\)/\U\1/' desert_summary.tsv -i
	Rscript my_scripts/plot_gene_desert_stats.r desert_summary.tsv gene_deserts/all_speices_length_distribution_gene_deserts.pdf
	
	#Intergenic
	cd $base_directory/SCFR-main/
	for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
	do
	cd gene_deserts
	rf=$(ls | grep "_only_intergenic_intergenic_regions.tsv" | grep $species)
	awk '!/zscore/ {print$3-$2}' $rf | python3 $base_directory/SCFR-main/my_scripts/get_stats.py | egrep "Count|^Minimum|^Maximum|^Mean|^Median|^Q1|^Q3" | awk -F ":" '{print$NF}' | tr -d " ," | paste -s -d " " | sed "s/^/$species /g"
	cd $base_directory/SCFR-main/
	done | sort -k1,1 -k2,2n | sed '1i species N min max mean q1 median q3' | tr " " "\t" > intergenic_region_summary.tsv
	sed -e 's/_only_intergenic_gene_deserts.bed//g' -e 's/sorangutan/Sumatran orangutan/g' -e 's/borangutan/Bornean orangutan/g' intergenic_region_summary.tsv -i
	sed '/species/! s/^\(.\)/\U\1/' intergenic_region_summary.tsv -i
	Rscript my_scripts/plot_intergenic_stats.r intergenic_region_summary.tsv gene_deserts/all_speices_length_distribution_intergenic_regions.pdf
	
	mv desert_summary.tsv intergenic_region_summary.tsv gene_deserts/

####################################################################################################################################################################################################################################################################################################################
#7.3. Identify SCFRs in gene deserts

cd $base_directory/SCFR-main/
mkdir gene_deserts/SCFR_overlap_gene_deserts

	for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
	do
	echo ">"$species
	for len in 0 100 500 1000 2500 5000 7500 10000
	do
	mkdir gene_deserts/SCFR_overlap_gene_deserts/$species
	#get overlaps & amount of overlap
	/media/aswin/programs/bedtools2-2.31.1/bin/bedtools intersect -a SCFR_lists/$len/$species"_SCFR_atleast_"$len".out" -b gene_deserts/fishers_test/$species/$species"_intronic_intergenic_gene_deserts.bed" -wo > gene_deserts/SCFR_overlap_gene_deserts/$species/$species"_intronic_intergenic_gene_deserts_scfr_overlaps.out"
	/media/aswin/programs/bedtools2-2.31.1/bin/bedtools intersect -a SCFR_lists/$len/$species"_SCFR_atleast_"$len".out" -b gene_deserts/fishers_test/$species/$species"_only_intergenic_gene_deserts.bed" -wo > gene_deserts/SCFR_overlap_gene_deserts/$species/$species"_only_intergenic_gene_deserts_overlaps.out"
	done
	done

#Summary of SCFRs in gene deserts at different length thresholds (16m3.864s mins)
	cd $base_directory/SCFR-main/
	time for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
	do
	for o in $(find gene_deserts/SCFR_overlap_gene_deserts/$species -name "*_overlaps.out")
	do
	len=$(echo $o | awk -F "/" '{print$NF}' | cut -f2 -d "_")
	#stats=$(awk '{print$NF}' gene_deserts/SCFR_overlap_gene_deserts/$species/$species"_"$len"_only_intergenic_gene_deserts_overlaps.out" | ministat -n | tail -1 | sed 's/^x //g' | sed 's/[ ]\+/ /g' | sed 's/^[ ]\+//g' | sed "s/^/$species $len /g")
	stats=$(awk '{print$NF}' gene_deserts/SCFR_overlap_gene_deserts/$species/$species"_"$len"_only_intergenic_gene_deserts_overlaps.out" | python3 $base_directory/SCFR-main/my_scripts/get_stats.py | egrep "Count|^Minimum|^Maximum|^Mean|^Median|^Q1|^Q3" | awk -F ":" '{print$NF}' | tr -d " ," | paste -s -d " " | sed "s/^/$species $len /g")
	if [[ $(echo "$stats" | awk '{print NF}') -lt 3 ]]; then stats=$(echo $species $len "0 0 0 0 0 0 0"); else :; fi
	echo $stats
	unset len stats
	done
	unset o
	done | sort -k1,1 -k2,2n | sed '1i Species Length_threshold overlapping_SCFR_count Min Max Mean Q1 median Q3' | column -t > gene_deserts/SCFR_overlap_gene_deserts/all_species_scfr_gene_deserts_overlap_summary

#Plot overlap length stats
	cd $base_directory/SCFR-main/gene_deserts/SCFR_overlap_gene_deserts
	for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
	do
	echo ">"$species
	grep "$species" all_species_scfr_gene_deserts_overlap_summary | awk '!($1="")' | sed 's/^[ ]\+//g' | sed '1i Length_threshold N min max mean q1 median q3' | tr " " "\t" > summary_"$species".tsv
	if  [[ "$species" == "human" ]]; then
	Rscript $base_directory/SCFR-main/my_scripts/Figure_2/plot_overlap_stats.R summary_"$species".tsv "$species"_scfr_gene_deserts_overlap_stats.pdf $species
	elif [[ "$species" == "gorilla" ]]; then
	Rscript $base_directory/SCFR-main/my_scripts/Figure_2/plot_overlap_stats_log_transformed.R summary_"$species".tsv "$species"_scfr_gene_deserts_overlap_stats.pdf $species
	else
	Rscript $base_directory/SCFR-main/my_scripts/Figure_2/plot_overlap_stats_extended.R summary_"$species".tsv "$species"_scfr_gene_deserts_overlap_stats.pdf $species
	fi
	rm summary_"$species".tsv
	done

#Plot percentage coverages
	cd $base_directory/SCFR-main
	for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
	do
	gs=$(awk '{a+=$2} END{print a}' genome_sizes/"$species".genome)
	gds=$(bedtools merge -i gene_deserts/"$species"_only_intergenic_gene_deserts.bed | awk '{print $3-$2}' | awk '{a+=$1} END{print a}')
	pgbygd=$(calc $gds / $gs | tr -d "~" | awk '{print$1*100}' | awk '{$NF+=0}1' CONVFMT="%.2f")
	msp0=$(readlink -f genes/$species/"$species"_scfr_atleast_0_merged.bed)
	msp5000=$(readlink -f genes/$species/"$species"_scfr_atleast_5000_merged.bed)
	gdp=$(readlink -f gene_deserts/"$species"_only_intergenic_gene_deserts.bed)
	ovl0=$(/media/aswin/programs/bedtools2-2.31.1/bin/bedtools intersect -a $msp0 -b $gdp -wo | awk '{a+=$NF} END{print a}')
	ovl5000=$(/media/aswin/programs/bedtools2-2.31.1/bin/bedtools intersect -a $msp5000 -b $gdp -wo | awk '{a+=$NF} END{print a}')
	pgdbys0=$(calc $ovl0 / $gds | tr -d "~" | awk '{print$1*100}' | awk '{$NF+=0}1' CONVFMT="%.2f")
	pgdbys5000=$(calc $ovl5000 / $gds | tr -d "~" | awk '{print$1*100}' | awk '{$NF+=0}1' CONVFMT="%.2f")
	echo $species $gs $gds $ovl0 $ovl5000 $pgbygd $pgdbys0 $pgdbys5000
	unset gs gds pgbygd msp0 msp5000 gdp ovl0 ovl5000 pgdbys0 pgdbys5000
	done | sed '1i species genome_size gene_desert_size overlap_with_unfiltered_scfr overlap_with_5000_scfr percent_genome_by_gene_deserts percent_gene_deserts_by_unfiltered_scfr percent_gene_deserts_by_5000_scfr' | sed 's/ /\t/g' > gene_deserts/gene_deserts_percent_covered_summary.tsv

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#7.4. Identify proto-genes or gene-like regions within SCFRs overlapping gene deserts

#7.4.1. Get unique ORF fasta of overlapping SCFRs (12.1167 mins)

cd $base_directory/SCFR-main/
start_time=$(date +%s)
for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
do
echo ">"$species
cd gene_deserts/SCFR_overlap_gene_deserts/$species
mkdir SCFR_fasta
#for each length thresholds
for o in $(find . -name "*_overlaps.out" | grep "5000_only" | grep -v "fourier")
do
#Get fasta only of the overlap file is non-empty
	if [[ -s "$o" ]]; then
	len=$(echo $o | awk -F "/" '{print$NF}' | cut -f2 -d "_")
	#for each scfr
	for scfr in $(awk '{if($4~"-") print$4"::"$1":"$2"-"$3"("substr($4,1,1)")"; else print$4"::"$1":"$2"-"$3"(+)"}' $o)
	do
	chr=$(echo $scfr | cut -f3 -d ":")
	name=$(echo $scfr | sed 's/^-/minus_/g' | cut -f1 -d "(" | tr ":-" "_" | sed 's/__/_/g')
	#myfasta -mfp $base_directory/SCFR-main/PCA/$species/$len/with_coding_region/$chr".fasta" "$scfr" > SCFR_fasta/$name".fa"
	myfasta -mfp $base_directory/SCFR-main/PCA/$species/$len/with_coding_region/$chr".fasta" "$scfr"
	unset chr names
	done > SCFR_fasta/$species"_"$len"_overlapping_scfrs.fa"
#Get canonical ORFs
	ORFfinder -in SCFR_fasta/$species"_"$len"_overlapping_scfrs.fa" -n false -s 0 -ml 600 | myfasta -comb > SCFR_fasta/$species"_"$len"_overlapping_scfrs_canonical_orf.fa"
#Get unique sequences
	if [[ $species == "borangutan" || $species == "sorangutan" ]]; then it="0.65"; else it="0.70"; fi
	$base_directory/SCFR-main/my_scripts/cd_hit_find_unique_sequences.sh SCFR_fasta/$species"_"$len"_overlapping_scfrs_canonical_orf.fa" SCFR_fasta/$species"_"$len"_overlapping_scfrs_canonical_orf_unique.fa" $it
	rm SCFR_fasta/$species"_"$len"_overlapping_scfrs_canonical_orf_unique.log" 
#Get non-canonical ORFs
	ORFfinder -in SCFR_fasta/$species"_"$len"_overlapping_scfrs.fa" -n false -s 1 -ml 600 | myfasta -comb > SCFR_fasta/$species"_"$len"_overlapping_scfrs_non_canonical_orf.fa"
#Get unique sequences
	$base_directory/SCFR-main/my_scripts/cd_hit_find_unique_sequences.sh SCFR_fasta/$species"_"$len"_overlapping_scfrs_non_canonical_orf.fa" SCFR_fasta/$species"_"$len"_overlapping_scfrs_non_canonical_orf_unique.fa" $it
	rm SCFR_fasta/$species"_"$len"_overlapping_scfrs_non_canonical_orf_unique.log"
#Remove canonical ORFs from non-canonical ORFs
	python3 $base_directory/SCFR-main/my_scripts/subtract_fasta.py SCFR_fasta/$species"_"$len"_overlapping_scfrs_canonical_orf_unique.fa" SCFR_fasta/$species"_"$len"_overlapping_scfrs_non_canonical_orf_unique.fa" > SCFR_fasta/$species"_"$len"_overlapping_scfrs_non_canonical_orf_unique_filtered.fa"
unset len scfr it
else :
fi
done
unset o
cd $base_directory/SCFR-main/
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#Summary: print number of overlaps & fasta & ORFs
	cd $base_directory/SCFR-main/gene_deserts/SCFR_overlap_gene_deserts
	for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
	do
	i0=$(wc -l "$species"/"$species"_5000_only_intergenic_gene_deserts_overlaps.out | awk '{print$1}')
	i1=$(grep ">" -c "$species"/SCFR_fasta/"$species"_5000_overlapping_scfrs.fa)
	i2=$(grep ">" -c "$species"/SCFR_fasta/"$species"_5000_overlapping_scfrs_canonical_orf_unique.fa)
	i3=$(grep ">" -c "$species"/SCFR_fasta/"$species"_5000_overlapping_scfrs_non_canonical_orf_unique_filtered.fa)
	echo $species $i0 $i1 $i2 $i3 | awk '{print$0,$4+$5}'
	unset i0 i1 i2 i3
	done | sed '1i species #_scfr_gene_deserts_overlaps #_scfr_fasta #_canonical_ORFs #_non_canonical_ORFs Total_ORFs' | tr " " "\t" > all_species_chosen_scfr_overlap_gene_deserts.tsv

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#7.4.2. Run nr blast locally

#Download nr blastdb
	cd /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison
	time ./v5_download_nr_database.sh &> stdout_v5_download_nr_database
	#cd /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/v5_nr_blastdb
	#time blastdbcmd -db /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/v5_nr_blastdb/nr -entry all -outfmt "%f" -out extracted_nr.fasta
	#time /media/aswin/programs/diamond makedb --in extracted_nr.fasta --db nr_diamond.dmnd --threads 32
	#/media/aswin/programs/diamond makedb --in nr.fasta --db nr_tax_db.dmnd --taxon-map prot.accession2taxid.gz --taxon-nodes taxid.map --threads 32

#For 6 species: 2451.38 mins / 40.8564 hours / 1.70235 days
#For gorilla: 2074.65 mins// 34.5775 hours / 1.44073 days

#Run blastp
	cd $base_directory/SCFR-main/
	start_time=$(date +%s)
	#Don't use gorilla here
	for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
	do
	echo ">"$species
	cd gene_deserts/SCFR_overlap_gene_deserts/
	mkdir $species/SCFR_fasta/nr_blast
	for path in $(find $species/SCFR_fasta -name "*_5000_overlapping_scfrs_canonical_orf_unique.fa" -type f)
	do
	can=$(echo $path)
	scfrcan=$(echo $can | awk -F "/" '{print$NF}')
	noncan=$(echo $can | sed 's/_canonical_orf_unique.fa/_non_canonical_orf_unique_filtered.fa/g')
	scfrnoncan=$(echo $noncan | awk -F "/" '{print$NF}')
	echo " - "$scfrcan
	#Run blast on canonical ORFs
	time /media/aswin/programs/ncbi-blast-2.16.0+/bin/blastp -task blastp-fast -query $can -db /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/v5_nr_blastdb/nr -evalue 0.0001 -max_target_seqs 20 -max_hsps 1 -qcov_hsp_perc 70 -num_threads 32 -outfmt "6 stitle qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand" \
	 | sed '1i Subject_Title\tQuery\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\n' > $species/SCFR_fasta/nr_blast/$scfrcan".outfmt6"
	#Run blast on non-canonical ORFs
	#time /media/aswin/programs/ncbi-blast-2.16.0+/bin/blastp -task blastp-fast -query $noncan -db /media/aswin/gene_loss/APOBEC1/bird_mammal_A1_comparison/v5_nr_blastdb/nr -evalue 0.0001 -max_target_seqs 20 -max_hsps 1 -qcov_hsp_perc 70 -num_threads 32 -outfmt "6 stitle qseqid sseqid qlen length qstart qend sstart send evalue bitscore score qcovs qcovhsp pident nident mismatch gaps sstrand" \
	 #| sed '1i Subject_Title\tQuery\tSubject\tQuery_length\tAlignment_length\tQ_start\tQ_end\tS_start\tS_end\tE_value\tBit_score\tRaw_score\t%_Query_covered_per_sub\t%_Query_covered_per_hsp\t%_ident\tMatches\tMismatches\tGaps\tStrand\n' > $species/SCFR_fasta/nr_blast/$scfrnoncan".outfmt6"
	unset can scfrcan noncan scfrnoncan
	done
	cd $base_directory/SCFR-main/
	done
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#Time taken: 
	#borangutan_5000_overlapping_scfrs_canonical_orf_unique.fa (985m37.837s), gibbon_5000_overlapping_scfrs_canonical_orf_unique.fa (49m0.796s), bonobo_5000_overlapping_scfrs_canonical_orf_unique.fa (150m23.588s)
	#chimpanzee_5000_overlapping_scfrs_canonical_orf_unique.fa (102m24.520s), human_5000_overlapping_scfrs_canonical_orf_unique.fa (104m32.931s)


#7.4.3. Identify proteins homologous to ORFs identified from SCFRs overlapping gene deserts

cd $base_directory/SCFR-main
start_time=$(date +%s)
for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
do
echo ">"$species
#All unique hits & their data
	cd gene_deserts/SCFR_overlap_gene_deserts/$species/SCFR_fasta/nr_blast/
	for hit in $(awk -F "\t" 'NR>1{print$3}' "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6 | awk NF | sort | uniq)
	do
	awk -F "\t" -v a="$hit" '$3==a' "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6 | sed 's/ /_/g' > temp_hit
	nh=$(wc -l temp_hit | awk '{print$1}')
	ev=$(sort -k10,10g temp_hit | awk 'NR==1{print$10}')
	qc=$(sort -k14,14nr temp_hit | awk 'NR==1{print$14}')
	pi=$(sort -k15,15nr temp_hit | awk 'NR==1{print$15}')
	echo $hit $nh $ev $qc $pi
	unset nh ev qc pi
	rm temp_hit
	done | sed '1i Accession #_hits lowest_evalue highest_qcov highest_pid' | column -t > "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits
	#Top 5 hits per query
	sed 's/ /_/g' "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6 | sort -k2,2V -k10,10g -k15,15nr -k14,14nr | awk 'a[$2]++<5' | awk '{print$3}' | grep -v "Subject" > "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.top5_hits
	#top 5 unique hits with lowest e-value
	grep -v "lowest_evalue" "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits | sort -k3,3g | head -5 | awk '{print$1}' > "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits_top_5_eval
	#top 5 unique hits with highest query coverage per hsp
	grep -v "lowest_evalue" "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits | sort -k4,4nr | head -5 | awk '{print$1}' > "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits_top_5_qcov
	#top 5 unique hits with highest percent identity
	grep -v "lowest_evalue" "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits | sort -k5,5nr | head -5 | awk '{print$1}' > "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits_top_5_pid
	#All unique XP ids
	awk '$1~"XP_"' "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits | awk '{print$1}' > "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits_XP
	cat "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits_top_5_eval "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits_top_5_qcov "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits_top_5_pid "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits_XP "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.top5_hits | awk NF | sort | uniq > "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.final_list
echo " - Total No of final hits:" $(wc -l "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.final_list)
unset hit 
#Fetch metadata from NCBI
	mkdir metadata
	#XP accessions
	for xp in $(awk -F "|" '/XP_/ {print$2}' "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.final_list)
	do
	datasets summary gene accession "$xp" --report gene > metadata/"$xp"_gene.json
	datasets summary gene accession ""$xp --report product > metadata/"$xp"_product.json
	f1=$(cat metadata/"$xp"_product.json | json2xml | xtract -pattern reports -def "-" -element gene_id symbol taxname common_name description | sed 's/ /_/g')
	f2=$(cat metadata/"$xp"_product.json | json2xml | xtract -pattern reports -sep "\n" -element genomic_accession_version | sort -u)
	f3=$(cat metadata/"$xp"_product.json | json2xml | xtract -pattern transcripts -def "-" -element accession_version genomic_range/begin genomic_range/end  | grep "$xp" | awk '{print$(NF-1),$NF}')
	f4=$(cat metadata/"$xp"_product.json | xtract -pattern protein -element accession_version length | grep $xp | awk '{print$2}')
	f5=$(cat metadata/"$xp"_gene.json | json2xml | xtract -pattern gene_ontology -def "-" -sep "," -element biological_processes/name cellular_components/name molecular_functions/name | tr " " "_")
	echo $xp $f1 $f2 $f3 $f4 $f5
	unset f1 f2 f3 f4 f5
	#rm a.json b.json
	done | sed '1i Accession gene_id symbol taxname common_name description chr start end protein_length bio_proc cell_comp mol_fun' | column -t > "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.final_list_xp_summary
#Non XP accessions
	for nonxp in $(awk -F "|" '{if($2!~"XP_") print$2}' "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.final_list)
	do
	esearch -db protein -query "$nonxp" | efetch -format all -mode xml > metadata/"$nonxp".json
	#if [[ ! -s metadata/"$nonxp".json ]]; then esearch -db protein -query "$nonxp" | efetch -format all -mode xml > metadata/"$nonxp".json; else :; fi
	f1=$(cat metadata/"$nonxp".json | xtract -pattern GBSeq -def "-" -element GBSeq_locus GBSeq_definition GBSeq_organism GBSeq_length | tr " " "_")
	f2=$(cat metadata/"$nonxp".json | xtract -pattern GBSeq_references -def "-" -sep "\n" -element GBReference_position | uniq | tr " " "_")
	f3=$(cat metadata/"$nonxp".json | xtract -pattern GBFeature -def "-" -element GBQualifier_value | grep CDD | awk -F "\t" '{print$2}' | sort | uniq -c | awk '{a = $1; $1=""; print $0,"("a")"}' | sed 's/^[ ]\+//g' | tr " " "_" | paste -s -d ",")
	echo $nonxp $f1 $f2 $f3
	unset f1 f2 f3
	#rm d.json
	done | sed '1i Accession Locus definition organism length location CDD' | column -t > "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.final_list_nonxp_summary
unset xp nonxp 
cd $base_directory/SCFR-main
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#8. Fishers test

#8.1. Prepare data
	#Merge SCFRs into genome-wide intervals
	time sort -k1,1 -k2,2n SCFR_all/human_SCFR_all.out > human_SCFR_all_sorted.out
	time bedtools merge -i human_SCFR_all_sorted.out | awk '{print$1,$2,$3}' OFS="\t" > human_SCFR_all_sorted_merged.bed
	#Keep only long SCFRs
	time awk '$3-$2 >= 1000' human_SCFR_all_sorted_merged.bed > human_SCFR_all_sorted_merged_atleast_1kb.bed
	#Fisher's test
	bedtools fisher -a human_SCFR_all_sorted_merged_atleast_1kb.bed -b gene_deserts/human_only_intergenic_gene_deserts.bed -g genome_sizes/human.genome

#8.2. Run fishers test (5 secs)
mkdir $base_directory/SCFR-main/gene_deserts/fishers_test
cd $base_directory/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
	echo " -"$species
	mkdir $base_directory/SCFR-main/gene_deserts/fishers_test/$species
	for tsv in $(find gene_deserts/ -maxdepth 1 -mindepth 1 -name "$species*.tsv")
	do
	prefix=$(echo $tsv | sed 's/.tsv//g' | awk -F "/" '{print$NF}')
	echo " -"$prefix
	awk '!/start/ {print$1,$2,$3,$4}' OFS="\t" $tsv > gene_deserts/fishers_test/$species/$prefix".bed"
	for scfr in $(find $base_directory/SCFR-main/SCFR_lists/ -name "${species}_SCFR_atleast_*.out" | egrep -v "atleast_0.out|atleast_100.out")
	do
	(
	name=$(echo $scfr | awk -F "/" '{print$NF}')
	/media/aswin/programs/bedtools2-2.31.1/bin/bedtools fisher -a $scfr -b gene_deserts/fishers_test/$species/$prefix".bed" -g $base_directory/SCFR-main/genome_sizes/$species".genome" > $base_directory/SCFR-main/gene_deserts/fishers_test/$species/$name"_"$prefix".out"
	) &
	done
	wait
	unset scfr name
	done
	unset tsv prefix
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#8.3. Summary of Fishers test
	cd $base_directory/SCFR-main/
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo " -"$species
	cd $base_directory/SCFR-main/gene_deserts/fishers_test/$species
	for out in $(ls *.out | grep "gene_deserts")
	do
	ab=$(echo $out | sed 's/^.*_SCFR_atleast_//g' | sed 's/\.out/ /g' | sed "s/$species//g" | sed 's/__//g')
	o1=$(cat $out | awk -F ":" '/# Number of/ {print$2}' | paste -s -d " " | sed 's/[ ]\+/ /g' | sed 's/^[ ]\+//g')
	o2=$(cat $out | grep "in -a" | awk -F "|" '{print$2,$3}' | paste -s -d " " | sed 's/[ ]\+/ /g' | sed 's/^[ ]\+//g')
	o3=$(cat $out | tail -1 | paste -s -d " " | sed 's/[ ]\+/ /g' | sed 's/^[ ]\+//g')
	echo $ab $o1 $o2 $o3
	unset ab o1 o2 o3
	done | sed '1i Query DB #Query_intervals #DB_intervals #Overlaps #Possible_intervals in_a_in_b in_a_not_in_b not_in_a_in_b not_in_a_not_in_b left_pvalue right_pvalue two_tail_pvalue ratio' | column -t > $species"_fisher_test_summary"
	unset out
	cd $base_directory/SCFR-main/
	done

#Combine all species summary
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	for fisher in $(find gene_deserts/fishers_test/ -name $species"_fisher_test_summary")
	do
	sed "s/^/$species /g" $fisher | grep -v "Possible_intervals"
	done
	done | sed '1i Species Query DB #Query_intervals #DB_intervals #Overlaps #Possible_intervals in_a_in_b in_a_not_in_b not_in_a_in_b not_in_a_not_in_b left_pvalue right_pvalue two_tail_pvalue ratio' | column -t > $base_directory/SCFR-main/gene_deserts/fishers_test/all_species_summary
	sort -k15,15V -k2,2n all_species_summary

cd $base_directory/SCFR-main
awk '$13<0.05' gene_deserts/fishers_test/all_species_summary > gene_deserts/fishers_test/scfr_significant_higher_overlaps_with_gene_deserts_than_expected

####################################################################################################################################################################################################################################################################################################################
#9. SCFRs at length thresholds

#9.1. Filter with different length cut offs (16.0167 mins)
	cd $base_directory/SCFR-main
	start_time=$(date +%s)
	for win in 0 100 500 1000 2500 5000 7500 10000
	do
	echo ">"$win
	mkdir -p $base_directory/SCFR-main/SCFR_lists/"$win"
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	  do
	  (
	  echo " -"$species
	  cds=$(find genes/$species/ -name "GCF_*.bed")
	  #Filter SCFRs longer than the window
	  awk -v a="$win" '{if($3-$2>=a) print$1,$2,$3,$4}' OFS="\t" SCFR_all/${species}_SCFR_all.out > SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out"
	  #Get SCFRs that don't overlap with coding genes
	  bedtools intersect -v -a SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out" -b $cds | awk '{print$1,$2,$3,$4,$3-$2}' OFS="\t" > SCFR_lists/${win}/${species}"_SCFR_atleast_"$win"_in_non_coding.bed"
	  ) &
	done
	wait
	unset species cds
	done
	end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
	echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#9.2. Get window-wise summary of SCFR count and length & it's overlap with whole genome & cds from all species (99.3667 mins)

cd $base_directory/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
genome_size=$(awk '{a+=$2} END{print a}' $base_directory/SCFR-main/genome_sizes/$species".genome")
#cds
cd $base_directory/SCFR-main/genes/$species
cds=$(ls | grep "GC.*.bed")
#merge cds
	sort -k1,1V -k2,2n $cds | /media/aswin/programs/bedtools2-2.31.1/bin/bedtools merge -i - > $species"_cds_merged.bed"
merged_cds_len=$(awk '{print$3-$2}' $species"_cds_merged.bed" | awk '{a+=$1} END{print a}')
#scfr
total_unfiltered_scfr_len=$(wc -l $base_directory/SCFR-main/SCFR_lists/0/$species"_SCFR_atleast_0.out" | awk '{print$1}')
cd $base_directory/SCFR-main/SCFR_lists
for win in 0 100 500 1000 2500 5000 7500 10000
do
cd $win
#get stats
	all_scfr=$(ls | grep $species"_SCFR_atleast_"$win".out")
	noncoding_scfr=$(ls | grep $species"_SCFR_atleast_"$win"_in_non_coding.bed")
	all_count=$(wc -l $all_scfr | awk '{print$1}')
	noncoding_count=$(wc -l $noncoding_scfr | awk '{print$1}')
	all_len=$(awk '{print$3-$2}' $all_scfr | awk '{a+=$1} END{print a}')
	noncoding_len=$(awk '{print$3-$2}' $noncoding_scfr | awk '{a+=$1} END{print a}')
	all_merged_len=$(/media/aswin/programs/bedtools2-2.31.1/bin/bedtools merge -i $all_scfr | awk '{print$3-$2}' | awk '{a+=$1} END{print a}')
	noncoding_merged_len=$(/media/aswin/programs/bedtools2-2.31.1/bin/bedtools merge -i $noncoding_scfr | awk '{print$3-$2}' | awk '{a+=$1} END{print a}')
#merge csfr
	sort -k1,1V -k2,2n $base_directory/SCFR-main/SCFR_lists/"$win"/"$all_scfr" | /media/aswin/programs/bedtools2-2.31.1/bin/bedtools merge -i - > $base_directory/SCFR-main/genes/$species/$species"_scfr_atleast_"$win"_merged.bed"
#get percentages
	per_unfiltered_scfr_by_filtered=$(calc $all_count / $total_unfiltered_scfr_len | tr -d "~" | awk '{print$1*100}')
	per_genome_by_scfr=$(calc $all_merged_len / $genome_size | tr -d "~" | awk '{print$1*100}')
	merged_scfr_len=$(awk '{print$3-$2}' $base_directory/SCFR-main/genes/$species/$species"_scfr_atleast_"$win"_merged.bed" | awk '{a+=$1} END{print a}')
	scfr_cds_overlap=$(/media/aswin/programs/bedtools2-2.31.1/bin/bedtools intersect -a $base_directory/SCFR-main/genes/$species/$species"_scfr_atleast_"$win"_merged.bed" -b $base_directory/SCFR-main/genes/$species/$species"_cds_merged.bed" | awk '{print$3-$2}' | awk '{a+=$1} END{print a}')
	per_scfr_by_coding=$(calc $scfr_cds_overlap / $all_merged_len | tr -d "~" | awk '{print$1*100}')
scfr=$(echo $species $win $all_count $noncoding_count $all_len $noncoding_len $all_merged_len $noncoding_merged_len $per_unfiltered_scfr_by_filtered $genome_size $per_genome_by_scfr $merged_cds_len $scfr_cds_overlap $per_scfr_by_coding)
echo $scfr
unset all_scfr noncoding_scfr all_count noncoding_count all_len noncoding_len all_merged_len noncoding_merged_len per_unfiltered_scfr_by_filtered per_genome_by_scfr merged_scfr_len scfr_cds_overlap per_scfr_by_coding 
cd $base_directory/SCFR-main/SCFR_lists
done
unset genome_size win cds merged_cds_len total_unfiltered_scfr_len
cd $base_directory/SCFR-main
done | sed '1i Species Window Total_No_SCFR Total_No_noncoding_SCFR Total_length_SCFR Total_length_noncoding_SCFR Merged_length_SCFR Merged_length_noncoding_SCFR Percent_unfiltered_by_filtered_SCFR Genome_size Percent_genome_by_SCFR Merged_CDS_length SCFR_CDS_overlap Percent_SCFR_by_coding' > $base_directory/SCFR-main/SCFR_summaries/window_wise_all_species_scfr_coding_stats_with_zero_window
wait
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e > $base_directory/SCFR-main/SCFR_summaries/runtime_window_wise_all_species_scfr_coding_stats

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#9.3. Plot SCFR count & length stats & percent in genome & cds 
	awk -F "," '{print$1,$2,(($3-$4)/$3*100),($4/$3*100),$3,$9,$11,$14}' window_wise_all_species_scfr_coding_stats.csv | awk 'BEGIN{FS=OFS=","} NR==1{$3="Percent_coding_SCFR"; $4="Percent_noncoding_SCFR"}1' | tr " " "," > filtered_window_wise_all_species_scfr_coding_stats.csv
	awk -F "," 'NR>1{print$1,$2,$3,(($3-$4)/$3*100),($4/$3*100),$9,$11,$14}' window_wise_all_species_scfr_coding_stats.csv | sed '1i Species,Window,Total_No_SCFR,Percent_coding_SCFR_count,Percent_noncoding_SCFR_count,Percent_unfiltered_by_filtered_SCFR,Percent_genome_by_SCFR,Percent_SCFR_by_coding' | sed 's/ /,/g' > filtered_window_wise_all_species_scfr_coding_stats.csv
	
	mkdir SCFR_summaries/plot_all_species_SCFR_count_length_percent_with_genome_and_cds
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo ">"$species
	cat SCFR_summaries/filtered_window_wise_all_species_scfr_coding_stats.csv | egrep "$species|Window" > "temp_"$species".csv"
	Rscript my_scripts/plot_all_species_SCFR_count_length_percent_with_genome_and_cds.R "temp_"$species".csv" SCFR_summaries/plot_all_species_SCFR_count_length_percent_with_genome_and_cds/$species"_plot_all_species_SCFR_count_length_percent_with_genome_and_cds.pdf"
	done

#9.4. Get coding gene length stats of species
	cd $base_directory/SCFR-main/Fourier_analysis/genes
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo ">"$species
	gene_length_stats=$(myfasta -l $species/GC*_cds.fa | awk '{print$NF}' | ministat -n | tail -1 | sed 's/^x //g' | sed 's/[ ]\+/ /g')
	echo $species $gene_length_stats
	unset species gene_length_stats
	done | column -t > gene_length_stats

####################################################################################################################################################################################################################################################################################################################
#11. Discrete Fourier Transform Analysis of SCFR

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#11.1. Run fourier analysis for SCFRs (for 5000, 7500 & 10000 - 72.8667 mins) (for 1000 & 2500 - )

cd $base_directory/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
for win in 1000 2500 5000 7500 10000
do
echo " -"$win
#For SCFRs of atleast a specific length
	mkdir -p Fourier_analysis/"$species"/"$win"/with_coding_region
	for chr in `cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out" | cut -f 1|sort -u`
	do
	cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win".out" | awk -v k=$chr '$1==k{print $0}' | sort -k1,1 -k2n,2 | awk '{if($4~"-") print$0,1,"-"; else print$0,1,"+"}' OFS="\t" \
	| /media/aswin/programs/bedtools2-2.31.1/bin/bedtools getfasta -fi chrs/"$species"/"$chr".fasta -bed stdin -name+ -s > Fourier_analysis/"$species"/"$win"/with_coding_region/"$chr".fasta
	done
#For SCFRs of atleast a specific length and with no overlap with coding region
	mkdir -p Fourier_analysis/"$species"/"$win"/without_coding_region
	for chr in `cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win"_in_non_coding.bed" | cut -f 1|sort -u`
	do
	cat SCFR_lists/${win}/${species}"_SCFR_atleast_"$win"_in_non_coding.bed" | awk -v k=$chr '$1==k{print $0}' | sort -k1,1 -k2n,2 | awk '{if($4~"-") print$0,"-"; else print$0,"+"}' OFS="\t" \
	| /media/aswin/programs/bedtools2-2.31.1/bin/bedtools getfasta -fi chrs/"$species"/"$chr".fasta -bed stdin -name+ -s > Fourier_analysis/"$species"/"$win"/without_coding_region/"$chr".fasta
	done
#Run Fourier analysis
	for type in with_coding_region without_coding_region
	do
	cd Fourier_analysis/"$species"/"$win"/"$type"/
	for scfr in $(find . -mindepth 1 -maxdepth 1 -name "*.fasta" | cut -f2 -d "/")
	do
	echo " -"$scfr
	python3 $base_directory/SCFR-main/Fourier_analysis/scfr_parallel_fft_motif_report_grouped.py -o "output_"$scfr -t 32 $scfr
	done
	unset scfr
	cd $base_directory/SCFR-main
	done
unset chr type
done
unset win
cd $base_directory/SCFR-main
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#11.2. Get chromosome-wise summary of fourier analysis (14.3667 mins)

cd $base_directory/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
  do
  echo ">"$species
for fourier in $(find Fourier_analysis/$species/ -name "output_*.fasta" -type d)
do
path=$(echo $fourier | awk -F "/" '!($NF="")' OFS="/")
folder=$(echo $fourier | awk -F "/" '{print$NF}')
cd $path
python3 $base_directory/SCFR-main/my_scripts/scfr_fourier_chromosome_wise_summary.py $folder --top 3 --cores 32
cd $base_directory/SCFR-main
done
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#Get species-wise summary
cd $base_directory/SCFR-main
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
cd Fourier_analysis/$species
for len in 5000 7500 10000
do
for tsv in $(find $len -name "summary.tsv" -type f)
do
type=$(echo $tsv | cut -f2 -d "/")
chr=$(echo $tsv | cut -f3 -d "/" | cut -f2,3 -d "_" | sed 's/\.fasta//g')
grep -v "Num_Raw_Peaks" $tsv | sed "s/^/$type\t$len\t$chr\t/g"
unset chr type
done
unset tsv
done | column -t > all_length_thresholds_fourier_summary
unset len
cd $base_directory/SCFR-main
done

#Filter SCFRs with 3 periodocity
cd $base_directory/SCFR-main
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
awk '$6~"0.33"' Fourier_analysis/$species/all_length_thresholds_fourier_summary > Fourier_analysis/$species/3_periodicity_scfrs
done

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Plot fourier frequency density

cd $base_directory/SCFR-main
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
cd Fourier_analysis/"$species"/
awk '$1=="without_coding_region"' all_length_thresholds_fourier_summary | sed '1i Coding_status Length_threshold chr SCFR_Name Num_Raw_Peaks Top_Peak_Count Frequencies Magnitudes Periods' | sed 's/[ ]\+/\t/g' > "$species"_without_coding_region_all_length_thresholds_fourier_summary.tsv
awk '$1=="without_coding_region" && $2==5000' all_length_thresholds_fourier_summary | sed '1i Coding_status Length_threshold chr SCFR_Name Num_Raw_Peaks Top_Peak_Count Frequencies Magnitudes Periods' | sed 's/[ ]\+/\t/g' > "$species"_5000_without_coding_region_all_length_thresholds_fourier_summary.tsv
awk '$1=="without_coding_region" && $2==7500' all_length_thresholds_fourier_summary | sed '1i Coding_status Length_threshold chr SCFR_Name Num_Raw_Peaks Top_Peak_Count Frequencies Magnitudes Periods' | sed 's/[ ]\+/\t/g' > "$species"_7500_without_coding_region_all_length_thresholds_fourier_summary.tsv
awk '$1=="without_coding_region" && $2==10000' all_length_thresholds_fourier_summary | sed '1i Coding_status Length_threshold chr SCFR_Name Num_Raw_Peaks Top_Peak_Count Frequencies Magnitudes Periods' | sed 's/[ ]\+/\t/g' > "$species"_10000_without_coding_region_all_length_thresholds_fourier_summary.tsv
awk '$1=="with_coding_region"' all_length_thresholds_fourier_summary | sed '1i Coding_status Length_threshold chr SCFR_Name Num_Raw_Peaks Top_Peak_Count Frequencies Magnitudes Periods' | sed 's/[ ]\+/\t/g' > "$species"_with_coding_region_all_length_thresholds_fourier_summary.tsv
awk '$1=="with_coding_region" && $2==5000' all_length_thresholds_fourier_summary | sed '1i Coding_status Length_threshold chr SCFR_Name Num_Raw_Peaks Top_Peak_Count Frequencies Magnitudes Periods' | sed 's/[ ]\+/\t/g' > "$species"_5000_with_coding_region_all_length_thresholds_fourier_summary.tsv
awk '$1=="with_coding_region" && $2==7500' all_length_thresholds_fourier_summary | sed '1i Coding_status Length_threshold chr SCFR_Name Num_Raw_Peaks Top_Peak_Count Frequencies Magnitudes Periods' | sed 's/[ ]\+/\t/g' > "$species"_7500_with_coding_region_all_length_thresholds_fourier_summary.tsv
awk '$1=="with_coding_region" && $2==10000' all_length_thresholds_fourier_summary | sed '1i Coding_status Length_threshold chr SCFR_Name Num_Raw_Peaks Top_Peak_Count Frequencies Magnitudes Periods' | sed 's/[ ]\+/\t/g' > "$species"_10000_with_coding_region_all_length_thresholds_fourier_summary.tsv
#plot
Rscript $base_directory/SCFR-main/my_scripts/plot_fourier_frequencies.R "$species"_without_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_without_coding_region_all_length_thresholds_fourier_summary.pdf without_coding
Rscript $base_directory/SCFR-main/my_scripts/plot_fourier_frequencies.R "$species"_5000_without_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_5000_without_coding_region_all_length_thresholds_fourier_summary.pdf without_coding
Rscript $base_directory/SCFR-main/my_scripts/plot_fourier_frequencies.R "$species"_7500_without_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_7500_without_coding_region_all_length_thresholds_fourier_summary.pdf without_coding
Rscript $base_directory/SCFR-main/my_scripts/plot_fourier_frequencies.R "$species"_10000_without_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_10000_without_coding_region_all_length_thresholds_fourier_summary.pdf without_coding
#plot 
Rscript $base_directory/SCFR-main/my_scripts/plot_fourier_frequencies.R "$species"_with_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_with_coding_region_all_length_thresholds_fourier_summary.pdf all
Rscript $base_directory/SCFR-main/my_scripts/plot_fourier_frequencies.R "$species"_5000_with_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_5000_with_coding_region_all_length_thresholds_fourier_summary.pdf all
Rscript $base_directory/SCFR-main/my_scripts/plot_fourier_frequencies.R "$species"_7500_with_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_7500_with_coding_region_all_length_thresholds_fourier_summary.pdf all
Rscript $base_directory/SCFR-main/my_scripts/plot_fourier_frequencies.R "$species"_10000_with_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_10000_with_coding_region_all_length_thresholds_fourier_summary.pdf all
cd $base_directory/SCFR-main
done

cd $base_directory/SCFR-main
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">"$species
mkdir -p Fourier_analysis/"$species"/table
cd Fourier_analysis/"$species"/table
Rscript $base_directory/SCFR-main/my_scripts/table_fourier_frequencies.R ../"$species"_without_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_without_coding_region_all_length_thresholds_fourier_summary.tsv
Rscript $base_directory/SCFR-main/my_scripts/table_fourier_frequencies.R ../"$species"_5000_without_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_5000_without_coding_region_all_length_thresholds_fourier_summary.tsv
Rscript $base_directory/SCFR-main/my_scripts/table_fourier_frequencies.R ../"$species"_7500_without_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_7500_without_coding_region_all_length_thresholds_fourier_summary.tsv
Rscript $base_directory/SCFR-main/my_scripts/table_fourier_frequencies.R ../"$species"_10000_without_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_10000_without_coding_region_all_length_thresholds_fourier_summary.tsv
#plot 
Rscript $base_directory/SCFR-main/my_scripts/table_fourier_frequencies.R ../"$species"_with_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_with_coding_region_all_length_thresholds_fourier_summary.tsv
Rscript $base_directory/SCFR-main/my_scripts/table_fourier_frequencies.R ../"$species"_5000_with_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_5000_with_coding_region_all_length_thresholds_fourier_summary.tsv
Rscript $base_directory/SCFR-main/my_scripts/table_fourier_frequencies.R ../"$species"_7500_with_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_7500_with_coding_region_all_length_thresholds_fourier_summary.tsv
Rscript $base_directory/SCFR-main/my_scripts/table_fourier_frequencies.R ../"$species"_10000_with_coding_region_all_length_thresholds_fourier_summary.tsv "$species"_10000_with_coding_region_all_length_thresholds_fourier_summary.tsv
for t in $(ls *region_all_length_thresholds_fourier_summary.tsv)
do
t1=$(echo $t | awk -F "_" '{print$1,$2,$3}')
awk -v t1="$t1" '{$1=t1; print}' $t | sed '1d'
done | awk '{if($3=="with") $3="all"; else if($2=="with") $3="all"; else if($2=="without") $3="without_coding"; else if($3=="without") $3="without_coding"; print}' \
| sed '1i species length_threshold SCFR_coding_content highest_peak second_highest_peak peak_difference dip modal_of_distribution' | sed 's/[ ]\+/\t/g' > "$species"_final_fourier_summary.tsv
cd $base_directory/SCFR-main
done
find Fourier_analysis/ -mindepth 3 -maxdepth 3 -name "*_final_fourier_summary.tsv" -type f | xargs cat | awk '!a[$0]++' > Fourier_analysis/fourier_frequency_table.tsv

mkdir -p $base_directory/SCFR-main/shreya/Fourier/
find Fourier_analysis/ -mindepth 2 -maxdepth 2 -name "*_region_all_length_thresholds_fourier_summary.pdf" -type f | xargs -n1 sh -c 'cp $0 $base_directory/SCFR-main/shreya/Fourier/'

#Look at fourier periodicity in SCFRs located inside gene deserts

####################################################################################################################################################################################################################################################################################################################
#12. Fourier analysis of genes

#12.1. Get GTF associated with genome
	mkdir $base_directory/SCFR-main/Fourier_analysis/genes
	cd $base_directory/SCFR-main/Fourier_analysis/genes
	cat $base_directory/SCFR-main/QC/genome_accessions | xargs -n2 bash -c 'paste <(echo $1 $0) <(datasets summary genome accession $0 --as-json-lines | json2xml | xtract -pattern root -def "-" -element accession)' | column -t > genome_gtf_accessions
	for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo ">"$species

#Download cds from ncbi (21.9667 mins)
cd $base_directory/SCFR-main/Fourier_analysis/genes
start_time=$(date +%s)
while read i
	do
	species=$(echo $i | awk '{print$1}')
	gtf=$(echo $i | awk '{print$NF}')
	echo ">"$species $gtf
	time datasets download genome accession $gtf --include cds --filename $gtf".zip"
	unzip -o $gtf".zip" -d cds
	mkdir $species
	mv cds/ncbi_dataset/data/$gtf/cds_from_genomic.fna $species/$gtf"_cds.fa"
	unset species gtf
done < <(grep borangutan genome_gtf_accessions)
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

rm -r cds

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#12.2. Check gene cds length stats

#The code for fourier don't work for very short sequences; hence check length stats of cds
cd $base_directory/SCFR-main/Fourier_analysis/genes
time for seq in $(find . -mindepth 2 -maxdepth 2 -name "GCF_*_cds.fa")
do
id=$(echo $seq | awk -F / '{print$2,$3}' | sed 's/_cds.fa//g')
stat=$(myfasta -l $seq | awk '{print$NF}' | ministat -n | tail -1 | sed 's/^x//g' | sed 's/[ ]\+/ /g' | sed 's/^ //g')
echo $id $stat
unset id stat
done | sed '1i Species annotation N Min Max Median Avg Stddev' | column -t > annotation_length_stats

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#12.3. Run fourier analysis of genes (474.25)

cd $base_directory/SCFR-main
start_time=$(date +%s)
for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo ">"$species
	cd Fourier_analysis/genes/$species
	gene=$(ls | grep "GC.*_cds.fa")
	if [[ $species == human ]]; then seqtk seq -L 30 $gene > $gene"_filtered.fa"; gene=$(ls | grep "_filtered.fa"); else :; fi
	python3 $base_directory/SCFR-main/Fourier_analysis/scfr_parallel_fft_motif_report_grouped.py -o "output_"$gene -t 32 $gene
	cd $base_directory/SCFR-main
	unset gene
done
end_time=$(date +%s) && elapsed_time=$((end_time - start_time))
echo -e "\n Total time taken:" && echo $elapsed_time | awk '{print"-days:",$NF/60/60/24,"\n","-hours:",$NF/60/60,"\n","-mins:",$NF/60,"\n","-secs:",$1}' | column -t | sed 's/^/   /g' && echo -e

#Summary chromosome-wise summary (113 mins)
cd $base_directory/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
	do
	echo ">"$species
	cd Fourier_analysis/genes/$species
	folder=$(ls | grep "output_GCF_.*_cds.fa")
	python3 $base_directory/SCFR-main/my_scripts/scfr_fourier_chromosome_wise_summary.py $folder --top 3 --cores 32
	cd $base_directory/SCFR-main
done

#Plot frequency Vs amplitude
cd $base_directory/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
input=$(find Fourier_analysis/genes/$species/output_GCF_*_cds.fa*/chromosome_wise_summary/ -name "summary.tsv")
output="Fourier_analysis/genes/$species/${species}_all_genes_freq_mag.tsv"
output2="Fourier_analysis/genes/$species/${species}_all_genes_positive_freq_mag.tsv"
echo -e ">"$species "\n -input: "$input "\n -output:"$output
my_scripts/extract_fourier_freq_mag_from_summary.sh $input $output
awk '$2!~"-"' $output > $output2
sed 's/^SCFR_Name/Gene/g' $output -i
sed 's/^SCFR_Name/Gene/g' $output2 -i
Rscript $base_directory/SCFR-main/my_scripts/plot_heatmap_fourier_peaks_genes.R $output2 "$species"_all_genes_positive_freq_mag.pdf
unset input output output2
done

Rscript $base_directory/SCFR-main/my_scripts/plot_heatmap_fourier_peaks_genes.R human_all_genes_positive_freq_mag.tsv human_all_genes_positive_freq_mag.pdf

cd $base_directory/SCFR-main/Fourier_analysis/genes
time find . -maxdepth 2 -mindepth 2 -name "*_all_genes_positive_freq_mag.tsv" | xargs -n1 sh -c 'grep -v "Magnitude" $0' | sed '1i Gene Frequency Magnitude' | sed 's/ /\t/g' > all_species_all_genes_positive_freq_mag.tsv
Rscript $base_directory/SCFR-main/my_scripts/plot_heatmap_fourier_peaks_genes.R all_species_all_genes_positive_freq_mag.tsv all_species_all_genes_positive_freq_mag.pdf

####################################################################################################################################################################################################################################################################################################################
#13. Identify proto-genes

#Get fourier peak freq of SCFRs overlapping gene deserts
cd $base_directory/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
  do
grep -f <(awk '{if($4~"-") print$1,$2,$3,"frame_"$4; else print$1,$2,$3,"frame"$4}' OFS="_" gene_deserts/SCFR_overlap_gene_deserts/$species/$species"_5000_only_intergenic_gene_deserts_overlaps.out" | tr -d "-" | tr "." "_") <(awk '$1=="with_coding_region" && $2=="5000"' Fourier_analysis/$species/all_length_thresholds_fourier_summary) > gene_deserts/SCFR_overlap_gene_deserts/$species/$species"_5000_only_intergenic_gene_deserts_overlaps_fourier.out"
done

#Run Fourier analysis in ORFs
cd $base_directory/SCFR-main
time for species in human bonobo chimpanzee gorilla borangutan sorangutan gibbon
do
echo ">" $species
cd gene_deserts/SCFR_overlap_gene_deserts/"$species"/SCFR_fasta
#Run fourier
python3 $base_directory/SCFR-main/my_scripts/orf_parallel_fft_motif_report_grouped.py "$species"_5000_overlapping_scfrs_canonical_orf_unique.fa -t 32 -o output_"$species"_5000_overlapping_scfrs_canonical_orf_unique.fa
python3 $base_directory/SCFR-main/my_scripts/orf_parallel_fft_motif_report_grouped.py "$species"_5000_overlapping_scfrs_non_canonical_orf_unique_filtered.fa -t 32 -o output_"$species"_5000_overlapping_scfrs_non_canonical_orf_unique_filtered.fa
#Get chromosomr-wise summary
_composite_exon_filteredpython3 $base_directory/SCFR-main/my_scripts/scfr_fourier_chromosome_wise_summary.py output_"$species"_5000_overlapping_scfrs_canonical_orf_unique.fa --top 3 --cores 32
python3 $base_directory/SCFR-main/my_scripts/scfr_fourier_chromosome_wise_summary.py output_"$species"_5000_overlapping_scfrs_non_canonical_orf_unique_filtered.fa --top 3 --cores 32
cd $base_directory/SCFR-main
done

#Summary
cd $base_directory/SCFR-main
for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
do
grep -v "Top_Peak_Count" gene_deserts/SCFR_overlap_gene_deserts/"$species"/SCFR_fasta/output_"$species"_5000_overlapping_scfrs_canonical_orf_unique.fa/chromosome_wise_summary/summary.tsv | sed "s/^/$species /g"
done | sed '1i Species SCFR_Name Num_Raw_Peaks Top_Peak_Count Frequencies Magnitudes Periods' | tr " " "\t" > gene_deserts/SCFR_overlap_gene_deserts/all_species_5000_overlapping_scfrs_canonical_orf_unique_fourier.tsv

for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
do
grep -v "Top_Peak_Count" gene_deserts/SCFR_overlap_gene_deserts/"$species"/SCFR_fasta/output_"$species"_5000_overlapping_scfrs_non_canonical_orf_unique_filtered.fa/chromosome_wise_summary/summary.tsv | sed "s/^/$species /g"
done | sed '1i Species SCFR_Name Num_Raw_Peaks Top_Peak_Count Frequencies Magnitudes Periods' | tr " " "\t" > gene_deserts/SCFR_overlap_gene_deserts/all_species_5000_overlapping_scfrs_non_canonical_orf_unique_filtered_fourier.tsv

awk -F "\t" '$5~/0\.33/' gene_deserts/SCFR_overlap_gene_deserts/all_species_5000_overlapping_scfrs_canonical_orf_unique_fourier.tsv
awk -F "\t" '$5~/0\.33/' gene_deserts/SCFR_overlap_gene_deserts/all_species_5000_overlapping_scfrs_non_canonical_orf_unique_filtered_fourier.tsv

for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan
do
all=$(wc -l gene_deserts/SCFR_overlap_gene_deserts/"$species"/SCFR_fasta/nr_blast/"$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits | awk '{print$1}')
xp=$(awk '$1~"XP_"' gene_deserts/SCFR_overlap_gene_deserts/"$species"/SCFR_fasta/nr_blast/"$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits | wc -l | awk '{print$1}')
echo $species $all $xp | awk '{print$0,($3/$2)*100}'
#sort -k3,3g gene_deserts/SCFR_overlap_gene_deserts/"$species"/SCFR_fasta/nr_blast/"$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.all_hits | grep "XP_" --color=always -z
unset all xp
done | less -SNR

for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan; do echo ">"$species; awk 'NR>1{print$6}' gene_deserts/SCFR_overlap_gene_deserts/"$species"/SCFR_fasta/nr_blast/"$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.final_list_xp_summary | sort | uniq -c; done | less
for species in human chimpanzee gorilla bonobo gibbon borangutan sorangutan; do echo ">"$species; awk 'NR>1{print$6}' gene_deserts/SCFR_overlap_gene_deserts/"$species"/SCFR_fasta/nr_blast/"$species"_5000_overlapping_scfrs_canonical_orf_unique.fa.outfmt6.final_list_xp_summary | awk NF | sed 's/uncharacterized_LOC.*/uncharacterized_LOC/g' | sort | uniq -c | sort -k1,1nr ; done | egrep ">|serine"


awk '{if($4~"-") print$1,$10,$11,$5"_exitron","1","-"; else print$1,$10,$11,$5"_exitron","1","+"}' OFS="\t" "$species"_results.exitron_candidates.txt | grep -v "#chrom" > "$species"_results.exitron_candidates.bed
time while read e
do
e1=$(echo $e | awk '{print$1}')
fa=$(find $base_directory/SCFR-main/chrs/$species/ -name "$e1*.fasta")
/media/aswin/programs/bedtools2-2.31.1/bin/bedtools getfasta -fi "$fa" -bed <(echo "$e") -name+ -s | sed '/^>/ s/(-)/_minus/g' | sed '/^>/ s/(+)/_plus/g' | tr "-" "_"
unset e1 fa
done < "$species"_results.exitron_candidates.bed > "$species"_results.exitron_candidates.fa


####################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################
