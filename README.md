# SCFR 

This GitHub repository contains the data and associated code for the paper **"The Landscape of Stop Codon-Free Regions in the Human Genome: A Reservoir of Proto-Genes"** 

Aswin S Soman<sup>1</sup> , G Shreyasree<sup>2</sup> , Aishwarya Dwivedi<sup>3</sup>, Gayathri S Pramod, Chirayu Sakarkar, Dipannita Bhattacharya<sup>4</sup>, Nagarjun Vijay <sup>1*</sup>

Computational Evolutionary Genomics Lab, Department of Biological Sciences, IISER Bhopal, Bhauri, Madhya Pradesh, India

Corresponding Author: nagarjun@iiserb.ac.in

# Folder Structure 
<ul>
  <li> <b>Supplementary Data:</b> 
    <ul>
      <li> Length_threshold_PCA_kmeans: contains data related to PCA loadings and the plots with SCFRs indicated by clade, cluster, coding status, GC content and repeat composition. Plots for every species are present within their respective folder. </li>
      <li> Exon_shadow: contains data of all exon shadows for 7 primate species. </li>
      <li> Gene deserts: contains data related to intergenic and intronic gene deserts, and intergenic along with distribution plot. Contains a subfolder 'chromosome_wide_plots' that has chromosome-wide analysis of gene desert distribution. </li>
      <li> QC: contains genome accessions and genome metadata </li>
      <li> genome_reports: contains chromosome-wise summary data of each species' genome </li>
      <li> genome_sizes: contains data concerning chromosome-wise size in each genome </li>
      <li> nr_blastp: contains results of running blastp using gene desert derived SCFR-ORFs against nr database </li>      
    </ul>
  </li>
  <li> Scripts: </li>
  <ul>
    <li> fourier_analysis: contains all scripts required for conducting discrete Fourier analysis </li>
    <li> gene_deserts: Python scripts for the identification of gene deserts </li>
    <li> PCA_kmeans: all scripts (Python + R) used for PCA-KMeans dimensional reduction and clustering analysis + scripts used for visualization of clusters </li>
    <li> Plots: scripts used in generation of main figures </li>
    <li> SCFR: scripts used in identification of SCFRs</li>
    <li> exon_shadow: files used in analysis of exon shadow </li>
  </ul>
</ul>

# Programming Languages Utilized
Major programming languages utilized: R, Shell scripts, Python

## R Packages

## Python Packages 

