import os
import re
import pandas as pd
from collections import defaultdict

def gene_enrichment_processor(files, final_path):

    def clean_pathway(pathway, species_name="Nomascus leucogenys"):
        if pd.isna(pathway):
            return pathway

        pathway = str(pathway).strip()

        # remove KEGG ID prefix (Path:nleXXXXX )
        pathway = re.sub(r"^[^ ]+\s+", "", pathway)

        # only for gibbon: remove "-Nomascus leucogenys (...)"
        if species_name in pathway:
            pathway = re.sub(rf"\s*-\s*{species_name}.*$", "", pathway)

        return pathway.strip()

    def parse_genes(gene_string):
        if pd.isna(gene_string) or gene_string == "":
            return set()
        return set(re.split(r"[|\s]+", str(gene_string).strip()))

    # property → species → gene set
    property_species_genes = defaultdict(dict)

    # property → species → enrichment quantity
    property_species_metric = defaultdict(dict)

    for species, path in files.items():
        if not os.path.exists(path):
            continue

        df = pd.read_csv(path)

        if "Pathway" not in df.columns or "Genes" not in df.columns:
            raise KeyError(f"Missing required columns in {path}")

        # Check for required columns for metric
        if not {"nGenes", "Pathway Genes", "Fold Enrichment"}.issubset(df.columns):
            raise KeyError(f"{path} missing nGenes, Pathway Genes, or Fold Enrichment")

        df["Pathway"] = df["Pathway"].apply(clean_pathway, species_name=species)

        for _, row in df.iterrows():
            prop = row["Pathway"]
            genes = parse_genes(row["Genes"])
            if genes:
                property_species_genes[prop][species] = genes

            # --- compute enrichment_quantity ---
            if pd.notna(row["nGenes"]) and pd.notna(row["Pathway Genes"]) and pd.notna(row["Fold Enrichment"]):
                if row["Fold Enrichment"] > 0 and row["nGenes"] > 0 and row["Pathway Genes"] > 0:
                    metric = ((row["nGenes"] / row["Pathway Genes"]) * row["Fold Enrichment"]) ** 0.5
                else:
                    metric = 0
            else:
                metric = 0

            property_species_metric[prop][species] = metric

    rows = []

    for prop, species_genes in property_species_genes.items():
        species_list = sorted(species_genes.keys())
        gene_sets = list(species_genes.values())

        # common genes
        common_genes = set.intersection(*gene_sets) if gene_sets else set()

        # gene → species map
        gene_species_map = defaultdict(set)
        for sp, genes in species_genes.items():
            for g in genes:
                gene_species_map[g].add(sp)

        # enrichment quantity per species
        metric_dict = property_species_metric[prop]

        rows.append({
            "property": prop,
            "n_species": len(species_list),
            "species": "; ".join(species_list),

            "n_Genes_per_species": "; ".join(
                f"{sp}:{len(species_genes[sp])}"
                for sp in species_list
            ),

            "total_pathway_genes_per_species": "; ".join(
                f"{sp}:{'|'.join(sorted(species_genes[sp]))}"
                for sp in species_list
            ),

            "genes_common_in_all_species": "|".join(sorted(common_genes)),
            "genes_common_in_all_species_count": len(common_genes),

            "genes_unique_to_one_species": "; ".join(
                f"{sp}:{'|'.join(sorted(species_genes[sp] - common_genes))}"
                for sp in species_list
            ),

            # gene → number of species
            "gene_species_count": "|".join(
                f"{g}:{len(gene_species_map[g])}"
                for g in sorted(gene_species_map)
            ),

            # gene → species list
            "gene_species_list": "|".join(
                f"{g}:{','.join(sorted(gene_species_map[g]))}"
                for g in sorted(gene_species_map)
            ),

            # --- enrichment quantity ---
            "enrichment_quantity_per_species": "; ".join(
                f"{sp}:{metric_dict.get(sp, 0):.4f}" for sp in species_list
            )
        })

    pd.DataFrame(rows).to_csv(final_path, index=False)

## FOR GENES WITH NO SHADOW - ALL GENES AS BACKGROUND 
genes_no_shadow_allgenes_bg = {
    "bonobo": "data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/genes_without_shadow_nobg_enrichment.csv",
    "chimpanzee": "data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    "gorilla": "data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/genes_without_shadow_nobg_enrichment.csv",
    "borangutan": "data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    "human": "data/exon_shadow/gene_enrichment/human_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    # "gibbon": "PATH/TO/gibbon.csv",
    "sorangutan": "data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
}

genes_no_shadow_allgenes_bg_path = r"gene_enrichment_data_processing/results/combined_list_genes_no_shadow_allgenes_bg.csv"

## GENES WITH NO SHADOW - NO GENES AS BACKGROUND 
genes_no_shadow_nobg_files = {
    "bonobo": "data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/genes_without_shadow_nobg_enrichment.csv",
    "chimpanzee": "data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    "gorilla": "data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/genes_without_shadow_nobg_enrichment.csv",
    "borangutan": "data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    "human": "data/exon_shadow/gene_enrichment/human_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
    # "gibbon": "PATH/TO/gibbon.csv",
    "sorangutan": "data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/genes_without_shadow_no_bg_enrichment.csv",
}
genes_no_shadow_nobg_path = "gene_enrichment_data_processing/results/combined_list_genes_without_shadow_no_bg.csv"

## NO DOWNSTREAM SHADOW - NO BACKGROUND GENES
no_ds_shadow_no_bg_files = {'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/no_downstream_no_bg_enrichment.csv', 
                            'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/no_downstream_no_bg_enrichment.csv', 
                            'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/no_downstream_shadow_no_bg_enrichment.csv', 
                            'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/no_downstream_no_bg_enrichment.csv', 
                            'gorillla': 'data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/no_downstream_no_bg_enrichment.csv', 
                            'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/no_downstream_no_bg_enrichment.csv', 
                            'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/no_downstream_no_bg_enrichment.csv'}
no_ds_shadow_no_bg_path = "gene_enrichment_data_processing/results/no_ds_shadow_no_bg.csv"

## NO UPSTREAM SHADOW ALL GENES BACKGROUND ENRICHMENT 
no_us_shadow_allgenes_bg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'gorillla': 'data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/no_upstream_allgenes_bg_enrichment.csv'
    } 
no_us_shadow_allgenes_bg_path = "gene_enrichment_data_processing/results/no_us_shadow_allgenes_bg.csv"

## NO UPSTREAM SHADOW NO BACKGROUND 
no_us_shadow_nobg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'gorillla': 'data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/no_upstream_no_bg_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/no_upstream_no_bg_enrichment.csv'
}
no_us_shadow_nobg_path = "gene_enrichment_data_processing/results/no_us_shadow_no_bg.csv"

## NO UPSTREAM SHADOW - ALL GENES WITH SHADOW BACKGROUND
no_us_shadow_geneswshadow_bg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'gorilla': 'data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/no_upstream_geneswshadow_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/no_upstream_geneswshadow_enrichment.csv'
}
no_us_shadow_geneswshadow_bg_path = "gene_enrichment_data_processing/results/no_us_shadow_geneswshadow_bg.csv"

## SYMMETRIC SHADOW ALL GENES BG ENRICHMENT 
sym_shadow_allgenes_bg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/symmetric_allgenes_bg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/symmetric_allgenes_bg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/symmetric_allgenes_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/symmetric_allgenes_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/symmetric_allgenes_bg_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/symmetric_allgenes_bg_enrichment.csv'
}
sym_shadow_allgenes_path = "gene_enrichment_data_processing/results/sym_shadow_allgenes_bg.csv"

## SYMMETRIC SHADOWS NO BACKGROUND ENRICHMENT 
sym_shadow_nobg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/symmetric_no_bg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/symmetric_no_bg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/symmetric_no_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/symmetric_no_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/symmetric_no_bg_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/symmetric_no_bg_enrichment.csv',
    'gorilla': 'data\exon_shadow/gene_enrichment/gorillla_gene_enrichment/shadow_symmetric_no_bg_enrichment.csv'
}
sym_shadow_nobg_path = "gene_enrichment_data_processing/results/sym_shadow_nobg.csv"

## SYMMETRIC SHADOW ALL GENES WITH SHADOW ENRICHMENT 
sym_shadow_geneswshadow_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/symmetric_geneswshadow_bg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/symmetric_geneswshadow_bg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/symmetric_geneswshadow_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/symmetric_geneswshadow_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/symmetric_geneswshadow_bg_enrichment.csv'
}
## not there in gorilla and sumatran orangutan 
sym_shadow_geneswshadow_path = "gene_enrichment_data_processing/results/sym_shadow_geneswshadow_bg.csv"

## GENES WITH SHADOW/ NO BACKGROUND SPECIFIED
geneswshadow_no_bg_files = {
    'bonobo': 'data/exon_shadow/gene_enrichment/bonobo_gene_enrichment/geneswshadow_nobg_enrichment.csv', 
    'borangutan': 'data/exon_shadow/gene_enrichment/borangutan_gene_enrichment/geneswshadow_nobg_enrichment.csv', 
    'chimpanzee': 'data/exon_shadow/gene_enrichment/chimpanzee_gene_enrichment/geneswshadow_no_bg_enrichment.csv', 
    'gibbon': 'data/exon_shadow/gene_enrichment/gibbon_gene_enrichment/geneswshadow_no_bg_enrichment.csv', 
    'gorilla': 'data/exon_shadow/gene_enrichment/gorillla_gene_enrichment/geneswshadow_no_bg_enrichment.csv', 
    'human': 'data/exon_shadow/gene_enrichment/human_gene_enrichment/geneswshadow_nobg_enrichment.csv', 
    'sorangutan': 'data/exon_shadow/gene_enrichment/sorangutan_gene_enrichment/geneswshadow_no_bg_enrichment.csv'
}
geneswshadow_no_bg_path = "gene_enrichment_data_processing/results/geneswsahdow_no_bg.csv"

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
## FUNCTION CALLS

## genes with no shadow: all genes as background 
gene_enrichment_processor(files=genes_no_shadow_allgenes_bg, final_path=genes_no_shadow_allgenes_bg_path)

## genes with no shadow: no background (default background used)
gene_enrichment_processor(files=genes_no_shadow_nobg_files, final_path=genes_no_shadow_nobg_path)

## no downstream shadow: no background (default background)
gene_enrichment_processor(files=no_ds_shadow_no_bg_files, final_path=no_ds_shadow_no_bg_path)

## no upstream shadow: all genes as background 
gene_enrichment_processor(files=no_us_shadow_allgenes_bg_files, final_path=no_us_shadow_allgenes_bg_path)

## no upstream shadow: no background (default background)
gene_enrichment_processor(files=no_us_shadow_nobg_files, final_path=no_us_shadow_nobg_path)

## no upstream shadow: all genes with shadow as background 
gene_enrichment_processor(files=no_us_shadow_geneswshadow_bg_files, final_path=no_us_shadow_geneswshadow_bg_path)

## symmetric shadow: all genes as background 
gene_enrichment_processor(files=sym_shadow_allgenes_bg_files, final_path=sym_shadow_allgenes_path)

## symmetric shadow: no background 
gene_enrichment_processor(files=sym_shadow_nobg_files, final_path=sym_shadow_nobg_path)

## symmetric shadow: all genes with shadow as background
gene_enrichment_processor(files=sym_shadow_geneswshadow_files, final_path=sym_shadow_geneswshadow_path)

## genes with shadow enrichment with no background 
gene_enrichment_processor(files=geneswshadow_no_bg_files, final_path=geneswshadow_no_bg_path)