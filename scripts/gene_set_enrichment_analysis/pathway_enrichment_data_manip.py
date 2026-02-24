import pandas as pd
import itertools

# ================================
# 1. USER INPUT (NO AUTODETECTION)
# ================================

# Each tuple: (csv_path, foreground_label, background_label)
FILES = [
    ("gene_enrichment_data_processing/results/combined_list_genes_no_shadow_allgenes_bg.csv", "genes_without_shadow", "all_genes"),
    ("gene_enrichment_data_processing/results/combined_list_genes_without_shadow_no_bg.csv", "genes_without_shadow", "no_bg"),
    ("gene_enrichment_data_processing/results/geneswsahdow_no_bg.csv", "genes_with_shadow", "no_bg"),
    ("gene_enrichment_data_processing/results/no_ds_shadow_no_bg.csv", "no_downstream", "no_bg"),
    ("gene_enrichment_data_processing/results/no_us_shadow_allgenes_bg.csv", "no_upstream", "all_genes"),
    ("gene_enrichment_data_processing/results/no_us_shadow_geneswshadow_bg.csv", "no_upstream", "genes_with_shadow"),
    ("gene_enrichment_data_processing/results/no_us_shadow_no_bg.csv", "no_upstream", "no_bg"),
    ("gene_enrichment_data_processing/results/sym_shadow_allgenes_bg.csv", "symmetric_shadow", "all_genes"),
    ("gene_enrichment_data_processing/results/sym_shadow_geneswshadow_bg.csv", "symmetric_shadow", "genes_with_shadow"),
    ("gene_enrichment_data_processing/results/sym_shadow_nobg.csv", "symmetric_shadow", "no_bg")
    # add all your files explicitly
]

FOREGROUNDS = [
    "symmetric_shadow",
    "no_upstream",
    "no_downstream",
    "genes_with_shadow",
    "genes_without_shadow"
]

BACKGROUNDS = [
    "no_bg",
    "all_genes",
    "genes_with_shadow"
]

dfs = []

# ================================
# 2. PROCESS FILES
# ================================
for path, fg, bg in FILES:
    df = pd.read_csv(path)

    # Clean pathway names if gibbon info is present
    df["property"] = df.apply(
    lambda row: str(row["property"]).split("-", 1)[0].strip()
    if "gibbon" in str(row["species"]).lower()
    else str(row["property"]),
    axis=1
    )

    # ================================
    # STANDARDIZE SPECIFIC PATHWAY NAMES
    # ================================

    df["property"] = df["property"].replace({
        "Pathways of neurodegeneration-multiple diseases":
            "Pathways of Neurodegeneration - Multiple Diseases",
        "Chemical carcinogenesis-reactive oxygen species":
            "Chemical Carcinogens - ROS"
    })


    required = {"property", "species", "enrichment_quantity_per_species"}
    missing_cols = required - set(df.columns)
    if missing_cols:
        raise ValueError(f"{path} missing columns: {missing_cols}")

    df["foreground"] = fg
    df["background"] = bg

    # Fix "gorillla" typo
    df["species"] = df["species"].str.replace(r"\bgorillla\b", "gorilla", regex=True)

    # Explode rows for multiple species in the 'species' column
    if df["species"].str.contains(";").any():
        df["species"] = df["species"].str.split(";")
        df = df.explode("species")

    df["species"] = df["species"].str.strip()

    # ================================
    # 2a. EXTRACT ENRICHMENT QUANTITY FOR THE SPECIES
    # ================================
    def extract_enrichment(row):
        text = row["enrichment_quantity_per_species"]
        if pd.isna(text) or text.strip() == "":
            return 0.0

        # Split by semicolon
        entries = text.split(";")
        q_dict = {}
        for e in entries:
            if ":" in e:
                sp, val = e.split(":")
                sp = sp.strip().replace("gorillla", "gorilla")
                try:
                    val = float(val.strip())
                except:
                    val = 0.0
                q_dict[sp] = val

        # Return value for this row’s species
        return q_dict.get(row["species"], 0.0)

    df["enrichment_quantity"] = df.apply(extract_enrichment, axis=1)

    dfs.append(df[[
        "property", "species", "foreground", "background", "enrichment_quantity"
    ]])

# ================================
# 3. CONCAT ALL CSVs
# ================================
df_present = pd.concat(dfs, ignore_index=True)

# ================================
# 4. BUILD FULL PERMUTATION SPACE
# ================================
all_pathways = df_present["property"].unique()
all_species = df_present["species"].unique()

full_grid = pd.DataFrame(
    itertools.product(
        all_pathways,
        all_species,
        FOREGROUNDS,
        BACKGROUNDS
    ),
    columns=["property", "species", "foreground", "background"]
)

# ================================
# 5. MERGE + FILL MISSING
# ================================
df_full = full_grid.merge(
    df_present,
    on=["property", "species", "foreground", "background"],
    how="left"
)

df_full["enrichment_quantity"] = df_full["enrichment_quantity"].fillna(0.0)

# ================================
# 6. HELPER COLUMNS
# ================================
df_full["fg_bg_pair"] = df_full["foreground"] + " | " + df_full["background"]
df_full["enriched"] = (df_full["enrichment_quantity"] > 0).astype(int)
df_full["enrichment_state"] = df_full["enriched"].map({1: "enriched", 0: "not_enriched"})

# ================================
# 7. EXPORT
# ================================
df_full.to_csv("pathway_enrichment/enrichment_long_format_with_quantity.csv", index=False)

print("DONE")
print("Rows:", df_full.shape[0])
print("Unique pathways:", df_full['property'].nunique())
print("Unique species:", df_full['species'].nunique())