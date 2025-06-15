import pandas as pd
import re

# Load data
df = pd.read_csv(
    "/Users/pragata/Codes/Projects/Oncolab_proj_1/USP14_proteomics_transcriptomics.csv"
)

# Expanded regex patterns for common omics databases
patterns = {
    "GEO": r"GSE\d{3,8}|GSM\d{3,8}|GDS\d{3,8}",  # GEO datasets/samples
    "SRA": r"SR[RPX]\d{3,8}|SAM[NED]\d{3,8}",  # SRA/NCBI BioProject
    "PRIDE": r"PXD\d{3,8}",  # Proteomics data
    "ProteomeXchange": r"PXD\d{3,8}|PMD\d{3,8}",
    "ENA": r"ERP\d{3,8}|ERX\d{3,8}",  # European Nucleotide Archive
    "ArrayExpress": r"E-[A-Z]{4}-\d{3,8}",  # Microarray data
    "DOI": r"10\.\d{4,9}/[-._;()/:A-Z0-9]+",  # Supplemental DOI links
}

# Combine patterns and compile for efficiency
combined_pattern = re.compile("|".join(patterns.values()), re.IGNORECASE)

# Check both 'Abstract' and 'Title' for dataset IDs
df["Has_Dataset"] = (
    df["Abstract"]
    .fillna("")
    .str.cat(df["Title"].fillna(""), sep=" ")
    .apply(lambda x: bool(combined_pattern.search(x)))
)

# Filter and save
filtered_df = df[df["Has_Dataset"]]
filtered_df.to_csv("USP14_filtered_articles.csv", index=False)

# Print summary
print(f"Found {len(filtered_df)} articles with dataset IDs out of {len(df)} total.")
print("Saved as 'USP14_filtered_articles.csv'")


# Optional: Extract the specific dataset IDs found
def extract_dataset_ids(text):
    return ", ".join(set(combined_pattern.findall(text)))


filtered_df["Dataset_IDs"] = (
    filtered_df["Abstract"]
    .fillna("")
    .str.cat(filtered_df["Title"].fillna(""), sep=" ")
    .apply(extract_dataset_ids)
)

# Save enriched data
filtered_df.to_csv("USP14_filtered_with_dataset_IDs.csv", index=False)
print("Also saved detailed version with extracted IDs.")
