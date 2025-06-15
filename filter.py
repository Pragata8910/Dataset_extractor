import pandas as pd
import re

df = pd.read_csv("pubmed_cancer_omics.csv")

patterns = {
    "GEO": r"GSE\d{3,6}|GSM\d{3,6}",
    "SRA": r"SRR\d{3,6}|SRP\d{3,6}",
    "PRIDE": r"PXD\d{3,6}",
    "ENA": r"ERP\d{3,6}|ERX\d{3,6}",
}

combined_pattern = re.compile("|".join(patterns.values()), re.IGNORECASE)

df["Has_Dataset"] = (
    df["Abstract"].fillna("").apply(lambda x: bool(combined_pattern.search(x)))
)
filtered_df = df[df["Has_Dataset"]]

filtered_df.to_csv("filtered_articles_with_datasets.csv", index=False)
print(f"Found {len(filtered_df)} articles with dataset IDs.")
print("Saved as filtered_articles_with_datasets.csv")
