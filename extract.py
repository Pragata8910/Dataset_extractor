import pandas as pd
import re

df = pd.read_csv("filtered_articles_with_datasets.csv")

patterns = {
    "GEO": r"(GSE\d{3,6}|GSM\d{3,6})",
    "SRA": r"(SRR\d{3,6}|SRP\d{3,6})",
    "PRIDE": r"(PXD\d{3,6})",
    "ENA": r"(ERP\d{3,6}|ERX\d{3,6})",
}


def extract_ids(text, pattern):
    matches = re.findall(pattern, text, re.IGNORECASE)
    return ", ".join(set(matches)) if matches else None


for dataset, pattern in patterns.items():
    df[dataset + "_IDs"] = (
        df["Abstract"].fillna("").apply(lambda x: extract_ids(x, pattern))
    )

df_filtered = df.dropna(subset=[col + "_IDs" for col in patterns.keys()], how="all")

df_filtered.to_csv("extracted_dataset_ids.csv", index=False)
print("Extracted dataset IDs and saved to extracted_dataset_ids.csv")
