import os
import GEOparse
import pandas as pd

df = pd.read_csv("extracted_dataset_ids.csv")

gse_ids = set()
for ids in df["GEO_IDs"].dropna():
    gse_ids.update(ids.split(", "))

print(f"Found {len(gse_ids)} unique GEO datasets.")

os.makedirs("GEO_Datasets", exist_ok=True)


def download_geo(gse_id):
    try:
        print(f"Downloading {gse_id} ...")
        gse = GEOparse.get_GEO(geo=gse_id, destdir="GEO_Datasets", silent=True)

        # Convert metadata dictionary to DataFrame
        meta_df = pd.DataFrame(list(gse.metadata.items()), columns=["Key", "Value"])
        meta_file = f"GEO_Datasets/{gse_id}_metadata.csv"
        meta_df.to_csv(meta_file, index=False)
        print(f"Saved metadata: {meta_file}")

        if gse.gpls:
            for gpl in gse.gpls:
                expr_file = f"GEO_Datasets/{gse_id}_{gpl}_expression.csv"
                gse.gpls[gpl].table.to_csv(expr_file, index=False)
                print(f"Saved expression data: {expr_file}")

    except Exception as e:
        print(f"Error downloading {gse_id}: {e}")


for gse_id in gse_ids:
    download_geo(gse_id)

print("All GEO datasets downloaded!")
