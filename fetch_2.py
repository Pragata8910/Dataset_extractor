from Bio import Entrez
import pandas as pd
from time import sleep
from Bio import Medline

Entrez.email = "pragata2004@gmail.com"

# 1. Cast the widest possible PubMed net
query = """
(USP14 OR "Ubiquitin-specific protease 14") 
AND 
(knockdown OR knockout OR KO OR "knock-out" OR silencing OR inhibition OR depletion)
AND
(proteom* OR transcriptom* OR "RNA-seq" OR "RNA sequencing" OR "mass spec*" OR microarray OR "gene expression" OR omics)
"""

print("Searching PubMed...")
handle = Entrez.esearch(db="pubmed", term=query, retmax=800, usehistory="y")
record = Entrez.read(handle)
pmids = record["IdList"]

print(f"Found {len(pmids)} articles.")


def fetch_details(id_list):
    ids = ",".join(id_list)
    fetch_handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
    records = list(Medline.parse(fetch_handle))
    fetch_handle.close()
    return records


all_records = []
for i in range(0, len(pmids), 200):
    batch = pmids[i : i + 200]
    print(f"Fetching articles {i + 1} to {i + len(batch)}...")
    records = fetch_details(batch)
    all_records.extend(records)
    sleep(0.5)


data = []
for r in all_records:
    data.append(
        {
            "PMID": r.get("PMID", ""),
            "Title": r.get("TI", ""),
            "Abstract": r.get("AB", ""),
            "Journal": r.get("JT", ""),
            "Publication Date": r.get("DP", ""),
            "DOI": r.get("AID", [""])[0].replace("[doi]", "").strip()
            if "AID" in r
            else "",
            "Databases": ", ".join(
                [x for x in r.get("OT", []) if "GSE" in x or "PRD" in x or "PXD" in x]
            ),
        }
    )

df = pd.DataFrame(data)
df.to_csv("USP14_proteomics_transcriptomics.csv", index=False)
print("Results saved to 'USP14_proteomics_transcriptomics.csv'.")
