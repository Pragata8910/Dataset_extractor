from Bio import Entrez
import pandas as pd
from time import sleep
from Bio import Medline


Entrez.email = "pragata2004@gmail.com"


query = '((prostate cancer) OR (cervical cancer)) AND (genomics OR proteomics OR "RNA sequencing")'


print("Searching PubMed...")
handle = Entrez.esearch(db="pubmed", term=query, retmax=2000)
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
    sleep(1)

data = []
for r in all_records:
    data.append(
        {
            "PMID": r.get("PMID", ""),
            "Title": r.get("TI", ""),
            "Abstract": r.get("AB", ""),
        }
    )

df = pd.DataFrame(data)
df.to_csv("pubmed_cancer_omics.csv", index=False)
print("Saved results to pubmed_cancer_omics.csv")
