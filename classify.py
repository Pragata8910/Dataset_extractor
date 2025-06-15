import pandas as pd
import re


def classify_articles(input_file, output_file):
    """
    Classify research articles based on cancer type, omics type, and sample type.

    Parameters:
    - input_file: Path to the CSV file containing article information
    - output_file: Path to save the output CSV with classifications
    """
    # Load the dataset
    df = pd.read_csv(input_file)

    # Initialize classification columns
    df["Cancer_Type"] = "Unknown"
    df["Omics_Type"] = "Unknown"
    df["Sample_Type"] = "Unknown"

    # Define patterns for classification
    cancer_patterns = {
        "Prostate": [
            r"\bprostate\s+cancer\b",
            r"\bPCa\b",
            r"\bprostate\s+carcinoma\b",
            r"\bPRAD\b",
        ],
        "Cervical": [
            r"\bcervical\s+cancer\b",
            r"\bCC\b",
            r"\bcervical\s+carcinoma\b",
            r"\bCSCC\b",
        ],
    }

    omics_patterns = {
        "Transcriptomics": [
            r"\btranscriptom",
            r"\bRNA-seq",
            r"\bexpression\s+profil",
            r"\bmRNA",
            r"\bgene\s+expression",
            r"\bscRNA-seq",
        ],
        "Proteomics": [
            r"\bproteom",
            r"\bprotein\s+profil",
            r"\bprotein\s+express",
            r"\bMass\s+spectrometry",
            r"\bLC-MS/MS",
            r"\bproteins",
        ],
    }

    sample_patterns = {
        "Blood": [
            r"\bblood\b",
            r"\bserum\b",
            r"\bplasma\b",
            r"\bperipheral\s+blood\b",
            r"\bcirculating\b",
            r"\bPBMC\b",
        ],
        "Tissue": [
            r"\btissue\b",
            r"\btumor\s+tissue\b",
            r"\btumour\s+tissue\b",
            r"\bcancer\s+tissue\b",
            r"\bbiopsy\b",
        ],
    }

    # Process each article
    for index, row in df.iterrows():
        abstract = str(row["Abstract"]).lower() if pd.notna(row["Abstract"]) else ""
        title = str(row["Title"]).lower() if pd.notna(row["Title"]) else ""

        # Combined text for analysis
        text = abstract + " " + title

        # Classify cancer type
        for cancer_type, patterns in cancer_patterns.items():
            if any(re.search(pattern, text, re.IGNORECASE) for pattern in patterns):
                df.at[index, "Cancer_Type"] = cancer_type
                break

        # Classify omics type
        for omics_type, patterns in omics_patterns.items():
            if any(re.search(pattern, text, re.IGNORECASE) for pattern in patterns):
                df.at[index, "Omics_Type"] = omics_type
                break

        # Classify sample type
        for sample_type, patterns in sample_patterns.items():
            if any(re.search(pattern, text, re.IGNORECASE) for pattern in patterns):
                df.at[index, "Sample_Type"] = sample_type
                break

    # Save the classified data
    df.to_csv(output_file, index=False)
    print(f"Classification complete. Results saved to {output_file}")

    # Print summary
    print("\nClassification Summary:")
    print("----------------------")
    print(f"Cancer Type Distribution:\n{df['Cancer_Type'].value_counts()}")
    print(f"\nOmics Type Distribution:\n{df['Omics_Type'].value_counts()}")
    print(f"\nSample Type Distribution:\n{df['Sample_Type'].value_counts()}")

    # Print articles where classification was not possible
    unknown_articles = df[
        (df["Cancer_Type"] == "Unknown")
        | (df["Omics_Type"] == "Unknown")
        | (df["Sample_Type"] == "Unknown")
    ]

    if not unknown_articles.empty:
        print("\nArticles with incomplete classification:")
        for _, row in unknown_articles.iterrows():
            print(f"PMID: {row['PMID']}, Title: {row['Title']}")
            print(
                f"  Cancer: {row['Cancer_Type']}, Omics: {row['Omics_Type']}, Sample: {row['Sample_Type']}\n"
            )

    return df


# Example usage
if __name__ == "__main__":
    input_file = "extracted_dataset_ids.csv"
    output_file = "classified_articles.csv"

    classified_df = classify_articles(input_file, output_file)
