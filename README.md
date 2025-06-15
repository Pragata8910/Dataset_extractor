# NCBI Dataset Pipeline

An automated pipeline for discovering, filtering, and extracting biomedical datasets from NCBI databases using keyword-based searches. This tool streamlines the process of dataset discovery by automating API interactions and generating structured CSV outputs.

##  Overview

This pipeline helps researchers quickly identify relevant datasets from NCBI's vast collection by:
- Searching datasets using custom keywords
- Filtering results based on specific criteria
- Extracting and classifying relevant dataset information
- Generating organized CSV files for further analysis

##  Project Structure

```
pubmed_env/
├── classify.py              # Dataset classification logic
├── extract.py              # Data extraction utilities
├── fetch.py                # NCBI API interaction (main)
├── fetch_2.py              # Alternative/enhanced fetching methods
├── filter.py               # Primary filtering logic
├── filter_2.py             # Additional filtering options
├── GEO.py                  # Gene Expression Omnibus specific handling
├── classified_articles_*.xlsx    # Classified article outputs
├── extracted_datas_*.xlsx        # Raw extracted data
├── filtered_articles_*.xlsx      # Filtered results
├── Final_classified_*.xlsx       # Final processed datasets
├── pubmed_cancer_*.xlsx          # Cancer-specific dataset results
└── USP14_proteomics_*.xlsx       # Proteomics-specific results
```

##  Features

- **Multi-source Support**: Works with PubMed, GEO, and other NCBI databases
- **Keyword-based Search**: Flexible search using custom keywords
- **Smart Filtering**: Multiple filtering strategies for refined results
- **Classification**: Automatic categorization of discovered datasets
- **Export Options**: Generates Excel/CSV files for easy data handling
- **Specialized Modules**: Domain-specific processing (cancer research, proteomics)

##  Requirements

```bash
pip install requests pandas openpyxl biopython
```

##  Usage

### Basic Dataset Fetching
```python
python fetch.py
```

### With Custom Keywords
```python
python classify.py
```

### Apply Filters
```python
python filter.py
```

### Extract Specific Data
```python
python extract.py
```

### GEO-specific Processing
```python
python GEO.py
```

##  Output Files

The pipeline generates several types of output files:

- **`extracted_datas_*.xlsx`**: Raw data pulled from NCBI APIs
- **`filtered_articles_*.xlsx`**: Results after applying filter criteria
- **`classified_articles_*.xlsx`**: Categorized datasets with metadata
- **`Final_classified_*.xlsx`**: Final processed results ready for analysis

##  Use Cases

- **Research Discovery**: Find datasets relevant to specific research topics
- **Meta-analysis Preparation**: Collect datasets for systematic reviews
- **Data Mining**: Extract structured information from biomedical literature
- **Workflow Automation**: Reduce manual effort in dataset curation

##  Configuration

Modify the keyword lists and filtering criteria in the respective Python files to customize searches for your specific research domain.

##  Contributing

Feel free to submit issues, fork the repository, and create pull requests for any improvements.


## 📧 Contact

Email: [pragata](pragata2004@gmail.com)

Phone No: [mobile](+91 9475170335)

**Note**: This tool is designed to work with publicly available NCBI data. Please ensure compliance with NCBI's usage policies and rate limits.
