# Gene Expression Association Rule Mining Using Apriori

An exploratory bioinformatics project that applies the Apriori algorithm to gene-expression data from the NCBI Gene Expression Omnibus (GEO) dataset **GSE10072**. The project processes gene-expression measurements, converts expression levels into binary transactions, and mines frequent gene combinations and association rules in tumor and normal sample groups.

## Project Overview

The workflow includes:

- Downloading and processing public GEO data.
- Mapping microarray probes to gene identifiers.
- Preparing gene-expression and sample metadata matrices.
- Selecting 200 highly variable genes for exploratory analysis.
- Converting gene expression into binary transactions using cohort-level median thresholds.
- Mining frequent itemsets and association rules separately in tumor and normal samples.
- Filtering rules by support, confidence, and lift.
- Generating visualizations and automated tests.

## Dataset

| Property | Description |
|---|---|
| GEO accession | GSE10072 |
| Data source | NCBI Gene Expression Omnibus |
| Platform | Affymetrix Human Genome U133A Array (GPL96) |
| Samples analyzed | 107 |
| Tumor samples | 58 |
| Normal samples | 49 |
| Genes selected for Apriori | 200 |

The processed expression matrix contains 12,547 mapped genes before selection of the 200 most variable genes.

Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10072

## Methodology

1. Obtain the GEO series matrix and platform annotation.
2. Prepare the expression matrix and map probes to genes.
3. Prepare sample labels and sample characteristics.
4. Select the 200 most variable genes.
5. Convert expression values to binary transactions using median expression thresholds.
6. Run Apriori separately for Tumor and Normal samples.
7. Generate association rules and filter them using selected thresholds.
8. Summarize the results, create figures, and run automated tests.

### Rule-filtering thresholds

- Minimum support: 0.25
- Minimum confidence: 0.70
- Minimum lift: 1.20

These thresholds are exploratory analysis settings, not validated biological significance cutoffs.

## Example Results

The analysis generates frequent itemsets and association rules for each sample group. Filtered rules are summarized by support, confidence, and lift.

Generated visualizations:

- `results/figures/filtered_rules_by_group.png`
- `results/figures/top_association_rules_by_lift.png`

Selected result tables include:

- `results/tables/selected_variable_genes.csv`
- `results/tables/top_rules_by_group.csv`
- `results/tables/rule_summary_by_group.csv`

The figures and rule tables are exploratory outputs. A high lift or confidence value does not, by itself, establish biological importance.

## Installation

Use Python 3.12 or another compatible Python 3 version.

```bash
git clone https://github.com/hanumesha-bioinfo/gene-expression-apriori.git
cd gene-expression-apriori

python3 -m venv venv
source venv/bin/activate

python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Reproducibility

Create and activate a Python environment, then install dependencies as described above.

### Run the newer GSE10072 workflow

From the repository root, execute these commands in order:

1. Download the GEO series record and expression matrix:

   `python download_geo_data.py`

2. Prepare the GPL96 probe-to-gene mapping:

   `python src/gene_expression_apriori/prepare_gene_mapping.py`

3. Build the gene-level expression matrix:

   `python src/gene_expression_apriori/prepare_expression_matrix.py`

4. Prepare sample metadata and tumor/normal labels:

   `python src/gene_expression_apriori/prepare_sample_metadata.py`

5. Select variable genes and create binary transactions:

   `python src/gene_expression_apriori/prepare_apriori_data.py`

6. Mine frequent itemsets and association rules:

   `python src/gene_expression_apriori/run_apriori_analysis.py`

7. Filter association rules:

   `python src/gene_expression_apriori/filter_association_rules.py`

8. Generate figures and compact summaries:

   `python src/gene_expression_apriori/visualize_apriori_results.py`

**Resource note:** Step 6 can generate millions of rules and very large intermediate CSV files. Ensure sufficient disk space and memory before running it. Step 7 reads the generated rules in chunks to reduce memory use.

The `inspect_sample_characteristics.py` and `check_patient_pairs.py` scripts are optional exploratory utilities. They require their corresponding generated metadata files.

The original `run_geo_apriori.py` is a separate, earlier workflow, not a replacement for the sequence above.

Raw GEO downloads, generated sample-level metadata, and large intermediate files are excluded from Git. The compact result tables and figures in this repository are outputs from a completed local run.

## Tests

Run the automated tests with:

```bash
python -m pytest -v
```

The current test suite checks the processed matrices, sample groups, filtered-rule thresholds, summary output, and figure generation.

## Limitations

- Association rules describe co-occurrence patterns; they do not establish causation or regulatory relationships between genes.
- The analysis uses a selected subset of variable genes and cohort-level median binarization, which can influence the resulting rules.
- Tumor and normal sample groups may differ in characteristics such as sex and patient composition.
- Some patients contribute both tumor and normal samples, but the exploratory Apriori analysis treats sample transactions independently.
- The number of rules depends strongly on the selected support, confidence, and lift thresholds.
- Findings require further biological validation and should not be interpreted as clinical conclusions.

## Technologies

Python, pandas, NumPy, GEOparse, mlxtend, Matplotlib, and pytest.

## Author

Hanumesha

GitHub: https://github.com/hanumesha-bioinfo
