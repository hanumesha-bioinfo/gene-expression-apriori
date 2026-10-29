from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results" / "tables"

EXPRESSION_PATH = PROCESSED / "gene_expression_matrix.csv"
METADATA_PATH = PROCESSED / "sample_metadata.csv"

expression = pd.read_csv(EXPRESSION_PATH, index_col=0)
metadata = pd.read_csv(METADATA_PATH)

# Validate sample alignment before analysis.
if expression.columns.duplicated().any():
    raise ValueError("Duplicate sample IDs in expression matrix.")

if metadata["GSM_ID"].duplicated().any():
    raise ValueError("Duplicate sample IDs in metadata.")

if set(expression.columns) != set(metadata["GSM_ID"]):
    raise ValueError("Expression samples and metadata samples do not match.")

metadata = metadata.set_index("GSM_ID").loc[expression.columns]

if not metadata["group"].isin(["Tumor", "Normal"]).all():
    raise ValueError("Unexpected or unknown sample groups found.")

# Select the 200 genes with the greatest variance across all samples.
gene_variances = expression.var(axis=1)
selected_genes = gene_variances.nlargest(200).index
selected_expression = expression.loc[selected_genes]

# Binarize each gene relative to its median across the full cohort.
# True means expression is above that gene's cohort median.
thresholds = selected_expression.median(axis=1)
binary = selected_expression.gt(thresholds, axis=0)

# Save a sample-by-gene binary matrix for reproducibility.
binary_samples = binary.T.astype(bool)
binary_samples.index.name = "GSM_ID"
binary_samples.to_csv(PROCESSED / "apriori_binary_matrix.csv")

# Save selected genes and thresholds for transparency.
gene_table = pd.DataFrame({
    "Gene Symbol": selected_expression.index,
    "Variance": gene_variances.loc[selected_expression.index].values,
    "Median Threshold": thresholds.loc[selected_expression.index].values,
})
gene_table.to_csv(RESULTS / "selected_variable_genes.csv", index=False)

# Save sample groups in the same order as the binary matrix.
sample_groups = metadata[["title", "group"]].copy()
sample_groups.index.name = "GSM_ID"
sample_groups.to_csv(PROCESSED / "apriori_sample_groups.csv")

print("=== APRIORI INPUT PREPARATION ===")
print("Samples:", binary_samples.shape[0])
print("Selected genes:", binary_samples.shape[1])
print("Tumor samples:", int((metadata["group"] == "Tumor").sum()))
print("Normal samples:", int((metadata["group"] == "Normal").sum()))
print("Missing binary values:", int(binary_samples.isna().sum().sum()))
print("Binary values (True = above cohort median)")
print(binary_samples.iloc[:5, :8].to_string())
print("\nSaved:")
print(" - data/processed/apriori_binary_matrix.csv")
print(" - data/processed/apriori_sample_groups.csv")
print(" - results/tables/selected_variable_genes.csv")
