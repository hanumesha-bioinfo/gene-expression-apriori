from io import StringIO
from pathlib import Path
import gzip

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
PROCESSED.mkdir(parents=True, exist_ok=True)

matrix_path = RAW / "GSE10072_series_matrix.txt.gz"
mapping_path = PROCESSED / "gpl96_clean_probe_gene_mapping.csv"

# Read only the expression table between GEO markers.
with gzip.open(matrix_path, "rt") as handle:
    table_lines = []
    inside_table = False

    for line in handle:
        line = line.rstrip("\n")

        if line.strip() == "!series_matrix_table_begin":
            inside_table = True
            continue

        if line.strip() == "!series_matrix_table_end":
            break

        if inside_table:
            table_lines.append(line)

if not table_lines:
    raise RuntimeError("Expression table markers were not found.")

probes = pd.read_csv(StringIO("\n".join(table_lines)), sep="\t")
probes = probes.rename(columns={probes.columns[0]: "ID"})

sample_columns = [c for c in probes.columns if c != "ID"]
probes[sample_columns] = probes[sample_columns].apply(
    pd.to_numeric, errors="coerce"
)

if probes["ID"].duplicated().any():
    raise ValueError("Duplicate probe IDs found in expression matrix.")

mapping = pd.read_csv(mapping_path)
mapping = mapping[["ID", "Gene Symbol"]].drop_duplicates()

# Keep only probes with an unambiguous gene-symbol mapping.
merged = probes.merge(mapping, on="ID", how="inner", validate="one_to_one")

if merged.empty:
    raise RuntimeError("No expression probes matched the annotation.")

# Multiple probes can map to one gene: summarize with the median.
gene_expression = (
    merged.drop(columns="ID")
    .groupby("Gene Symbol", sort=True)[sample_columns]
    .median()
)

gene_expression.index.name = "Gene Symbol"
output_path = PROCESSED / "gene_expression_matrix.csv"
gene_expression.to_csv(output_path)

print("=== GENE-LEVEL EXPRESSION MATRIX ===")
print("Matched probes:", merged["ID"].nunique())
print("Unique genes:", gene_expression.shape[0])
print("Samples:", gene_expression.shape[1])
print("Missing values:", int(gene_expression.isna().sum().sum()))
print("Expression range:",
      round(float(gene_expression.min().min()), 3),
      "to",
      round(float(gene_expression.max().max()), 3))
print("Saved:", output_path)
print("\nPreview:")
print(gene_expression.iloc[:5, :5].to_string())
