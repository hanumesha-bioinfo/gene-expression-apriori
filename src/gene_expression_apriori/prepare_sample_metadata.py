from pathlib import Path
import GEOparse
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
PROCESSED.mkdir(parents=True, exist_ok=True)

series_path = RAW / "GSE10072.txt"

# Load the locally downloaded GEO series metadata.
gse = GEOparse.get_GEO(filepath=str(series_path), silent=True)

records = []

for gsm_id, sample in gse.gsms.items():
    metadata = sample.metadata
    title = metadata.get("title", [""])[0]

    # Determine tissue group from the sample title.
    title_lower = title.lower()

    if "tumor" in title_lower:
        group = "Tumor"
    elif "normal" in title_lower:
        group = "Normal"
    else:
        group = "Unknown"

    records.append({
        "GSM_ID": gsm_id,
        "title": title,
        "group": group,
    })

metadata_df = pd.DataFrame(records)

expression_path = PROCESSED / "gene_expression_matrix.csv"
expression_df = pd.read_csv(expression_path, index_col=0)
expression_samples = set(expression_df.columns)

metadata_df = metadata_df[
    metadata_df["GSM_ID"].isin(expression_samples)
].copy()

missing_samples = expression_samples - set(metadata_df["GSM_ID"])
unknown_samples = metadata_df.loc[
    metadata_df["group"] == "Unknown", "GSM_ID"
].tolist()

if missing_samples:
    raise ValueError(
        f"Expression samples missing from metadata: {sorted(missing_samples)}"
    )

if unknown_samples:
    print("WARNING: Samples with unrecognized tissue labels:", unknown_samples)

output_path = PROCESSED / "sample_metadata.csv"
metadata_df.to_csv(output_path, index=False)

print("=== SAMPLE METADATA SUMMARY ===")
print("Samples matched:", len(metadata_df))
print("\nGroup counts:")
print(metadata_df["group"].value_counts().to_string())
print("\nSample preview:")
print(metadata_df.head(10).to_string(index=False))
print("\nSaved:", output_path)
