from pathlib import Path
import GEOparse
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"

gse = GEOparse.get_GEO(
    filepath=str(RAW / "GSE10072.txt"),
    silent=True,
)

records = []

for gsm_id, sample in gse.gsms.items():
    meta = sample.metadata
    row = {
        "GSM_ID": gsm_id,
        "title": meta.get("title", [""])[0],
    }

    # Preserve GEO's original characteristic labels and values.
    for characteristic in meta.get("characteristics_ch1", []):
        if ":" in characteristic:
            key, value = characteristic.split(":", 1)
            row[key.strip()] = value.strip()

    records.append(row)

characteristics = pd.DataFrame(records)

metadata_path = PROCESSED / "sample_metadata.csv"
groups = pd.read_csv(metadata_path)

merged = groups.merge(
    characteristics,
    on=["GSM_ID", "title"],
    how="left",
    validate="one_to_one",
)

output_path = PROCESSED / "sample_characteristics.csv"
merged.to_csv(output_path, index=False)

print("=== AVAILABLE GEO CHARACTERISTICS ===")
print("Samples:", len(merged))
print("Columns:")
for column in merged.columns:
    print(" -", column)

for column in merged.columns:
    if column.lower() in {"gender", "sex", "tissue", "disease state"}:
        print(f"\nDistribution for {column}:")
        print(pd.crosstab(merged["group"], merged[column], dropna=False))

print("\nSaved:", output_path)
print("\nPreview:")
print(merged.head(8).to_string(index=False))
