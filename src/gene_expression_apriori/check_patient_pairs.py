from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "processed"

metadata = pd.read_csv(DATA / "sample_metadata.csv")
characteristics = pd.read_csv(DATA / "sample_characteristics.csv")

# Extract the patient identifier from titles such as:
# "Lung Tumor_GT00006" and "Normal Lung_GT00006".
metadata["patient_id"] = metadata["title"].str.extract(
    r"(GT\d+)", expand=False
)

if metadata["patient_id"].isna().any():
    print("WARNING: Some samples have no patient ID.")
    print(metadata.loc[metadata["patient_id"].isna(), ["GSM_ID", "title"]])

patient_groups = (
    metadata.groupby("patient_id")["group"]
    .agg(lambda values: sorted(set(values)))
)

paired_patients = patient_groups[
    patient_groups.apply(lambda groups: {"Tumor", "Normal"}.issubset(groups))
].index

print("=== PATIENT PAIRING ===")
print("Total samples:", len(metadata))
print("Unique patient IDs:", metadata["patient_id"].nunique())
print("Patients with both tumor and normal samples:", len(paired_patients))
print("Patients with tumor samples only:",
      int(patient_groups.apply(lambda x: x == ["Tumor"]).sum()))
print("Patients with normal samples only:",
      int(patient_groups.apply(lambda x: x == ["Normal"]).sum()))

print("\n=== SAMPLES PER PATIENT ===")
print(metadata.groupby(["patient_id", "group"]).size().unstack(fill_value=0).head(12))

print("\n=== CLINICAL CHARACTERISTICS BY GROUP ===")
for column in ["Gender", "Cigarette Smoking Status", "Stage"]:
    if column in characteristics.columns:
        print(f"\n{column}:")
        print(pd.crosstab(
            characteristics["group"],
            characteristics[column],
            dropna=False
        ))

print("\n=== PAIRED PATIENT EXAMPLES ===")
print(
    metadata[metadata["patient_id"].isin(paired_patients)]
    .sort_values(["patient_id", "group"])
    [["patient_id", "GSM_ID", "title", "group"]]
    .head(20)
    .to_string(index=False)
)
