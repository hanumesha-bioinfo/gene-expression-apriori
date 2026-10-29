
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "processed"

metadata_path = DATA / "sample_metadata.csv"
characteristics_path = DATA / "sample_characteristics.csv"

metadata = pd.read_csv(metadata_path)
clinical = pd.read_csv(characteristics_path)

required = {"GSM_ID", "title", "group"}
missing = required - set(metadata.columns)
if missing:
    raise ValueError(f"Missing metadata columns: {sorted(missing)}")

if metadata["GSM_ID"].duplicated().any():
    raise ValueError("Duplicate GSM_ID values in sample metadata.")

metadata["patient_id"] = metadata["title"].str.extract(
    r"(GT\d+)", expand=False
)

if metadata["patient_id"].isna().any():
    bad = metadata.loc[
        metadata["patient_id"].isna(), ["GSM_ID", "title"]
    ]
    raise ValueError(f"Could not extract patient IDs:\n{bad}")

valid_groups = {"Tumor", "Normal"}
unexpected = set(metadata["group"].dropna()) - valid_groups
if unexpected:
    raise ValueError(f"Unexpected sample groups: {sorted(unexpected)}")

patient_groups = metadata.groupby("patient_id")["group"].agg(
    lambda values: set(values)
)

paired_ids = patient_groups[
    patient_groups.apply(lambda groups: valid_groups.issubset(groups))
].index

tumor_only = patient_groups[
    patient_groups.apply(lambda groups: groups == {"Tumor"})
].index

normal_only = patient_groups[
    patient_groups.apply(lambda groups: groups == {"Normal"})
].index

print("=== SAMPLE STRUCTURE VALIDATION ===")
print(f"Total samples: {len(metadata)}")
print(f"Unique patients: {metadata['patient_id'].nunique()}")
print(f"Tumor samples: {(metadata['group'] == 'Tumor').sum()}")
print(f"Normal samples: {(metadata['group'] == 'Normal').sum()}")
print(f"Patients with both groups: {len(paired_ids)}")
print(f"Tumor-only patients: {len(tumor_only)}")
print(f"Normal-only patients: {len(normal_only)}")

print("\n=== AGE SUMMARY BY SAMPLE GROUP ===")
if "Age at Diagnosis" in clinical.columns:
    age = pd.to_numeric(clinical["Age at Diagnosis"], errors="coerce")
    age_summary = clinical.assign(age_numeric=age).groupby("group")[
        "age_numeric"
    ].agg(["count", "mean", "std", "median", "min", "max"])
    print(age_summary.round(2).to_string())
else:
    print("Age at Diagnosis column not found.")

print("\n=== CLINICAL MISSINGNESS ===")
clinical_columns = [
    "Gender",
    "Age at Diagnosis",
    "Cigarette Smoking Status",
    "Stage",
]
for column in clinical_columns:
    if column in clinical.columns:
        print(f"{column}: {clinical[column].isna().sum()} missing values")

print("\n=== GROUP DISTRIBUTIONS ===")
for column in ["Gender", "Cigarette Smoking Status", "Stage"]:
    if column in clinical.columns:
        print(f"\n{column}:")
        print(
            pd.crosstab(
                clinical["group"], clinical[column], dropna=False
            ).to_string()
        )

print("\nValidation completed. No source files were modified.")
