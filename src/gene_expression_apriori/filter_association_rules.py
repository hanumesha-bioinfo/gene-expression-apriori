from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"

SOURCE = TABLES / "association_rules_by_group.csv"
FILTERED = TABLES / "association_rules_filtered.csv"
TOP_RULES = TABLES / "top_rules_by_group.csv"

if not SOURCE.exists():
    raise FileNotFoundError(f"Missing input file: {SOURCE}")

# Exploratory filtering thresholds.
MIN_SUPPORT = 0.25
MIN_CONFIDENCE = 0.70
MIN_LIFT = 1.20

total_rows = 0
kept_rows = 0
first_chunk = True
top_chunks = []

for chunk in pd.read_csv(SOURCE, chunksize=100_000):
    total_rows += len(chunk)

    metric_columns = ["support", "confidence", "lift"]
    for column in metric_columns:
        chunk[column] = pd.to_numeric(chunk[column], errors="coerce")

    valid = np.isfinite(chunk[metric_columns]).all(axis=1)

    filtered = chunk.loc[
        valid
        & (chunk["support"] >= MIN_SUPPORT)
        & (chunk["confidence"] >= MIN_CONFIDENCE)
        & (chunk["lift"] >= MIN_LIFT)
    ].copy()

    if not filtered.empty:
        filtered.to_csv(
            FILTERED,
            mode="w" if first_chunk else "a",
            header=first_chunk,
            index=False,
        )
        first_chunk = False
        kept_rows += len(filtered)

        # Retain a small ranked sample from each chunk.
        top_chunks.append(
            filtered.sort_values(
                ["lift", "confidence", "support"],
                ascending=False,
            ).groupby("group", group_keys=False).head(100)
        )

if first_chunk:
    # Create an empty output with the source columns if no rules qualify.
    columns = pd.read_csv(SOURCE, nrows=0).columns
    pd.DataFrame(columns=columns).to_csv(FILTERED, index=False)

if top_chunks:
    candidates = pd.concat(top_chunks, ignore_index=True)
    candidates = candidates.drop_duplicates(
        subset=["group", "antecedents", "consequents"]
    )

    top = (
        candidates.sort_values(
            ["lift", "confidence", "support"],
            ascending=False,
        )
        .groupby("group", group_keys=False)
        .head(20)
    )
    top.to_csv(TOP_RULES, index=False)
else:
    top = pd.DataFrame()

print("=== RULE FILTERING SUMMARY ===")
print("Original rules processed:", f"{total_rows:,}")
print("Rules passing thresholds:", f"{kept_rows:,}")
print("Thresholds:")
print(f"  Support >= {MIN_SUPPORT}")
print(f"  Confidence >= {MIN_CONFIDENCE}")
print(f"  Lift >= {MIN_LIFT}")
print("\nSaved:", FILTERED)
print("Saved:", TOP_RULES)

if not top.empty:
    print("\nTop filtered rules by group:")
    print(
        top[
            ["group", "antecedents", "consequents",
             "support", "confidence", "lift"]
        ].to_string(index=False)
    )
else:
    print("\nNo rules passed all three thresholds.")
