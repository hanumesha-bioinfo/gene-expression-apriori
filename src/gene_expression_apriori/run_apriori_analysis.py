from pathlib import Path

import pandas as pd
from mlxtend.frequent_patterns import apriori, association_rules

ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

binary = pd.read_csv(
    PROCESSED / "apriori_binary_matrix.csv",
    index_col="GSM_ID",
)
groups = pd.read_csv(
    PROCESSED / "apriori_sample_groups.csv",
    index_col="GSM_ID",
)

# Ensure sample order matches before analysis.
if set(binary.index) != set(groups.index):
    raise ValueError("Binary matrix and sample group IDs do not match.")

groups = groups.loc[binary.index]

all_itemsets = []
all_rules = []

for group_name in ["Tumor", "Normal"]:
    sample_ids = groups.index[groups["group"] == group_name]
    group_data = binary.loc[sample_ids].astype(bool)

    print(f"\n=== {group_name.upper()} ANALYSIS ===")
    print("Samples:", len(group_data))
    print("Genes:", group_data.shape[1])

    itemsets = apriori(
        group_data,
        min_support=0.20,
        use_colnames=True,
        max_len=3,
        low_memory=True,
    )

    if itemsets.empty:
        print("No frequent itemsets found at this support threshold.")
        continue

    itemsets["itemsets"] = itemsets["itemsets"].apply(
        lambda items: " | ".join(sorted(items))
    )
    itemsets["group"] = group_name
    all_itemsets.append(itemsets)

    # Convert itemset strings back to frozensets for rule generation.
    frequent = apriori(
        group_data,
        min_support=0.20,
        use_colnames=True,
        max_len=3,
        low_memory=True,
    )

    if len(frequent) >= 2:
        rules = association_rules(
            frequent,
            metric="confidence",
            min_threshold=0.60,
        )

        if not rules.empty:
            rules["antecedents"] = rules["antecedents"].apply(
                lambda items: " | ".join(sorted(items))
            )
            rules["consequents"] = rules["consequents"].apply(
                lambda items: " | ".join(sorted(items))
            )
            rules["group"] = group_name
            all_rules.append(rules)

    print("Frequent itemsets:", len(frequent))
    print(
        "Association rules:",
        len(rules) if len(frequent) >= 2 else 0,
    )

if all_itemsets:
    itemset_df = pd.concat(all_itemsets, ignore_index=True)
    itemset_df.to_csv(TABLES / "frequent_itemsets_by_group.csv", index=False)
    print("\nSaved:", TABLES / "frequent_itemsets_by_group.csv")
else:
    print("\nNo frequent itemsets were generated.")

if all_rules:
    rules_df = pd.concat(all_rules, ignore_index=True)
    rules_df = rules_df.sort_values(
        ["group", "lift", "confidence"],
        ascending=[True, False, False],
    )
    rules_df.to_csv(TABLES / "association_rules_by_group.csv", index=False)
    print("Saved:", TABLES / "association_rules_by_group.csv")

    print("\nTop rules by lift:")
    print(
        rules_df[
            ["group", "antecedents", "consequents",
             "support", "confidence", "lift"]
        ].head(10).to_string(index=False)
    )
else:
    print("No association rules met the selected thresholds.")
