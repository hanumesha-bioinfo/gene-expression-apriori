from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"

RULES_FILE = TABLES / "association_rules_filtered.csv"
TOP_RULES_FILE = TABLES / "top_rules_by_group.csv"

FIGURES.mkdir(parents=True, exist_ok=True)

rules = pd.read_csv(RULES_FILE)
required = {"group", "support", "confidence", "lift"}

missing = required - set(rules.columns)
if missing:
    raise ValueError(
        f"Missing required columns in {RULES_FILE.name}: {sorted(missing)}"
    )

rules["support"] = pd.to_numeric(rules["support"], errors="coerce")
rules["confidence"] = pd.to_numeric(rules["confidence"], errors="coerce")
rules["lift"] = pd.to_numeric(rules["lift"], errors="coerce")

rules = rules.dropna(subset=["group", "support", "confidence", "lift"])

# 1. Compare the number of filtered rules per group.
counts = rules.groupby("group").size().sort_values(ascending=False)

fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(counts.index.astype(str), counts.values)

ax.set_title("Filtered Association Rules by Sample Group")
ax.set_xlabel("Sample group")
ax.set_ylabel("Number of filtered rules")
ax.bar_label(bars, fmt="%.0f", padding=3)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()

path1 = FIGURES / "filtered_rules_by_group.png"
fig.savefig(path1, dpi=300, bbox_inches="tight")
plt.close(fig)

# 2. Summarize key metrics for each group.
summary = (
    rules.groupby("group")
    .agg(
        filtered_rule_count=("lift", "size"),
        mean_support=("support", "mean"),
        mean_confidence=("confidence", "mean"),
        median_lift=("lift", "median"),
        maximum_lift=("lift", "max"),
    )
    .reset_index()
)

summary.to_csv(TABLES / "rule_summary_by_group.csv", index=False)

# 3. Plot the top five rules by lift per group.
# Use the filtered results so the chart is reproducible from one source.
top = (
    rules.sort_values(["group", "lift"], ascending=[True, False])
    .groupby("group", group_keys=False)
    .head(5)
    .copy()
)

def rule_label(row):
    antecedents = str(row.get("antecedents", ""))
    consequents = str(row.get("consequents", ""))
    if antecedents != "nan" and consequents != "nan":
        return f"{antecedents} → {consequents}"
    return f"Rule {row.name}"

top["rule_label"] = top.apply(rule_label, axis=1)

groups = list(top["group"].drop_duplicates())
fig, axes = plt.subplots(
    nrows=max(len(groups), 1),
    ncols=1,
    figsize=(12, max(5, 3.5 * len(groups))),
    squeeze=False,
)

for ax, group in zip(axes.flat, groups):
    subset = (
        top[top["group"] == group]
        .sort_values("lift", ascending=True)
    )
    bars = ax.barh(subset["rule_label"], subset["lift"])
    ax.set_title(f"Top Association Rules by Lift — {group}")
    ax.set_xlabel("Lift")
    ax.set_ylabel("Association rule")
    ax.bar_label(bars, fmt="%.2f", padding=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=8)

fig.tight_layout()

path2 = FIGURES / "top_association_rules_by_lift.png"
fig.savefig(path2, dpi=300, bbox_inches="tight")
plt.close(fig)

print("=== VISUALIZATION COMPLETE ===")
print("Filtered rules used:", len(rules))
print("\nRule summary:")
print(summary.to_string(index=False))
print("\nSaved figures:")
print(path1.relative_to(ROOT))
print(path2.relative_to(ROOT))
print("\nSaved summary:")
print((TABLES / "rule_summary_by_group.csv").relative_to(ROOT))
