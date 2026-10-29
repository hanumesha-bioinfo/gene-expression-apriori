from pathlib import Path
import csv
import gzip
import io
import re
import sys

import numpy as np
import pandas as pd
from mlxtend.frequent_patterns import apriori, association_rules
from mlxtend.preprocessing import TransactionEncoder


# ------------------------------------------------------------
# 1. Project paths and analysis settings
# ------------------------------------------------------------

PROJECT_DIR = Path(__file__).resolve().parent
RAW_DIR = PROJECT_DIR / "data" / "raw"
RESULTS_DIR = PROJECT_DIR / "results"

MATRIX_FILE = RAW_DIR / "GSE10072_series_matrix.txt.gz"

TOP_VARIABLE_PROBES = 20
MIN_SUPPORT = 0.30
MIN_CONFIDENCE = 0.70


# ------------------------------------------------------------
# 2. Read sample metadata and the expression table
# ------------------------------------------------------------

def read_geo_series_matrix(path):
    """Read GEO sample titles and the expression table from a
    compressed GEO Series Matrix file."""

    if not path.exists():
        raise FileNotFoundError(
            f"Expression matrix not found: {path}\n"
            "Download GSE10072_series_matrix.txt.gz first."
        )

    metadata = {}
    table_lines = []
    inside_table = False

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")

            if line == "!series_matrix_table_begin":
                inside_table = True
                continue

            if line == "!series_matrix_table_end":
                inside_table = False
                break

            if inside_table:
                table_lines.append(line)
                continue

            if line.startswith("!Sample_"):
                fields = next(csv.reader([line], delimiter="\t", quotechar='"'))
                if len(fields) > 1:
                    metadata[fields[0]] = fields[1:]

    if not table_lines:
        raise ValueError(
            "Could not find the expression table in the series matrix."
        )

    expression = pd.read_csv(
        io.StringIO("\n".join(table_lines)),
        sep="\t",
        index_col=0,
        quotechar='"',
        low_memory=False,
    )

    expression.index = expression.index.astype(str)
    expression.columns = expression.columns.astype(str)

    # Ensure every expression measurement is numeric.
    expression = expression.apply(pd.to_numeric, errors="coerce")
    expression = expression.replace([np.inf, -np.inf], np.nan)
    expression = expression.dropna(axis=0, how="all")

    # Drop probes that contain missing values in any sample.
    expression = expression.dropna(axis=0, how="any")

    if expression.empty:
        raise ValueError("No usable numeric expression values were found.")

    return expression, metadata


# ------------------------------------------------------------
# 3. Identify tumor and normal samples using sample metadata
# ------------------------------------------------------------

def classify_sample(metadata_text):
    """Classify a sample from its title and available metadata.
    Normal/non-tumor terms are checked before tumor terms."""

    text = str(metadata_text).lower()
    text = re.sub(r"[_/]+", " ", text)

    normal_pattern = (
        r"\bnormal\b|"
        r"\bnon[\s-]?tumou?r\b|"
        r"\bhealthy\b|"
        r"\bcontrol\b|"
        r"\badjacent normal\b"
    )

    tumor_pattern = (
        r"\btumou?r\b|"
        r"\bcancer\b|"
        r"\bcarcinoma\b|"
        r"\bmalignan(?:t|cy)\b"
    )

    if re.search(normal_pattern, text):
        return "Normal"

    if re.search(tumor_pattern, text):
        return "Tumor"

    return "Unknown"


def build_sample_labels(metadata, sample_ids):
    """Use GEO sample titles and other sample annotations."""

    titles = metadata.get("!Sample_title", [])
    label_by_sample = {}

    for position, sample_id in enumerate(sample_ids):
        pieces = []

        if position < len(titles):
            pieces.append(str(titles[position]))

        # Include other sample metadata as fallback evidence.
        for key, values in metadata.items():
            if key == "!Sample_title":
                continue
            if position < len(values):
                pieces.append(str(values[position]))

        label_by_sample[sample_id] = classify_sample(" ".join(pieces))

    labels = pd.Series(label_by_sample, name="group")
    labels.index.name = "sample_id"

    return labels


# ------------------------------------------------------------
# 4. Run association rule mining
# ------------------------------------------------------------

def mine_rules(transactions, min_support, min_confidence):
    """Encode transactions and calculate frequent itemsets and rules."""

    encoder = TransactionEncoder()
    encoded_array = encoder.fit(transactions).transform(transactions)

    encoded = pd.DataFrame(
        encoded_array,
        columns=encoder.columns_,
        dtype=bool,
    )

    frequent = apriori(
        encoded,
        min_support=min_support,
        use_colnames=True,
        max_len=2,
    )

    if frequent.empty:
        return encoded, frequent, pd.DataFrame()

    # Rules require itemsets containing at least two items.
    if not (frequent["itemsets"].map(len) >= 2).any():
        return encoded, frequent, pd.DataFrame()

    rules = association_rules(
        frequent,
        metric="confidence",
        min_threshold=min_confidence,
    )

    return encoded, frequent, rules


def save_mining_results(prefix, transactions, labels, probe_ids):
    """Mine rules and save readable CSV outputs."""

    if len(transactions) < 2:
        print(f"Skipping {prefix}: not enough samples.")
        return

    encoded, frequent, rules = mine_rules(
        transactions,
        MIN_SUPPORT,
        MIN_CONFIDENCE,
    )

    if not frequent.empty:
        frequent = frequent.copy()
        frequent["itemsets"] = frequent["itemsets"].apply(
            lambda items: ", ".join(sorted(items))
        )
        frequent.to_csv(
            RESULTS_DIR / f"{prefix}_frequent_itemsets.csv",
            index=False,
        )
    else:
        pd.DataFrame(columns=["support", "itemsets"]).to_csv(
            RESULTS_DIR / f"{prefix}_frequent_itemsets.csv",
            index=False,
        )

    if not rules.empty:
        rules = rules.copy()
        rules["antecedents"] = rules["antecedents"].apply(
            lambda items: ", ".join(sorted(items))
        )
        rules["consequents"] = rules["consequents"].apply(
            lambda items: ", ".join(sorted(items))
        )
        rules = rules.sort_values(
            ["confidence", "lift", "support"],
            ascending=False,
        )
        rules.to_csv(
            RESULTS_DIR / f"{prefix}_association_rules.csv",
            index=False,
        )
    else:
        pd.DataFrame(
            columns=[
                "antecedents", "consequents", "support",
                "confidence", "lift",
            ]
        ).to_csv(
            RESULTS_DIR / f"{prefix}_association_rules.csv",
            index=False,
        )

    print(
        f"{prefix}: {len(transactions)} samples, "
        f"{len(probe_ids)} probes, "
        f"{len(frequent)} frequent itemsets, "
        f"{len(rules)} association rules"
    )


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 65)
    print("GSE10072 GENE-EXPRESSION ASSOCIATION RULE MINING")
    print("=" * 65)

    # Load data.
    expression, metadata = read_geo_series_matrix(MATRIX_FILE)

    print(f"\nExpression matrix: {expression.shape[0]} probes x "
          f"{expression.shape[1]} samples")

    # Classify sample groups.
    labels = build_sample_labels(metadata, expression.columns)
    labels = labels.reindex(expression.columns)
    labels.index.name = "sample_id"

    labels.to_csv(RESULTS_DIR / "sample_labels.csv")

    print("\nSample groups:")
    print(labels.value_counts(dropna=False).to_string())

    known = labels.isin(["Tumor", "Normal"])
    if known.sum() < 2:
        raise ValueError(
            "Could not classify enough tumor/normal samples. "
            "Inspect GEO sample titles and metadata before continuing."
        )

    # Keep only samples that can be classified confidently.
    expression = expression.loc[:, known]
    labels = labels.loc[known]

    # Save the cleaned expression matrix.
    expression.to_csv(RESULTS_DIR / "expression_matrix.csv")

    # Select probes with the largest variance across the retained samples.
    n_probes = min(TOP_VARIABLE_PROBES, expression.shape[0])
    variances = expression.var(axis=1).sort_values(ascending=False)
    selected_probes = variances.head(n_probes).index.tolist()

    selected_expression = expression.loc[selected_probes]

    # Save the selected probes and their variability.
    probe_summary = pd.DataFrame({
        "probe_id": selected_probes,
        "variance": variances.loc[selected_probes].values,
        "mean_expression": selected_expression.mean(axis=1).values,
        "tumor_mean": selected_expression.loc[:, labels == "Tumor"].mean(axis=1).values,
        "normal_mean": selected_expression.loc[:, labels == "Normal"].mean(axis=1).values,
    })

    probe_summary["tumor_minus_normal_mean"] = (
        probe_summary["tumor_mean"] - probe_summary["normal_mean"]
    )

    probe_summary.to_csv(
        RESULTS_DIR / "selected_variable_probes.csv",
        index=False,
    )

    # Convert each sample into items such as 1007_s_at_HIGH or
    # 1007_s_at_LOW, using each probe's median across samples.
    medians = selected_expression.median(axis=1)
    transactions = []

    for sample_id in selected_expression.columns:
        items = []

        for probe_id in selected_probes:
            value = selected_expression.at[probe_id, sample_id]
            state = "HIGH" if value >= medians[probe_id] else "LOW"
            items.append(f"{probe_id}_{state}")

        transactions.append(items)

    # Save transactions as readable rows.
    transaction_frame = pd.DataFrame({
        "sample_id": selected_expression.columns,
        "group": labels.reindex(selected_expression.columns).values,
        "items": [";".join(items) for items in transactions],
    })
    transaction_frame.to_csv(
        RESULTS_DIR / "transactions.csv",
        index=False,
    )

    # Mine across all classified samples.
    save_mining_results(
        "all_samples",
        transactions,
        labels,
        selected_probes,
    )

    # Mine tumor and normal groups separately to explore group-specific
    # co-expression patterns. Each group is mined independently.
    sample_ids = list(selected_expression.columns)

    for group_name in ["Tumor", "Normal"]:
        group_ids = [
            sample_id for sample_id in sample_ids
            if labels.loc[sample_id] == group_name
        ]

        group_transactions = [
            transactions[sample_ids.index(sample_id)]
            for sample_id in group_ids
        ]

        save_mining_results(
            group_name.lower(),
            group_transactions,
            labels.loc[group_ids],
            selected_probes,
        )

    print("\nAnalysis complete. Output files are in:")
    print(RESULTS_DIR)

    print("\nImportant interpretation notes:")
    print("- HIGH/LOW means relative to each probe's median in this dataset.")
    print("- Probe IDs are not necessarily gene symbols.")
    print("- Association rules indicate co-occurrence, not causation.")
    print("- These exploratory rules are not statistical proof of biomarkers.")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"\nERROR: {error}", file=sys.stderr)
        sys.exit(1)
