
"""Evaluate stricter Apriori thresholds from the existing baseline rule table.

Reads the existing rule CSV in chunks. Does not rerun Apriori or modify
existing result files. Outputs compact threshold summaries and top rules.
"""

from pathlib import Path

import numpy as np
import pandas as pd


INPUT = Path("results/tables/association_rules_by_group.csv")
OUTPUT_DIR = Path("results/sensitivity_thresholds")
CHUNK_SIZE = 100_000
TOP_N = 20

SCENARIOS = [
    ("baseline", 0.20, 0.60),
    ("support_025", 0.25, 0.60),
    ("confidence_070", 0.20, 0.70),
    ("confidence_080", 0.20, 0.80),
    ("support_025_confidence_070", 0.25, 0.70),
]

COLUMNS = [
    "antecedents",
    "consequents",
    "support",
    "confidence",
    "lift",
    "group",
]


def main():
    """Summarize stricter thresholds using the existing baseline rules."""

    if not INPUT.exists():
        raise FileNotFoundError(
            f"Input rule table not found: {INPUT}"
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    summary_path = OUTPUT_DIR / "threshold_summary.csv"
    top_path = OUTPUT_DIR / "top_rules_by_threshold.csv"

    # Never overwrite either existing output. This prevents mixed or stale
    # summary and top-rule files from being mistaken for one analysis run.
    existing_outputs = [
        path for path in (summary_path, top_path) if path.exists()
    ]
    if existing_outputs:
        paths = ", ".join(str(path) for path in existing_outputs)
        raise FileExistsError(
            f"Output file(s) already exist: {paths}. "
            "Move or back up the existing files before regenerating outputs."
        )

    summary_counts = {}
    baseline_counts = {}
    top_rules = {}

    for scenario, min_support, min_confidence in SCENARIOS:
        for group in ("Normal", "Tumor"):
            key = (scenario, group)
            summary_counts[key] = 0
            top_rules[key] = []

    # Process the large CSV in chunks; do not load the entire table at once.
    for chunk in pd.read_csv(INPUT, chunksize=CHUNK_SIZE):
        missing = set(COLUMNS) - set(chunk.columns)
        if missing:
            raise ValueError(
                "Input CSV is missing required columns: "
                + ", ".join(sorted(missing))
            )

        chunk["support"] = pd.to_numeric(
            chunk["support"], errors="coerce"
        )
        chunk["confidence"] = pd.to_numeric(
            chunk["confidence"], errors="coerce"
        )
        chunk["lift"] = pd.to_numeric(
            chunk["lift"], errors="coerce"
        )

        for group in ("Normal", "Tumor"):
            group_rows = chunk.loc[chunk["group"] == group]

            baseline_counts[group] = (
                baseline_counts.get(group, 0) + len(group_rows)
            )

            for scenario, min_support, min_confidence in SCENARIOS:
                key = (scenario, group)

                eligible = group_rows.loc[
                    (group_rows["support"] >= min_support)
                    & (group_rows["confidence"] >= min_confidence)
                ]

                summary_counts[key] += len(eligible)

                # Keep only finite-lift candidates for the top-rule table.
                candidates = eligible.loc[
                    np.isfinite(eligible["lift"])
                ].copy()

                if candidates.empty:
                    continue

                candidates = candidates.sort_values(
                    ["lift", "confidence", "support"],
                    ascending=[False, False, False],
                ).head(TOP_N)

                top_rules[key].append(candidates[COLUMNS].copy())

    # Build a compact summary for every scenario and group.
    summary_rows = []

    for group in ("Normal", "Tumor"):
        baseline_count = baseline_counts.get(group, 0)

        for scenario, min_support, min_confidence in SCENARIOS:
            rule_count = summary_counts[(scenario, group)]

            retained_percent = (
                100.0 * rule_count / baseline_count
                if baseline_count
                else 0.0
            )

            summary_rows.append(
                {
                    "scenario": scenario,
                    "group": group,
                    "min_support": min_support,
                    "min_confidence": min_confidence,
                    "rule_count": rule_count,
                    "baseline_rule_count": baseline_count,
                    "retained_percent_of_baseline": round(
                        retained_percent, 2
                    ),
                }
            )

    summary = pd.DataFrame(summary_rows)

    # Prepare the top rules for each scenario and group.
    top_rows = []

    for scenario, min_support, min_confidence in SCENARIOS:
        for group in ("Normal", "Tumor"):
            key = (scenario, group)
            frames = top_rules[key]

            if not frames:
                continue

            rules = pd.concat(frames, ignore_index=True)
            rules = rules.sort_values(
                ["lift", "confidence", "support"],
                ascending=[False, False, False],
            ).head(TOP_N)

            rules.insert(0, "scenario", scenario)
            rules["min_support"] = min_support
            rules["min_confidence"] = min_confidence

            top_rows.append(rules)

    # Write both outputs only after the analysis and summaries are ready.
    summary.to_csv(summary_path, index=False)

    output_columns = [
        "scenario",
        "antecedents",
        "consequents",
        "support",
        "confidence",
        "lift",
        "group",
        "min_support",
        "min_confidence",
    ]

    if top_rows:
        top_table = pd.concat(top_rows, ignore_index=True)
        top_table = top_table[output_columns]
    else:
        top_table = pd.DataFrame(columns=output_columns)

    top_table.to_csv(top_path, index=False)

    print("Sensitivity analysis completed.")
    print(f"Summary: {summary_path}")
    print(f"Top rules: {top_path}")
    print(f"Summary rows: {len(summary)}")
    print(f"Top-rule rows: {len(top_table)}")
    print()
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
