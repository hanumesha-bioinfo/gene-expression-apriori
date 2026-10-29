"""Interactive research dashboard for gene-expression association rules."""

from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st


ROOT = Path(__file__).resolve().parents[1]
SUMMARY_FILE = ROOT / "results" / "sensitivity_thresholds" / "threshold_summary.csv"
RULES_FILE = ROOT / "results" / "sensitivity_thresholds" / "top_rules_by_threshold.csv"

st.set_page_config(
    page_title="Gene Expression Rule Explorer",
    page_icon="🧬",
    layout="wide",
)

SCENARIO_LABELS = {
    "baseline": "Baseline (support 0.20, confidence 0.60)",
    "support_025": "Support 0.25",
    "confidence_070": "Confidence 0.70",
    "confidence_080": "Confidence 0.80",
    "support_025_confidence_070": "Support 0.25 + confidence 0.70",
}


@st.cache_data
def load_data():
    """Load only the compact sensitivity outputs, not the large rules CSV."""
    summary = pd.read_csv(SUMMARY_FILE)
    rules = pd.read_csv(RULES_FILE)

    required_summary = {
        "scenario", "group", "min_support", "min_confidence",
        "rule_count", "baseline_rule_count",
        "retained_percent_of_baseline",
    }
    required_rules = {
        "scenario", "antecedents", "consequents", "support",
        "confidence", "lift", "group", "min_support", "min_confidence",
    }

    missing_summary = required_summary - set(summary.columns)
    missing_rules = required_rules - set(rules.columns)

    if missing_summary:
        raise ValueError(
            f"Summary CSV is missing columns: {sorted(missing_summary)}"
        )
    if missing_rules:
        raise ValueError(
            f"Top-rules CSV is missing columns: {sorted(missing_rules)}"
        )

    for column in [
        "min_support", "min_confidence", "rule_count",
        "baseline_rule_count", "retained_percent_of_baseline",
    ]:
        summary[column] = pd.to_numeric(summary[column], errors="coerce")

    for column in [
        "support", "confidence", "lift", "min_support", "min_confidence",
    ]:
        rules[column] = pd.to_numeric(rules[column], errors="coerce")

    return summary, rules


st.title("🧬 Gene Expression Association Rule Explorer")
st.caption(
    "Exploratory Apriori analysis of tumor and normal gene-expression groups"
)

if not SUMMARY_FILE.exists() or not RULES_FILE.exists():
    st.error("Required sensitivity-analysis CSV files were not found.")
    st.write("Expected files:")
    st.code(f"{SUMMARY_FILE}\n{RULES_FILE}")
    st.info(
        "Run the existing sensitivity-analysis workflow first. "
        "This dashboard does not rerun Apriori."
    )
    st.stop()

try:
    summary, rules = load_data()
except (ValueError, OSError, pd.errors.ParserError) as exc:
    st.error(f"Could not load dashboard data: {exc}")
    st.stop()

st.info(
    "Scientific interpretation: association rules describe co-occurrence "
    "patterns in the binarized expression data. They do not establish "
    "causation, statistical significance, or clinical validity."
)

# Overview metrics
baseline = summary[summary["scenario"] == "baseline"]
normal_baseline = baseline[baseline["group"] == "Normal"]
tumor_baseline = baseline[baseline["group"] == "Tumor"]

normal_count = (
    int(normal_baseline["rule_count"].iloc[0])
    if not normal_baseline.empty else 0
)
tumor_count = (
    int(tumor_baseline["rule_count"].iloc[0])
    if not tumor_baseline.empty else 0
)

col1, col2, col3 = st.columns(3)
col1.metric("Baseline rules — Normal", f"{normal_count:,}")
col2.metric("Baseline rules — Tumor", f"{tumor_count:,}")
col3.metric("Top-rule records loaded", f"{len(rules):,}")

st.caption(
    "The rule counts above come from the full baseline rule table. "
    "The explorer below contains only the saved top-rule sample "
    "(up to 20 rules per scenario and group)."
)

st.divider()

# Threshold sensitivity
st.header("1. Threshold sensitivity")

metric_choice = st.radio(
    "Chart metric",
    ["Rule count", "Retained percentage"],
    horizontal=True,
)

value_column = (
    "rule_count"
    if metric_choice == "Rule count"
    else "retained_percent_of_baseline"
)
axis_title = (
    "Number of rules"
    if metric_choice == "Rule count"
    else "Retained rules (% of baseline)"
)

chart_summary = summary.copy()
chart_summary["scenario_label"] = chart_summary["scenario"].map(
    SCENARIO_LABELS
).fillna(chart_summary["scenario"])

fig = px.bar(
    chart_summary,
    x="scenario_label",
    y=value_column,
    color="group",
    barmode="group",
    title=f"{metric_choice}: Normal vs Tumor",
    labels={
        "scenario_label": "Threshold scenario",
        value_column: axis_title,
        "group": "Expression group",
    },
    hover_data=["min_support", "min_confidence"],
)
fig.update_layout(xaxis_tickangle=-20)
st.plotly_chart(fig, width='stretch')

with st.expander("View threshold summary table"):
    st.dataframe(
        summary.sort_values(["scenario", "group"]),
        width='stretch',
        hide_index=True,
    )

st.divider()

# Interactive explorer
st.header("2. Explore saved top rules")

available_scenarios = [
    s for s in SCENARIO_LABELS if s in rules["scenario"].unique()
]
available_groups = sorted(rules["group"].dropna().unique().tolist())

with st.sidebar:
    st.header("Rule filters")

    chosen_scenarios = st.multiselect(
        "Threshold scenarios",
        options=available_scenarios,
        default=["baseline"] if "baseline" in available_scenarios
        else available_scenarios[:1],
        format_func=lambda s: SCENARIO_LABELS.get(s, s),
    )

    chosen_groups = st.multiselect(
        "Expression group",
        options=available_groups,
        default=available_groups,
    )

    gene_query = st.text_input(
        "Gene search",
        placeholder="e.g. TP53 or TFPI2",
        help="Searches both antecedents and consequents.",
    ).strip().upper()

    if not rules.empty:
        support_min = float(rules["support"].min())
        support_max = float(rules["support"].max())
        confidence_min = float(rules["confidence"].min())
        confidence_max = float(rules["confidence"].max())
        lift_min = float(rules["lift"].min())
        lift_max = float(rules["lift"].max())

        support_range = st.slider(
            "Support range",
            min_value=support_min,
            max_value=support_max,
            value=(support_min, support_max),
            format="%.3f",
        )
        confidence_range = st.slider(
            "Confidence range",
            min_value=confidence_min,
            max_value=confidence_max,
            value=(confidence_min, confidence_max),
            format="%.3f",
        )
        lift_floor = st.number_input(
            "Minimum lift",
            min_value=0.0,
            max_value=max(100.0, lift_max),
            value=max(0.0, lift_min),
            step=0.1,
        )
    else:
        support_range = (0.0, 1.0)
        confidence_range = (0.0, 1.0)
        lift_floor = 0.0

filtered = rules.copy()

if chosen_scenarios:
    filtered = filtered[filtered["scenario"].isin(chosen_scenarios)]
else:
    filtered = filtered.iloc[0:0]

if chosen_groups:
    filtered = filtered[filtered["group"].isin(chosen_groups)]
else:
    filtered = filtered.iloc[0:0]

filtered = filtered[
    filtered["support"].between(*support_range)
    & filtered["confidence"].between(*confidence_range)
    & (filtered["lift"] >= lift_floor)
]

if gene_query:
    gene_mask = (
        filtered["antecedents"].fillna("").str.upper().str.contains(
            gene_query, regex=False
        )
        | filtered["consequents"].fillna("").str.upper().str.contains(
            gene_query, regex=False
        )
    )
    filtered = filtered[gene_mask]

st.write(f"**Matching saved rules:** {len(filtered):,}")

if filtered.empty:
    st.warning(
        "No saved rules match these filters. Try selecting more scenarios "
        "or widening the metric ranges."
    )
else:
    sort_by = st.selectbox(
        "Sort rules by",
        ["lift", "confidence", "support"],
        index=0,
    )
    descending = st.checkbox("Highest values first", value=True)

    display_rules = filtered.sort_values(
        sort_by, ascending=not descending
    )

    fig_scatter = px.scatter(
        display_rules,
        x="support",
        y="confidence",
        color="group",
        size="lift",
        hover_data=["antecedents", "consequents", "scenario", "lift"],
        title="Support vs confidence",
        labels={
            "support": "Support",
            "confidence": "Confidence",
            "group": "Expression group",
        },
    )
    st.plotly_chart(fig_scatter, width='stretch')

    st.dataframe(
        display_rules[
            [
                "scenario", "group", "antecedents", "consequents",
                "support", "confidence", "lift",
            ]
        ],
        width='stretch',
        hide_index=True,
    )

    csv_bytes = display_rules.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download filtered rules as CSV",
        data=csv_bytes,
        file_name="filtered_gene_expression_rules.csv",
        mime="text/csv",
    )


st.divider()

# Dedicated top-rule comparison
st.header("3. Top rules by scenario")
st.caption(
    "Explore the saved top-rule sample for one threshold scenario and "
    "expression group. These are not all rules from the original dataset."
)

available_top_scenarios = [
    scenario
    for scenario in SCENARIO_LABELS
    if scenario in rules["scenario"].dropna().unique()
]

if available_top_scenarios:
    top_col1, top_col2 = st.columns(2)

    with top_col1:
        top_scenario = st.selectbox(
            "Choose threshold scenario",
            options=available_top_scenarios,
            format_func=lambda value: SCENARIO_LABELS.get(value, value),
            key="top_rule_scenario",
        )

    scenario_rules = rules[rules["scenario"] == top_scenario].copy()
    scenario_groups = sorted(
        scenario_rules["group"].dropna().unique().tolist()
    )

    with top_col2:
        top_group = st.selectbox(
            "Choose expression group",
            options=scenario_groups,
            key="top_rule_group",
        )

    selected_top_rules = scenario_rules[
        scenario_rules["group"] == top_group
    ].copy()

    metric_labels = {
        "lift": "Lift",
        "confidence": "Confidence",
        "support": "Support",
    }

    metric_to_rank = st.selectbox(
        "Rank rules by",
        options=list(metric_labels),
        format_func=lambda value: metric_labels[value],
        index=0,
        key="top_rule_rank_metric",
    )

    selected_top_rules = selected_top_rules.dropna(
        subset=[metric_to_rank]
    ).sort_values(metric_to_rank, ascending=False).head(20)

    if selected_top_rules.empty:
        st.warning("No saved top rules are available for this selection.")
    else:
        selected_top_rules["rule"] = (
            selected_top_rules["antecedents"].fillna("")
            + " → "
            + selected_top_rules["consequents"].fillna("")
        )

        chart_rules = selected_top_rules.sort_values(
            metric_to_rank, ascending=True
        )

        fig_top_rules = px.bar(
            chart_rules,
            x=metric_to_rank,
            y="rule",
            orientation="h",
            title=(
                f"Top saved {top_group} rules — "
                f"{metric_labels[metric_to_rank]}"
            ),
            labels={
                metric_to_rank: metric_labels[metric_to_rank],
                "rule": "Association rule",
            },
            hover_data=["support", "confidence", "lift"],
        )
        fig_top_rules.update_layout(
            yaxis_title="Antecedents → consequent",
            yaxis={"categoryorder": "total ascending"},
        )
        st.plotly_chart(fig_top_rules, width="stretch")

        st.dataframe(
            selected_top_rules[
                [
                    "scenario",
                    "group",
                    "antecedents",
                    "consequents",
                    "support",
                    "confidence",
                    "lift",
                ]
            ],
            width="stretch",
            hide_index=True,
        )

        st.download_button(
            "Download these top rules as CSV",
            data=selected_top_rules.drop(columns=["rule"]).to_csv(
                index=False
            ).encode("utf-8"),
            file_name=f"top_rules_{top_scenario}_{top_group.lower()}.csv",
            mime="text/csv",
            key="download_top_scenario_rules",
        )
else:
    st.warning("No recognized threshold scenarios were found in the rules file.")

st.divider()
st.header("4. Interpretation and limitations")
st.markdown(
    """
    - **Support** measures how frequently an itemset occurs in the transactions.
    - **Confidence** measures the proportion of antecedent-containing
      transactions that also contain the consequent.
    - **Lift** compares observed co-occurrence with the expectation under
      independence; lift above 1 indicates positive association in this dataset.
    - Threshold comparisons show how many baseline rules meet stricter
      thresholds. They do not test statistical significance.
    - The explorer uses a saved top-rule sample, not every rule in the
      original large CSV. A rule absent here may still exist in the full results.
    - Tumor and Normal rule counts should not be interpreted as a direct
      measure of biological importance. Patient pairing, confounding, and
      external validation require separate analyses.
    """
)
st.caption("Data source: local project sensitivity-analysis outputs.")
