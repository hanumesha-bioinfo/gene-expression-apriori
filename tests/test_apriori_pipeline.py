from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"


def test_expression_matrix_has_expected_samples_and_no_missing_values():
    matrix = pd.read_csv(
        DATA / "gene_expression_matrix.csv",
        index_col=0,
    )

    assert matrix.shape[1] == 107
    assert matrix.shape[0] > 0
    assert not matrix.isna().any().any()


def test_sample_metadata_has_expected_groups():
    metadata = pd.read_csv(DATA / "sample_metadata.csv")

    assert len(metadata) == 107
    assert metadata["GSM_ID"].is_unique
    assert metadata["group"].value_counts().to_dict() == {
        "Tumor": 58,
        "Normal": 49,
    }


def test_binary_matrix_contains_only_binary_values():
    matrix = pd.read_csv(
        DATA / "apriori_binary_matrix.csv",
        index_col=0,
    )

    assert matrix.shape == (107, 200)
    assert not matrix.isna().any().any()

    values = set(np.unique(matrix.to_numpy()))
    assert values.issubset({0, 1, False, True})


def test_apriori_sample_groups_match_binary_matrix_samples():
    binary = pd.read_csv(
        DATA / "apriori_binary_matrix.csv",
        index_col=0,
    )
    groups = pd.read_csv(DATA / "apriori_sample_groups.csv")

    assert len(groups) == 107
    assert groups["GSM_ID"].is_unique
    assert set(groups["GSM_ID"]) == set(binary.index.astype(str))
    assert set(groups["group"]) == {"Tumor", "Normal"}


def test_filtered_rules_meet_selected_thresholds():
    rules = pd.read_csv(TABLES / "association_rules_filtered.csv")

    assert len(rules) > 0
    assert {"group", "support", "confidence", "lift"}.issubset(rules.columns)
    assert rules["group"].isin(["Tumor", "Normal"]).all()

    for column in ["support", "confidence", "lift"]:
        assert np.isfinite(rules[column]).all()

    assert rules["support"].ge(0.25).all()
    assert rules["confidence"].ge(0.70).all()
    assert rules["lift"].ge(1.20).all()
    assert rules["support"].le(1).all()
    assert rules["confidence"].le(1).all()
    assert rules["lift"].gt(0).all()


def test_visualization_outputs_exist():
    assert (FIGURES / "filtered_rules_by_group.png").is_file()
    assert (FIGURES / "top_association_rules_by_lift.png").is_file()
    assert (TABLES / "rule_summary_by_group.csv").is_file()


def test_rule_summary_contains_both_groups():
    summary = pd.read_csv(TABLES / "rule_summary_by_group.csv")

    assert set(summary["group"]) == {"Tumor", "Normal"}
    assert summary["filtered_rule_count"].gt(0).all()
