import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "src"
    / "gene_expression_apriori"
    / "sensitivity_thresholds.py"
)

SPEC = importlib.util.spec_from_file_location(
    "sensitivity_thresholds", SCRIPT
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_threshold_scenarios_are_stricter_than_or_equal_to_baseline():
    baseline_support = 0.20
    baseline_confidence = 0.60

    assert MODULE.SCENARIOS[0] == (
        "baseline", baseline_support, baseline_confidence
    )

    for _, support, confidence in MODULE.SCENARIOS:
        assert support >= baseline_support
        assert confidence >= baseline_confidence


def test_thresholds_never_increase_rule_count():
    rules = pd.DataFrame({
        "support": [0.20, 0.25, 0.30, 0.40],
        "confidence": [0.60, 0.70, 0.80, 0.90],
    })

    baseline_count = len(
        rules.loc[
            (rules["support"] >= 0.20)
            & (rules["confidence"] >= 0.60)
        ]
    )

    for _, min_support, min_confidence in MODULE.SCENARIOS:
        count = len(
            rules.loc[
                (rules["support"] >= min_support)
                & (rules["confidence"] >= min_confidence)
            ]
        )
        assert count <= baseline_count


def test_top_rule_limit_is_small():
    assert MODULE.TOP_N == 20
    assert MODULE.CHUNK_SIZE == 100_000
