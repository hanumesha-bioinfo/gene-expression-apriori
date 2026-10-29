
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "src"
    / "gene_expression_apriori"
    / "validate_sample_structure.py"
)


def test_sample_structure_validation_runs_successfully():
    result = subprocess.run(
        [sys.executable, str(SCRIPT)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "Total samples: 107" in result.stdout
    assert "Unique patients: 74" in result.stdout
    assert "Patients with both groups: 33" in result.stdout
    assert "Tumor-only patients: 25" in result.stdout
    assert "Normal-only patients: 16" in result.stdout
    assert "Validation completed." in result.stdout


def test_age_summary_is_reported():
    result = subprocess.run(
        [sys.executable, str(SCRIPT)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "AGE SUMMARY BY SAMPLE GROUP" in result.stdout
    assert "Normal" in result.stdout
    assert "Tumor" in result.stdout
    assert "65.59" in result.stdout
    assert "66.40" in result.stdout
