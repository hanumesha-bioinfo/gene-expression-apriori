from pathlib import Path

import GEOparse
import pandas as pd

root = Path(__file__).resolve().parents[2]
raw_dir = root / "data" / "raw"
output_dir = root / "data" / "processed"
results_dir = root / "results" / "tables"

output_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)

gpl = GEOparse.get_GEO(
    geo="GPL96",
    destdir=str(raw_dir),
)

annotation = gpl.table[
    ["ID", "Gene Symbol", "Gene Title", "ENTREZ_GENE_ID"]
].copy()

annotation["Gene Symbol"] = (
    annotation["Gene Symbol"].fillna("").astype(str).str.strip()
)

has_symbol = annotation["Gene Symbol"].ne("")
is_ambiguous = annotation["Gene Symbol"].str.contains(
    r"\s///\s", regex=True
)

# Keep only probes with one unambiguous gene symbol.
clean = annotation.loc[
    has_symbol & ~is_ambiguous
].copy()

clean["Gene Symbol"] = clean["Gene Symbol"].str.upper()

# Save every probe-to-gene annotation for traceability.
annotation.to_csv(
    output_dir / "gpl96_probe_annotation.csv",
    index=False,
)

# One row per probe, restricted to unambiguous gene mappings.
clean.to_csv(
    output_dir / "gpl96_clean_probe_gene_mapping.csv",
    index=False,
)

summary = pd.DataFrame(
    [
        ("total_platform_probes", len(annotation)),
        ("probes_without_gene_symbol", int((~has_symbol).sum())),
        ("probes_with_multiple_symbols", int(is_ambiguous.sum())),
        ("unambiguous_mapped_probes", len(clean)),
        ("unique_mapped_gene_symbols", clean["Gene Symbol"].nunique()),
    ],
    columns=["metric", "value"],
)

summary.to_csv(results_dir / "gene_mapping_summary.csv", index=False)

print("\n=== GENE MAPPING SUMMARY ===")
print(summary.to_string(index=False))

print("\nSaved:")
print(output_dir / "gpl96_probe_annotation.csv")
print(output_dir / "gpl96_clean_probe_gene_mapping.csv")
print(results_dir / "gene_mapping_summary.csv")
