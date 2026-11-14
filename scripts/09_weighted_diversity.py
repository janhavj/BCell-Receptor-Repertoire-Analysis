from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import entropy

DATA = Path("data/processed/bcr_qc_filtered.tsv")
RESULTS = Path("results")
RESULTS.mkdir(exist_ok=True)

df = pd.read_csv(DATA, sep="\t")

# Validate abundance weights before using them
weights = pd.to_numeric(
    df["duplicate_count"],
    errors="coerce"
)

if weights.isna().any() or (weights <= 0).any():
    raise ValueError("Invalid duplicate_count values detected")

df["abundance_weight"] = weights

def calculate_metrics(counts):
    counts = counts[counts > 0]
    frequencies = counts / counts.sum()

    return {
        "richness": len(counts),
        "shannon": entropy(frequencies),
        "simpson": 1 - np.sum(frequencies ** 2),
        "largest_clone_fraction": frequencies.max()
    }

results = []

for sample, group in df.groupby("sample_id"):

    # Each annotated record contributes one observation
    record_counts = group.groupby("clone_id").size()

    # Each record contributes its duplicate_count weight
    weighted_counts = (
        group.groupby("clone_id")["abundance_weight"]
        .sum()
    )

    for method, counts in [
        ("record_based", record_counts),
        ("duplicate_count_weighted", weighted_counts)
    ]:
        metrics = calculate_metrics(counts)

        results.append({
            "sample_id": sample,
            "method": method,
            "total_observations": counts.sum(),
            **metrics
        })

summary = pd.DataFrame(results)

print("=" * 65)
print("RECORD-BASED VS ABUNDANCE-WEIGHTED DIVERSITY")
print("=" * 65)

print(
    summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)

summary.to_csv(
    RESULTS / "weighted_diversity_comparison.csv",
    index=False
)

print("\nResults saved successfully!")