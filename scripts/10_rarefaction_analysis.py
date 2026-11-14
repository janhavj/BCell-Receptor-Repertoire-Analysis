from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import entropy

DATA = Path("data/processed/bcr_qc_filtered.tsv")
RESULTS = Path("results")
RESULTS.mkdir(exist_ok=True)

df = pd.read_csv(DATA, sep="\t")

baseline = df[df["sample_id"] == "-1h"]
day7 = df[df["sample_id"] == "+7d"]

depth = min(len(baseline), len(day7))
iterations = 1000
rng = np.random.default_rng(42)

def metrics(clone_ids):
    counts = pd.Series(clone_ids).value_counts()
    frequencies = counts / counts.sum()

    return {
        "richness": len(counts),
        "shannon": entropy(frequencies),
        "simpson": 1 - np.sum(frequencies ** 2),
        "largest_clone_fraction": frequencies.max()
    }

results = []

for sample, group in [
    ("-1h", baseline),
    ("+7d", day7)
]:
    ids = group["clone_id"].to_numpy()

    for iteration in range(iterations):
        selected = rng.choice(
            ids,
            size=depth,
            replace=False
        )

        result = metrics(selected)

        results.append({
            "sample_id": sample,
            "iteration": iteration + 1,
            **result
        })

rarefaction = pd.DataFrame(results)

summary = (
    rarefaction.groupby("sample_id")
    .agg(
        richness_mean=("richness", "mean"),
        richness_sd=("richness", "std"),
        shannon_mean=("shannon", "mean"),
        shannon_sd=("shannon", "std"),
        simpson_mean=("simpson", "mean"),
        largest_clone_mean=(
            "largest_clone_fraction", "mean"
        )
    )
)

print("=" * 65)
print("EQUAL-DEPTH BCR REPERTOIRE COMPARISON")
print("=" * 65)

print(f"Subsampling depth: {depth}")
print(f"Iterations: {iterations}")

print("\nRESULTS")
print(summary.to_string(float_format=lambda x: f"{x:.4f}"))

rarefaction.to_csv(
    RESULTS / "rarefaction_iterations.csv",
    index=False
)

summary.to_csv(
    RESULTS / "rarefaction_summary.csv"
)

print("\nRarefaction completed successfully!")