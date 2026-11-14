import numpy as np
import pandas as pd
from scipy.stats import entropy
from pathlib import Path

DATA = Path("data/processed/bcr_qc_filtered.tsv")
RESULTS = Path("results")
RESULTS.mkdir(exist_ok=True)

df = pd.read_csv(DATA, sep="\t")

print("=" * 60)
print("BCR CLONOTYPE DIVERSITY ANALYSIS")
print("=" * 60)

summary = []
clone_tables = []

for sample, group in df.groupby("sample_id"):

    counts = group["clone_id"].value_counts()
    frequencies = counts / counts.sum()

    richness = len(counts)
    shannon = entropy(frequencies)
    simpson = 1 - np.sum(frequencies ** 2)
    largest_clone = frequencies.max()

    summary.append({
        "sample_id": sample,
        "n_records": len(group),
        "clone_richness": richness,
        "shannon_entropy": shannon,
        "simpson_diversity": simpson,
        "largest_clone_fraction": largest_clone
    })

    clone_table = pd.DataFrame({
        "clone_id": counts.index,
        "sequence_count": counts.values,
        "frequency": frequencies.values
    })

    clone_table.insert(0, "sample_id", sample)
    clone_tables.append(clone_table)

    print(f"\nSAMPLE: {sample}")
    print(f"Productive records: {len(group)}")
    print(f"Clone richness: {richness}")
    print(f"Shannon entropy: {shannon:.4f}")
    print(f"Simpson diversity: {simpson:.4f}")
    print(f"Largest-clone fraction: {largest_clone:.4f}")

summary_df = pd.DataFrame(summary)

summary_df.to_csv(
    RESULTS / "repertoire_diversity_summary.csv",
    index=False
)

pd.concat(clone_tables).to_csv(
    RESULTS / "clonotype_frequencies.csv",
    index=False
)

print("\nAnalysis completed successfully!")