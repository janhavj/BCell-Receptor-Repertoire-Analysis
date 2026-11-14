from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

RESULTS = Path("results")
FIGURES = Path("figures")
FIGURES.mkdir(exist_ok=True)

summary = pd.read_csv(
    RESULTS / "repertoire_diversity_summary.csv"
)

clones = pd.read_csv(
    RESULTS / "clonotype_frequencies.csv"
)

sns.set_theme(style="whitegrid", context="paper")

sample_order = ["-1h", "+7d"]

# Figure 1: Diversity metrics
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

metrics = [
    ("clone_richness", "Clone richness"),
    ("shannon_entropy", "Shannon entropy"),
    ("simpson_diversity", "Simpson diversity")
]

for ax, (column, title) in zip(axes, metrics):
    sns.barplot(
        data=summary,
        x="sample_id",
        y=column,
        order=sample_order,
        ax=ax
    )
    ax.set_title(title)
    ax.set_xlabel("Timepoint")
    ax.set_ylabel(title)

plt.tight_layout()
plt.savefig(
    FIGURES / "01_repertoire_diversity.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Figure 2: Top 20 clones per sample
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for ax, sample in zip(axes, sample_order):
    subset = (
        clones[clones["sample_id"] == sample]
        .nlargest(20, "sequence_count")
        .copy()
    )

    subset["clone_id"] = subset["clone_id"].astype(str)

    sns.barplot(
        data=subset,
        x="clone_id",
        y="sequence_count",
        ax=ax
    )

    ax.set_title(f"Top 20 clones: {sample}")
    ax.set_xlabel("Clone ID")
    ax.set_ylabel("Number of records")
    ax.tick_params(axis="x", rotation=90)

plt.tight_layout()
plt.savefig(
    FIGURES / "02_top_clonotypes.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Figure 3: Clone-size distribution
plt.figure(figsize=(8, 5))

sns.histplot(
    data=clones,
    x="sequence_count",
    hue="sample_id",
    hue_order=sample_order,
    bins=30,
    element="step",
    common_norm=False
)

plt.yscale("log")
plt.xlabel("Records per clone")
plt.ylabel("Number of clones (log scale)")
plt.title("BCR clone-size distribution")

plt.tight_layout()
plt.savefig(
    FIGURES / "03_clone_size_distribution.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("All three figures saved successfully!")

for file in sorted(FIGURES.glob("*.png")):
    print(file)