from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

RESULTS = Path("results")
FIGURES = Path("figures")
FIGURES.mkdir(exist_ok=True)

diversity = pd.read_csv(
    RESULTS / "weighted_diversity_comparison.csv"
)

v_usage = pd.read_csv(
    RESULTS / "v_gene_usage.csv"
)

j_usage = pd.read_csv(
    RESULTS / "j_gene_usage.csv"
)

order = ["-1h", "+7d"]

record_diversity = diversity[
    diversity["method"] == "record_based"
]

weighted_diversity = diversity[
    diversity["method"] == "duplicate_count_weighted"
]

sns.set_theme(style="whitegrid", context="talk")

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# Panel A: Shannon diversity
sns.barplot(
    data=record_diversity,
    x="sample_id",
    y="shannon",
    order=order,
    ax=axes[0, 0]
)
axes[0, 0].set_title("A. Shannon diversity")
axes[0, 0].set_xlabel("Timepoint")
axes[0, 0].set_ylabel("Shannon entropy")

# Panel B: Weighted dominant-clone fraction
sns.barplot(
    data=weighted_diversity,
    x="sample_id",
    y="largest_clone_fraction",
    order=order,
    ax=axes[0, 1]
)
axes[0, 1].set_title("B. Dominant-clone fraction (weighted)")
axes[0, 1].set_xlabel("Timepoint")
axes[0, 1].set_ylabel("Fraction of duplicate-count weight")

# Panel C: IGHV3-49 usage
v_subset = (
    v_usage[v_usage["v_gene"] == "IGHV3-49"]
    .set_index("sample_id")[["frequency"]]
    .reindex(order)
    .fillna({"frequency": 0})
    .reset_index()
)

sns.barplot(
    data=v_subset,
    x="sample_id",
    y="frequency",
    order=order,
    ax=axes[1, 0]
)
axes[1, 0].set_title("C. IGHV3-49 usage")
axes[1, 0].set_xlabel("Timepoint")
axes[1, 0].set_ylabel("Fraction of records")

# Panel D: IGHJ5 usage
# Panel D: IGHJ5 usage
j_subset = (
    j_usage[j_usage["j_gene"] == "IGHJ5"]
    .set_index("sample_id")[["frequency"]]
    .reindex(order)
    .fillna({"frequency": 0})
    .reset_index()
)

sns.barplot(
    data=j_subset,
    x="sample_id",
    y="frequency",
    order=order,
    ax=axes[1, 1]
)
axes[1, 1].set_title("D. IGHJ5 usage")
axes[1, 1].set_xlabel("Timepoint")
axes[1, 1].set_ylabel("Fraction of records")

plt.suptitle(
    "Human BCR Repertoire Analysis: Baseline vs Day 7",
    fontsize=17,
    fontweight="bold"
)

plt.tight_layout()

output = FIGURES / "05_integrated_bcr_summary.png"

plt.savefig(
    output,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Integrated figure saved:", output)