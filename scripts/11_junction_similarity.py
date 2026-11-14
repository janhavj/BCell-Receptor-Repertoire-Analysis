from pathlib import Path
from Bio import Align
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

DATA = Path("data/processed/bcr_qc_filtered.tsv")
RESULTS = Path("results")
FIGURES = Path("figures")

RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

df = pd.read_csv(DATA, sep="\t")

# Focus on the +7-day IGHV3-49 repertoire
day7 = df[df["sample_id"] == "+7d"].copy()

day7 = day7[
    day7["v_call"].str.contains(
        r"\bIGHV3-49\*",
        regex=True,
        na=False
    )
]

# Select the 20 most frequent nucleotide junctions
top = (
    day7["junction"]
    .value_counts()
    .head(20)
)

sequences = top.index.tolist()
counts = top.values

# Global alignment with unit mismatch/gap penalties
aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 0
aligner.mismatch_score = -1
aligner.open_gap_score = -1
aligner.extend_gap_score = -1

n = len(sequences)
distance_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(i + 1, n):
        score = aligner.score(
            sequences[i],
            sequences[j]
        )

        # With the scoring scheme above, the negative
        # alignment score is an edit-distance equivalent.
        distance = -score / max(
            len(sequences[i]),
            len(sequences[j])
        )

        distance_matrix[i, j] = distance
        distance_matrix[j, i] = distance

labels = [
    f"J{i+1} (n={count})"
    for i, count in enumerate(counts)
]

distance_df = pd.DataFrame(
    distance_matrix,
    index=labels,
    columns=labels
)

distance_df.to_csv(
    RESULTS / "junction_distance_matrix.csv"
)

plt.figure(figsize=(11, 9))

sns.heatmap(
    distance_df,
    cmap="viridis",
    square=True,
    linewidths=0.3,
    cbar_kws={
        "label": "Normalised edit distance"
    }
)

plt.title(
    "Top 20 IGHV3-49 nucleotide junctions (+7 days)"
)

plt.tight_layout()

plt.savefig(
    FIGURES / "04_junction_similarity_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

upper_triangle = distance_matrix[
    np.triu_indices(n, k=1)
]

print("=" * 60)
print("BCR JUNCTION SIMILARITY ANALYSIS")
print("=" * 60)

print("IGHV3-49 records:", len(day7))
print("Unique junctions:", day7["junction"].nunique())
print("Junctions compared:", n)

print(
    "Mean normalised distance:",
    round(upper_triangle.mean(), 4)
)

print(
    "Minimum distance:",
    round(upper_triangle.min(), 4)
)

print(
    "Maximum distance:",
    round(upper_triangle.max(), 4)
)

print("\nSimilarity heatmap saved successfully!")