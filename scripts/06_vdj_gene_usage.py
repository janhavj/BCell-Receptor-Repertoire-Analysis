from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

DATA = Path("data/processed/bcr_qc_filtered.tsv")
RESULTS = Path("results")
FIGURES = Path("figures")

df = pd.read_csv(DATA, sep="\t")

RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

# Extract standard immunoglobulin gene names.
# Example: "Homsap IGHV3-11*05 F" -> "IGHV3-11"
def extract_gene(value):
    if pd.isna(value):
        return None

    parts = str(value).split()

    for part in parts:
        if part.startswith(("IGHV", "IGHD", "IGHJ")):
            return part.split("*")[0]

    return None

df["v_gene"] = df["v_call"].apply(extract_gene)
df["j_gene"] = df["j_call"].apply(extract_gene)

print("=" * 60)
print("IGHV AND IGHJ GENE USAGE ANALYSIS")
print("=" * 60)

for gene_col in ["v_gene", "j_gene"]:
    print(f"\n{gene_col.upper()} USAGE")

    usage = (
        df.groupby("sample_id")[gene_col]
        .value_counts(normalize=True)
        .rename("frequency")
        .reset_index()
    )

    usage.to_csv(
        RESULTS / f"{gene_col}_usage.csv",
        index=False
    )

    for sample in ["-1h", "+7d"]:
        subset = (
            usage[usage["sample_id"] == sample]
            .sort_values("frequency", ascending=False)
            .head(10)
        )

        print(f"\nSample: {sample}")
        print(subset.to_string(index=False))

        plt.figure(figsize=(9, 5))

        sns.barplot(
            data=subset,
            x="frequency",
            y=gene_col,
            color="steelblue"
        )

        plt.title(f"Top 10 {gene_col} assignments: {sample}")
        plt.xlabel("Fraction of annotated records")
        plt.ylabel("Gene")
        plt.tight_layout()

        filename = (
            f"{gene_col}_usage_"
            f"{sample.replace('-', 'minus').replace('+', 'plus')}.png"
        )

        plt.savefig(
            FIGURES / filename,
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()

print("\nIGHV/IGHJ gene usage analysis completed!")