import pandas as pd

df = pd.read_csv(
    "data/processed/bcr_qc_filtered.tsv",
    sep="\t"
)

print("=" * 60)
print("BCR ABUNDANCE AND CLONE CONSISTENCY CHECK")
print("=" * 60)

print("\n1. DUPLICATE COUNT DISTRIBUTION")
print(
    df["duplicate_count"]
    .value_counts(dropna=False)
    .head(20)
    .to_string()
)

print("\n2. DUPLICATE COUNT SUMMARY BY SAMPLE")
print(
    df.groupby("sample_id")["duplicate_count"]
    .agg(["count", "sum", "min", "median", "max"])
    .to_string()
)

print("\n3. CLONE IDs WITH MULTIPLE JUNCTIONS")

clone_junctions = (
    df.groupby(["sample_id", "clone_id"])["junction"]
    .nunique()
)

print(
    "Total sample-clone groups:",
    len(clone_junctions)
)

print(
    "Groups with >1 junction:",
    int((clone_junctions > 1).sum())
)

print("\n4. JUNCTIONS ASSIGNED TO MULTIPLE CLONE IDs")

junction_clones = (
    df.groupby(["sample_id", "junction"])["clone_id"]
    .nunique()
)

print(
    "Junctions with >1 clone ID:",
    int((junction_clones > 1).sum())
)

print("\n5. CLONE SIZE DISTRIBUTION")

clone_sizes = (
    df.groupby(["sample_id", "clone_id"])
    .size()
    .rename("records")
    .reset_index()
)

for sample, group in clone_sizes.groupby("sample_id"):
    print(f"\nSample: {sample}")
    print(
        group["records"]
        .describe()
        .to_string()
    )

    print(
        "Singleton clones:",
        int((group["records"] == 1).sum())
    )

    print(
        "Clones with >=10 records:",
        int((group["records"] >= 10).sum())
    )

print("\nAbundance and consistency checks completed!")