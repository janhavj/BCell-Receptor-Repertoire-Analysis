import pandas as pd

df = pd.read_csv(
    "data/processed/bcr_qc_filtered.tsv",
    sep="\t"
)

def extract_gene(value):
    if pd.isna(value):
        return None

    for part in str(value).split():
        if part.startswith(("IGHV", "IGHD", "IGHJ")):
            return part.split("*")[0]

    return None

df["v_gene"] = df["v_call"].apply(extract_gene)
df["j_gene"] = df["j_call"].apply(extract_gene)

day7 = df[df["sample_id"] == "+7d"].copy()

print("=" * 60)
print("GENE USAGE VALIDATION")
print("=" * 60)

# IGHV3-49 records
v_subset = day7[day7["v_gene"] == "IGHV3-49"]

print("\n1. IGHV3-49 AT +7 DAYS")
print("Records:", len(v_subset))
print("Distinct clone IDs:", v_subset["clone_id"].nunique())
print("Distinct nucleotide junctions:", v_subset["junction"].nunique())

print("\n2. TOP CLONES ASSOCIATED WITH IGHV3-49")
print(
    v_subset["clone_id"]
    .value_counts()
    .head(15)
    .to_string()
)

print("\n3. IGHJ USAGE WITHIN IGHV3-49 RECORDS")
print(
    v_subset["j_gene"]
    .value_counts()
    .to_string()
)

# Combined V/J usage
print("\n4. TOP V/J GENE PAIRS AT +7 DAYS")
pairs = (
    day7.groupby(["v_gene", "j_gene"])
    .size()
    .sort_values(ascending=False)
    .head(15)
)

print(pairs.to_string())

# Examine nucleotide junction diversity
print("\n5. MOST FREQUENT NUCLEOTIDE JUNCTIONS")
print(
    day7["junction"]
    .value_counts()
    .head(10)
    .to_string()
)

# Cross-timepoint clone ID validation
print("\n6. SHARED CLONE IDs ACROSS TIMEPOINTS")

shared_ids = set(
    df.loc[df["sample_id"] == "-1h", "clone_id"]
).intersection(
    df.loc[df["sample_id"] == "+7d", "clone_id"]
)

for clone_id in sorted(shared_ids):
    subset = df[df["clone_id"] == clone_id]

    print(f"\nClone ID: {clone_id}")

    for sample, group in subset.groupby("sample_id"):
        print(
            f"  {sample}: {len(group)} records, "
            f"{group['junction'].nunique()} distinct junctions"
        )

    junction_sets = {
        sample: set(group["junction"])
        for sample, group in subset.groupby("sample_id")
    }

    overlap = junction_sets["-1h"] & junction_sets["+7d"]

    print("  Shared nucleotide junctions:", len(overlap))

print("\nValidation completed!")