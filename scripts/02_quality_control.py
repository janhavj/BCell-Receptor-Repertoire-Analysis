from pathlib import Path
import pandas as pd

INPUT = Path("data/raw/bcr_repertoire.tsv.gz")
OUTPUT = Path("data/processed")
RESULTS = Path("results")

OUTPUT.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT, sep="\t", compression="gzip")

print("=" * 60)
print("BCR REPERTOIRE QUALITY CONTROL")
print("=" * 60)

print(f"\nInitial records: {len(df):,}")

# Check productive flag
print("\n1. PRODUCTIVITY")
print(df["productive"].value_counts(dropna=False))

# Check immunoglobulin locus
print("\n2. RECEPTOR LOCUS")
print(df["locus"].value_counts(dropna=False))

# Check sample distribution
print("\n3. SAMPLE DISTRIBUTION")
print(df["sample_id"].value_counts(dropna=False))

# Check nucleotide junctions
print("\n4. JUNCTION QUALITY")
print("Missing junctions:", df["junction"].isna().sum())
print("Missing clone IDs:", df["clone_id"].isna().sum())

# Check junction length consistency
observed_length = df["junction"].str.len()
reported_length = pd.to_numeric(
    df["junction_length"], errors="coerce"
)

mismatch = (
    observed_length.notna()
    & reported_length.notna()
    & (observed_length != reported_length)
)

print("Junction length mismatches:", mismatch.sum())

# Standardise productive values
productive = (
    df["productive"]
    .astype("string")
    .str.strip()
    .str.upper()
    .isin(["T", "TRUE", "1"])
)

# Filter records suitable for analysis
qc = df.loc[
    productive
    & df["locus"].eq("IGH")
    & df["junction"].notna()
    & df["clone_id"].notna()
    & ~mismatch
].copy()

print("\n5. QC RESULTS")
print(f"Records before QC: {len(df):,}")
print(f"Records after QC: {len(qc):,}")
print(f"Records excluded: {len(df) - len(qc):,}")

print("\n6. POST-QC SAMPLE COUNTS")
print(qc["sample_id"].value_counts())

# Save processed records
qc.to_csv(
    OUTPUT / "bcr_qc_filtered.tsv",
    sep="\t",
    index=False
)

summary = pd.DataFrame({
    "metric": [
        "initial_records",
        "retained_records",
        "excluded_records",
        "junction_length_mismatches"
    ],
    "value": [
        len(df),
        len(qc),
        len(df) - len(qc),
        int(mismatch.sum())
    ]
})

summary.to_csv(
    RESULTS / "qc_summary.csv",
    index=False
)

print("\nQC files saved successfully!")