import pandas as pd

# Load processed BCR repertoire data
file_path = "data/raw/bcr_repertoire.tsv.gz"

df = pd.read_csv(
    file_path,
    sep="\t",
    compression="gzip",
    low_memory=False
)

print("\n" + "=" * 60)
print("B-CELL RECEPTOR REPERTOIRE DATASET INSPECTION")
print("=" * 60)

# 1. Dataset dimensions
print("\n1. DATASET SHAPE")
print(f"Rows: {df.shape[0]:,}")
print(f"Columns: {df.shape[1]}")

# 2. Column names
print("\n2. COLUMN NAMES")
for col in df.columns:
    print(f" - {col}")

# 3. First records
print("\n3. FIRST FIVE RECORDS")
print(df.head().to_string())

# 4. Sample information
print("\n4. SAMPLE INFORMATION")
if "sample_id" in df.columns:
    print(df["sample_id"].value_counts(dropna=False))
else:
    print("No sample_id column found.")

# 5. Immunogenetic annotation columns
print("\n5. IMMUNOGENETIC ANNOTATIONS")

columns_to_check = [
    "v_call",
    "d_call",
    "j_call",
    "junction_aa",
    "clone_id",
    "c_call"
]

for col in columns_to_check:
    if col in df.columns:
        print(f"{col}: {df[col].notna().sum():,} non-missing")
    else:
        print(f"{col}: NOT PRESENT")

# 6. Missing values
print("\n6. TOP 15 COLUMNS WITH MISSING VALUES")
print(
    df.isna()
      .sum()
      .sort_values(ascending=False)
      .head(15)
)

# 7. Duplicate records
print("\n7. DUPLICATE ROWS")
print(df.duplicated().sum())

print("\nInspection completed successfully!")