import pandas as pd

df = pd.read_csv(
    "data/processed/bcr_qc_filtered.tsv",
    sep="\t"
)

print("=" * 60)
print("CLONE IDENTIFIER INSPECTION")
print("=" * 60)

# Number of distinct clone IDs per sample
print("\n1. UNIQUE CLONE IDs PER SAMPLE")
print(
    df.groupby("sample_id")["clone_id"]
      .nunique()
)

# Check whether clone IDs appear in multiple samples
clone_samples = (
    df.groupby("clone_id")["sample_id"]
      .nunique()
)

shared_ids = clone_samples[clone_samples > 1]

print("\n2. CLONE IDs PRESENT IN BOTH SAMPLES")
print(f"Number of shared IDs: {len(shared_ids)}")

# Number of records per clone
clone_sizes = (
    df.groupby(["sample_id", "clone_id"])
      .size()
      .reset_index(name="sequence_count")
)

print("\n3. TOP 10 SAMPLE-CLONE GROUPS")
print(
    clone_sizes.sort_values(
        "sequence_count",
        ascending=False
    ).head(10).to_string(index=False)
)

# Inspect V and J assignments within clone groups
print("\n4. V/J GENE CONSISTENCY")

gene_consistency = (
    df.groupby(["sample_id", "clone_id"])
      .agg(
          n_sequences=("sequence_id", "count"),
          n_v_genes=("v_call_genotyped", "nunique"),
          n_j_genes=("j_call", "nunique"),
          n_junctions=("junction", "nunique")
      )
      .reset_index()
)

print(
    gene_consistency.sort_values(
        "n_sequences",
        ascending=False
    ).head(10).to_string(index=False)
)

# Save results
clone_sizes.to_csv(
    "results/clone_size_summary.csv",
    index=False
)

gene_consistency.to_csv(
    "results/clone_gene_consistency.csv",
    index=False
)

print("\nClone ID inspection completed!")