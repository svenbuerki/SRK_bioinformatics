#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns

# ----------------------------
# User input
# ----------------------------
hap_key_tsv = input("Enter haplotype key TSV (from collapse_protein_haplotypes.py): ").strip()
aa_fasta = input("Enter aligned protein haplotype FASTA: ").strip()

# ----------------------------
# Load haplotype key
# ----------------------------
hap_df = pd.read_csv(hap_key_tsv, sep="\t")
print("Columns in the TSV:", hap_df.columns.tolist())

# Default column names
ind_col = 'Original_sequence_ID'
hap_col = 'Protein_haplotype'

# Prompt user if columns not found
if ind_col not in hap_df.columns:
    ind_col = input("Enter the column name for individual IDs: ").strip()
if hap_col not in hap_df.columns:
    hap_col = input("Enter the column name for haplotypes: ").strip()

# ----------------------------
# Build presence/absence matrix
# ----------------------------
individuals = hap_df[ind_col].tolist()
unique_haplotypes = sorted(hap_df[hap_col].unique().tolist())

binary_matrix = pd.DataFrame(0, index=individuals, columns=unique_haplotypes)

for _, row in hap_df.iterrows():
    binary_matrix.at[row[ind_col], row[hap_col]] = 1

# ----------------------------
# Determine optimal K using silhouette or elbow method (optional)
# Here we just pick sqrt(n_samples) as heuristic
# ----------------------------
n_individuals = binary_matrix.shape[0]
k = max(2, int(np.sqrt(n_individuals)))
print(f"Using k={k} for KMeans clustering")

# ----------------------------
# KMeans clustering
# ----------------------------
kmeans = KMeans(n_clusters=k, random_state=42)
clusters = kmeans.fit_predict(binary_matrix)

binary_matrix['Cluster'] = clusters

# ----------------------------
# PCA for visualization
# ----------------------------
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(binary_matrix.drop('Cluster', axis=1))
pca_df = pd.DataFrame(pca_coords, columns=['PC1', 'PC2'])
pca_df['Cluster'] = clusters
pca_df['Individual'] = binary_matrix.index

# ----------------------------
# Plot PCA colored by clusters
# ----------------------------
plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Cluster', palette='tab10', s=80)
plt.title('PCA of Individuals Colored by KMeans Clusters')
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
plt.legend(title='Cluster')
plt.tight_layout()
plt.savefig('PCA_clusters.png', dpi=300)
plt.show()

# ----------------------------
# Output cluster assignments
# ----------------------------
cluster_out = 'individual_clusters.tsv'
pca_df[['Individual','Cluster']].to_csv(cluster_out, sep='\t', index=False)
print(f"\nCluster assignments saved to: {cluster_out}")

# ----------------------------
# Optional: save binary matrix with clusters
# ----------------------------
binary_matrix.to_csv('binary_matrix_with_clusters.tsv', sep='\t')
print("Binary presence/absence matrix saved with cluster info: binary_matrix_with_clusters.tsv")

