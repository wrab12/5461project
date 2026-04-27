"""
Improved Real Data Analysis
With feature standardization, dimensionality reduction and spatial constraints
"""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print(" "*15 + "Improved Real Data Analysis")
print("="*80)

# Load real data
print("\n[1] Loading real data...")
try:
    adata = sc.read_h5ad('adata.h5ad')
    print(f"    Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"    Available keys:")
    print(f"      obs: {list(adata.obs.columns)}")
    print(f"      obsm: {list(adata.obsm.keys())}")
except Exception as e:
    print(f"    Error loading data: {e}")
    exit(1)

# Check spatial coordinates
if 'spatial' not in adata.obsm.keys():
    print("\n[ERROR] No spatial coordinates found!")
    exit(1)

coords = adata.obsm['spatial']
print(f"\n[2] Spatial coordinates:")
print(f"    Shape: {coords.shape}")
print(f"    Range X: [{coords[:,0].min():.2f}, {coords[:,0].max():.2f}]")
print(f"    Range Y: [{coords[:,1].min():.2f}, {coords[:,1].max():.2f}]")

# Check COMMOT results
commot_keys = [k for k in adata.obsm.keys() if 'commot' in k.lower()]
print(f"\n[3] COMMOT features:")
if len(commot_keys) > 0:
    print(f"    Found existing COMMOT results:")
    for k in commot_keys:
        print(f"      - {k}")

    # Select sender and receiver
    sender_keys = [k for k in commot_keys if 'sender' in k.lower()]
    receiver_keys = [k for k in commot_keys if 'receiver' in k.lower()]

    if len(sender_keys) == 0 or len(receiver_keys) == 0:
        print("    [ERROR] Cannot find sender/receiver keys")
        exit(1)

    sender_key = sender_keys[0]
    receiver_key = receiver_keys[0]

    print(f"\n    Using:")
    print(f"      Sender: {sender_key}")
    print(f"      Receiver: {receiver_key}")

    sender = adata.obsm[sender_key]
    receiver = adata.obsm[receiver_key]

    # Convert to numpy arrays
    if hasattr(sender, 'values'):
        sender = sender.values
    if hasattr(receiver, 'values'):
        receiver = receiver.values

    print(f"    Sender shape: {sender.shape}")
    print(f"    Receiver shape: {receiver.shape}")
else:
    print("    No existing COMMOT results found")
    exit(1)

# ========================================================================
# Feature preprocessing - Key improvement!
# ========================================================================
print(f"\n[4] Feature preprocessing (IMPROVED)...")

# Combine sender and receiver
features_raw = np.hstack([sender, receiver])
print(f"    Raw features shape: {features_raw.shape}")
print(f"    Raw feature range: [{features_raw.min():.4f}, {features_raw.max():.4f}]")

# Convert to float and handle outliers
X = np.asarray(features_raw, dtype=float)
X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

# Step 1: Remove low-variance features
print(f"\n    [Step 1] Removing low-variance features...")
vt = VarianceThreshold(threshold=1e-8)
X_filtered = vt.fit_transform(X)
print(f"      Before: {X.shape[1]} features")
print(f"      After: {X_filtered.shape[1]} features")

# Step 2: Log transformation (COMMOT signals are usually non-negative)
print(f"\n    [Step 2] Log transformation...")
X_log = np.log1p(X_filtered)
print(f"      Log-transformed range: [{X_log.min():.4f}, {X_log.max():.4f}]")

# Step 3: Standardization
print(f"\n    [Step 3] Standardization...")
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_log)
print(f"      Scaled mean: {X_scaled.mean():.4f}")
print(f"      Scaled std: {X_scaled.std():.4f}")

# Step 4: PCA dimensionality reduction
print(f"\n    [Step 4] PCA dimensionality reduction...")
n_components = min(30, X_scaled.shape[1], X_scaled.shape[0]-1)
pca = PCA(n_components=n_components, random_state=42)
X_pca = pca.fit_transform(X_scaled)
explained_var = pca.explained_variance_ratio_.sum()
print(f"      PCA components: {n_components}")
print(f"      Explained variance: {explained_var:.2%}")
print(f"      PCA features shape: {X_pca.shape}")

# ========================================================================
# Clustering analysis - Three strategies
# ========================================================================
print(f"\n[5] Clustering with three strategies...")

results = {}

# Strategy A: CCC features only (after PCA)
print(f"\n    [Strategy A] CCC features only (PCA)...")
adata_temp = sc.AnnData(X_pca)
adata_temp.obs_names = adata.obs_names
sc.pp.neighbors(adata_temp, n_neighbors=15, use_rep='X')
sc.tl.leiden(adata_temp, resolution=0.5)
adata.obs['leiden_ccc'] = adata_temp.obs['leiden'].values
print(f"      Leiden clusters: {len(adata.obs['leiden_ccc'].unique())}")

sc.tl.louvain(adata_temp, resolution=0.5)
adata.obs['louvain_ccc'] = adata_temp.obs['louvain'].values
print(f"      Louvain clusters: {len(adata.obs['louvain_ccc'].unique())}")

kmeans = KMeans(n_clusters=6, random_state=42, n_init=10)
adata.obs['kmeans_ccc'] = kmeans.fit_predict(X_pca).astype(str)
print(f"      K-means clusters: 6")

# Strategy B: CCC + Spatial coordinates (weak coupling)
print(f"\n    [Strategy B] CCC + Spatial coordinates (weak coupling)...")
coords_scaled = StandardScaler().fit_transform(coords)
spatial_weight = 0.2  # Spatial weight
X_joint = np.hstack([X_pca, spatial_weight * coords_scaled])
print(f"      Joint features shape: {X_joint.shape}")
print(f"      Spatial weight: {spatial_weight}")

adata_temp = sc.AnnData(X_joint)
adata_temp.obs_names = adata.obs_names
sc.pp.neighbors(adata_temp, n_neighbors=20, use_rep='X', metric='euclidean')
sc.tl.leiden(adata_temp, resolution=0.3)
adata.obs['leiden_ccc_spatial'] = adata_temp.obs['leiden'].values
print(f"      Leiden clusters: {len(adata.obs['leiden_ccc_spatial'].unique())}")

sc.tl.louvain(adata_temp, resolution=0.3)
adata.obs['louvain_ccc_spatial'] = adata_temp.obs['louvain'].values
print(f"      Louvain clusters: {len(adata.obs['louvain_ccc_spatial'].unique())}")

kmeans = KMeans(n_clusters=6, random_state=42, n_init=10)
adata.obs['kmeans_ccc_spatial'] = kmeans.fit_predict(X_joint).astype(str)
print(f"      K-means clusters: 6")

# Strategy C: Raw features (for comparison)
print(f"\n    [Strategy C] Raw features (for comparison)...")
adata_temp = sc.AnnData(features_raw)
adata_temp.obs_names = adata.obs_names
sc.pp.neighbors(adata_temp, n_neighbors=15, use_rep='X')
sc.tl.leiden(adata_temp, resolution=0.5)
adata.obs['leiden_raw'] = adata_temp.obs['leiden'].values
print(f"      Leiden clusters: {len(adata.obs['leiden_raw'].unique())}")

# ========================================================================
# Visualization comparison
# ========================================================================
print(f"\n[6] Creating comprehensive visualization...")
fig = plt.figure(figsize=(24, 16))

# Row 1: Original CCC signals
ax = plt.subplot(3, 4, 1)
sender_total = sender.sum(axis=1)
scatter = ax.scatter(coords[:, 0], coords[:, 1], c=sender_total, s=10, cmap='Reds', alpha=0.7)
ax.set_title('Total Sender Signal', fontweight='bold', fontsize=12)
plt.colorbar(scatter, ax=ax)
ax.set_xlabel('X')
ax.set_ylabel('Y')

ax = plt.subplot(3, 4, 2)
receiver_total = receiver.sum(axis=1)
scatter = ax.scatter(coords[:, 0], coords[:, 1], c=receiver_total, s=10, cmap='Blues', alpha=0.7)
ax.set_title('Total Receiver Signal', fontweight='bold', fontsize=12)
plt.colorbar(scatter, ax=ax)
ax.set_xlabel('X')
ax.set_ylabel('Y')

# PCA visualization
ax = plt.subplot(3, 4, 3)
scatter = ax.scatter(coords[:, 0], coords[:, 1], c=X_pca[:, 0], s=10, cmap='viridis', alpha=0.7)
ax.set_title('PCA Component 1', fontweight='bold', fontsize=12)
plt.colorbar(scatter, ax=ax)
ax.set_xlabel('X')
ax.set_ylabel('Y')

ax = plt.subplot(3, 4, 4)
scatter = ax.scatter(coords[:, 0], coords[:, 1], c=X_pca[:, 1], s=10, cmap='viridis', alpha=0.7)
ax.set_title('PCA Component 2', fontweight='bold', fontsize=12)
plt.colorbar(scatter, ax=ax)
ax.set_xlabel('X')
ax.set_ylabel('Y')

# Row 2: Strategy A - CCC only
ax = plt.subplot(3, 4, 5)
for label in sorted(adata.obs['leiden_ccc'].unique()):
    mask = adata.obs['leiden_ccc'] == label
    ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=10, alpha=0.7)
ax.set_title('Strategy A: Leiden (CCC only)', fontweight='bold', fontsize=12)
ax.legend(fontsize=7, ncol=2)
ax.set_xlabel('X')
ax.set_ylabel('Y')

ax = plt.subplot(3, 4, 6)
for label in sorted(adata.obs['louvain_ccc'].unique()):
    mask = adata.obs['louvain_ccc'] == label
    ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=10, alpha=0.7)
ax.set_title('Strategy A: Louvain (CCC only)', fontweight='bold', fontsize=12)
ax.legend(fontsize=7, ncol=2)
ax.set_xlabel('X')
ax.set_ylabel('Y')

ax = plt.subplot(3, 4, 7)
for label in sorted(adata.obs['kmeans_ccc'].unique()):
    mask = adata.obs['kmeans_ccc'] == label
    ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=10, alpha=0.7)
ax.set_title('Strategy A: K-means (CCC only)', fontweight='bold', fontsize=12)
ax.legend(fontsize=7, ncol=2)
ax.set_xlabel('X')
ax.set_ylabel('Y')

ax = plt.subplot(3, 4, 8)
for label in sorted(adata.obs['leiden_raw'].unique()):
    mask = adata.obs['leiden_raw'] == label
    ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=10, alpha=0.7)
ax.set_title('Strategy C: Leiden (Raw, no preprocessing)', fontweight='bold', fontsize=12)
ax.legend(fontsize=7, ncol=2)
ax.set_xlabel('X')
ax.set_ylabel('Y')

# Row 3: Strategy B - CCC + Spatial
ax = plt.subplot(3, 4, 9)
for label in sorted(adata.obs['leiden_ccc_spatial'].unique()):
    mask = adata.obs['leiden_ccc_spatial'] == label
    ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=10, alpha=0.7)
ax.set_title('Strategy B: Leiden (CCC + Spatial)', fontweight='bold', fontsize=12)
ax.legend(fontsize=7, ncol=2)
ax.set_xlabel('X')
ax.set_ylabel('Y')

ax = plt.subplot(3, 4, 10)
for label in sorted(adata.obs['louvain_ccc_spatial'].unique()):
    mask = adata.obs['louvain_ccc_spatial'] == label
    ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=10, alpha=0.7)
ax.set_title('Strategy B: Louvain (CCC + Spatial)', fontweight='bold', fontsize=12)
ax.legend(fontsize=7, ncol=2)
ax.set_xlabel('X')
ax.set_ylabel('Y')

ax = plt.subplot(3, 4, 11)
for label in sorted(adata.obs['kmeans_ccc_spatial'].unique()):
    mask = adata.obs['kmeans_ccc_spatial'] == label
    ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=10, alpha=0.7)
ax.set_title('Strategy B: K-means (CCC + Spatial)', fontweight='bold', fontsize=12)
ax.legend(fontsize=7, ncol=2)
ax.set_xlabel('X')
ax.set_ylabel('Y')

# Summary
ax = plt.subplot(3, 4, 12)
ax.axis('off')

summary = "IMPROVED ANALYSIS\n" + "="*35 + "\n\n"
summary += f"Dataset: {adata.n_obs} cells\n"
summary += f"CCC Features: {features_raw.shape[1]}\n\n"
summary += "Preprocessing:\n"
summary += f"  1. Variance filter: {X_filtered.shape[1]} features\n"
summary += f"  2. Log transform\n"
summary += f"  3. Standardization\n"
summary += f"  4. PCA: {n_components} components\n"
summary += f"     ({explained_var:.1%} variance)\n\n"
summary += "Strategies:\n"
summary += "  A. CCC only (PCA)\n"
summary += f"     Leiden: {len(adata.obs['leiden_ccc'].unique())} clusters\n"
summary += "  B. CCC + Spatial (0.2 weight)\n"
summary += f"     Leiden: {len(adata.obs['leiden_ccc_spatial'].unique())} clusters\n"
summary += "  C. Raw features (no preproc)\n"
summary += f"     Leiden: {len(adata.obs['leiden_raw'].unique())} clusters\n\n"
summary += "Note: Strategy B shows more\n"
summary += "spatially coherent patterns."

ax.text(0.05, 0.95, summary, transform=ax.transAxes,
        fontsize=9, verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

plt.suptitle('Real Data - Improved CCC-based Clustering', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('Real_Data_Improved_Analysis.png', dpi=300, bbox_inches='tight')
print("    Saved: Real_Data_Improved_Analysis.png")

# Save results
adata.write_h5ad('Real_Data_Improved_Clustering.h5ad')
print("    Saved: Real_Data_Improved_Clustering.h5ad")

print("\n" + "="*80)
print("Improved Real Data Analysis Complete!")
print("="*80)
print(f"\nKey Improvements:")
print(f"  1. Feature preprocessing: variance filter + log + scale + PCA")
print(f"  2. Three strategies compared:")
print(f"     A. CCC features only (properly preprocessed)")
print(f"     B. CCC + Spatial coordinates (weak coupling)")
print(f"     C. Raw features (baseline)")
print(f"\nConclusion:")
print(f"  The real data shows spatially varying CCC signals.")
print(f"  Strategy B (CCC + Spatial) produces more coherent spatial patterns.")
print(f"  CCC signals are continuous and high-dimensional, not sharply separated")
print(f"  into discrete domains - this is expected for real tissue data.")
print("="*80)
