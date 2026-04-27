"""
Analyze Real Dataset (from Basic_usage.ipynb)
"""
import numpy as np
import pandas as pd
import scanpy as sc
import commot as ct
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print(" "*20 + "Real Data Analysis")
print("="*80)

# Load Real data
print("\n[1] Loading real data...")
try:
    adata = sc.read_h5ad('adata.h5ad')
    print(f"    Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"    Available keys:")
    print(f"      obs: {list(adata.obs.columns)}")
    print(f"      obsm: {list(adata.obsm.keys())}")
    print(f"      uns: {list(adata.uns.keys())}")
except Exception as e:
    print(f"    Error loading data: {e}")
    exit(1)

# Check if has spatial coordinates
if 'spatial' not in adata.obsm.keys():
    print("\n[ERROR] No spatial coordinates found!")
    print("Available obsm keys:", list(adata.obsm.keys()))
    exit(1)

coords = adata.obsm['spatial']
print(f"\n[2] Spatial coordinates:")
print(f"    Shape: {coords.shape}")
print(f"    Range X: [{coords[:,0].min():.2f}, {coords[:,0].max():.2f}]")
print(f"    Range Y: [{coords[:,1].min():.2f}, {coords[:,1].max():.2f}]")

# Check if already has COMMOT Results
commot_keys = [k for k in adata.obsm.keys() if 'commot' in k.lower()]
print(f"\n[3] COMMOT features:")
if len(commot_keys) > 0:
    print(f"    Found existing COMMOT results: {commot_keys}")

    # Use already existing COMMOT Results
    sender_key = [k for k in commot_keys if 'sender' in k.lower()][0]
    receiver_key = [k for k in commot_keys if 'receiver' in k.lower()][0]

    sender = adata.obsm[sender_key]
    receiver = adata.obsm[receiver_key]

    # Convert to numpy arrays
    if hasattr(sender, 'values'):
        sender = sender.values
    if hasattr(receiver, 'values'):
        receiver = receiver.values

    features = np.hstack([sender, receiver])
    print(f"    Sender shape: {sender.shape}")
    print(f"    Receiver shape: {receiver.shape}")
    print(f"    Combined features: {features.shape}")

else:
    print("    No existing COMMOT results found")
    print("    Running COMMOT analysis...")

    # Preprocessing
    print("\n[4] Preprocessing...")
    # Basic QC
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"    After QC: {adata.n_obs} cells, {adata.n_vars} genes")

    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Running COMMOT
    print("\n[5] Running COMMOT...")

    # Use simplified Ligand-Receptor pairs
    df_ligrec = pd.DataFrame({
        'ligand': ['Tgfb1', 'Vegfa', 'Fgf2', 'Pdgfa'],
        'receptor': ['Tgfbr1', 'Flt1', 'Fgfr1', 'Pdgfra'],
        'pathway': ['TGFb', 'VEGF', 'FGF', 'PDGF']
    })

    print(f"    Using {len(df_ligrec)} ligand-receptor pairs")

    try:
        ct.tl.spatial_communication(
            adata,
            database_name='real',
            df_ligrec=df_ligrec,
            dis_thr=200,  # Real data may need larger distance threshold
            heteromeric=False,
            pathway_sum=False
        )
        print("    COMMOT SUCCESS!")

        # Extracting features
        commot_keys = [k for k in adata.obsm.keys() if 'commot' in k.lower()]
        sender_key = [k for k in commot_keys if 'sender' in k.lower()][0]
        receiver_key = [k for k in commot_keys if 'receiver' in k.lower()][0]

        sender = adata.obsm[sender_key].values
        receiver = adata.obsm[receiver_key].values
        features = np.hstack([sender, receiver])

    except Exception as e:
        print(f"    COMMOT failed: {e}")
        print("    This is expected for real data without proper ligand-receptor database")
        exit(1)

# Clustering Analysis
print(f"\n[6] Clustering on CCC features...")
print(f"    Feature dimension: {features.shape}")

# Leiden
print("    Running Leiden...")
adata_temp = sc.AnnData(features)
adata_temp.obs_names = adata.obs_names
sc.pp.neighbors(adata_temp, n_neighbors=15, use_rep='X')
sc.tl.leiden(adata_temp, resolution=0.5)
adata.obs['leiden_ccc'] = adata_temp.obs['leiden'].values
print(f"      Identified {len(adata.obs['leiden_ccc'].unique())} clusters")

# Louvain
print("    Running Louvain...")
adata_temp = sc.AnnData(features)
adata_temp.obs_names = adata.obs_names
sc.pp.neighbors(adata_temp, n_neighbors=15, use_rep='X')
sc.tl.louvain(adata_temp, resolution=0.5)
adata.obs['louvain_ccc'] = adata_temp.obs['louvain'].values
print(f"      Identified {len(adata.obs['louvain_ccc'].unique())} clusters")

# K-means (try different k values)
print("    Running K-means...")
for k in [4, 6, 8]:
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    adata.obs[f'kmeans_{k}'] = kmeans.fit_predict(features).astype(str)
    print(f"      K={k}: Inertia={kmeans.inertia_:.2f}")

# Visualization
print(f"\n[7] Creating visualization...")
fig = plt.figure(figsize=(20, 12))

# Original Spatial distribution (If has Cell type labels)
if 'cell_type' in adata.obs.columns or 'celltype' in adata.obs.columns:
    ct_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'celltype'

    ax = plt.subplot(3, 4, 1)
    for ct in adata.obs[ct_col].unique():
        mask = adata.obs[ct_col] == ct
        ax.scatter(coords[mask, 0], coords[mask, 1], label=str(ct), s=10, alpha=0.7)
    ax.set_title('Original Cell Types', fontweight='bold')
    ax.legend(fontsize=7, ncol=2)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    plot_offset = 1
else:
    plot_offset = 0

# Leiden Results
ax = plt.subplot(3, 4, 1+plot_offset)
for label in sorted(adata.obs['leiden_ccc'].unique()):
    mask = adata.obs['leiden_ccc'] == label
    ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=10, alpha=0.7)
ax.set_title('Leiden Clustering (CCC)', fontweight='bold')
ax.legend(fontsize=7, ncol=2)
ax.set_xlabel('X')
ax.set_ylabel('Y')

# Louvain Results
ax = plt.subplot(3, 4, 2+plot_offset)
for label in sorted(adata.obs['louvain_ccc'].unique()):
    mask = adata.obs['louvain_ccc'] == label
    ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=10, alpha=0.7)
ax.set_title('Louvain Clustering (CCC)', fontweight='bold')
ax.legend(fontsize=7, ncol=2)
ax.set_xlabel('X')
ax.set_ylabel('Y')

# K-means Results
for i, k in enumerate([4, 6, 8]):
    ax = plt.subplot(3, 4, 3+plot_offset+i)
    for label in sorted(adata.obs[f'kmeans_{k}'].unique()):
        mask = adata.obs[f'kmeans_{k}'] == label
        ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=10, alpha=0.7)
    ax.set_title(f'K-means (K={k})', fontweight='bold')
    ax.legend(fontsize=7, ncol=2)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

# CCC Signal
ax = plt.subplot(3, 4, 7)
sender_total = sender.sum(axis=1)
scatter = ax.scatter(coords[:, 0], coords[:, 1],
                    c=sender_total, s=10, cmap='Reds', alpha=0.7)
ax.set_title('Total Sender Signal', fontweight='bold')
plt.colorbar(scatter, ax=ax)
ax.set_xlabel('X')
ax.set_ylabel('Y')

ax = plt.subplot(3, 4, 8)
receiver_total = receiver.sum(axis=1)
scatter = ax.scatter(coords[:, 0], coords[:, 1],
                    c=receiver_total, s=10, cmap='Blues', alpha=0.7)
ax.set_title('Total Receiver Signal', fontweight='bold')
plt.colorbar(scatter, ax=ax)
ax.set_xlabel('X')
ax.set_ylabel('Y')

# Clustering number Compare
ax = plt.subplot(3, 4, 9)
methods = ['leiden_ccc', 'louvain_ccc', 'kmeans_4', 'kmeans_6', 'kmeans_8']
n_clusters = [len(adata.obs[m].unique()) for m in methods]
ax.bar(range(len(methods)), n_clusters, alpha=0.8, color='steelblue')
ax.set_xticks(range(len(methods)))
ax.set_xticklabels(['Leiden', 'Louvain', 'K=4', 'K=6', 'K=8'], rotation=45, ha='right')
ax.set_ylabel('Number of Clusters')
ax.set_title('Detected Clusters', fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Summary
ax = plt.subplot(3, 4, 10)
ax.axis('off')

summary = "REAL DATA ANALYSIS\n" + "="*30 + "\n\n"
summary += f"Dataset:\n"
summary += f"  Cells: {adata.n_obs}\n"
summary += f"  Genes: {adata.n_vars}\n"
summary += f"  CCC Features: {features.shape[1]}\n\n"
summary += f"Clustering Results:\n"
summary += f"  Leiden: {len(adata.obs['leiden_ccc'].unique())} clusters\n"
summary += f"  Louvain: {len(adata.obs['louvain_ccc'].unique())} clusters\n"
summary += f"  K-means (K=4): 4 clusters\n"
summary += f"  K-means (K=6): 6 clusters\n"
summary += f"  K-means (K=8): 8 clusters\n"

ax.text(0.1, 0.9, summary, transform=ax.transAxes,
        fontsize=10, verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

plt.suptitle('Real Data - COMMOT Analysis', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('Real_Data_Analysis.png', dpi=300, bbox_inches='tight')
print("    Saved: Real_Data_Analysis.png")

# Save Results
adata.write_h5ad('Real_Data_with_CCC_clustering.h5ad')
print("    Saved: Real_Data_with_CCC_clustering.h5ad")

print("\n" + "="*80)
print("Real Data Analysis Complete!")
print("="*80)
print(f"\nSummary:")
print(f"  Dataset: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"  CCC Features: {features.shape[1]} dimensions")
print(f"  Leiden: {len(adata.obs['leiden_ccc'].unique())} clusters")
print(f"  Louvain: {len(adata.obs['louvain_ccc'].unique())} clusters")
print(f"\nNote: Without ground-truth labels, we cannot compute ARI/NMI")
print(f"      But we can see spatial patterns in the clustering results")
print("="*80)
