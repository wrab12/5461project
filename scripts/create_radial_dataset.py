"""
Create Radial Dataset
Design: Center-to-periphery structure with directional CCC patterns
- CellType_0: Core center
- CellType_1: Inner ring
- CellType_2: Middle ring
- CellType_3: Outer ring
"""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print(" "*20 + "Creating Radial Dataset")
print("="*80)

print("\n[1] Generating radial spatial pattern...")
print("    Design: Center-to-periphery structure")
print("    - CellType_0: Core (radius 0-15)")
print("    - CellType_1: Inner ring (radius 15-35)")
print("    - CellType_2: Middle ring (radius 35-55)")
print("    - CellType_3: Outer ring (radius 55-75)")

from optimized_synthetic_data import OptimizedSpatialCCCSimulator

class RadialSimulator(OptimizedSpatialCCCSimulator):
    def generate_spatial_coordinates(self):
        """Generate radial distribution"""
        n_per_type = self.n_cells // self.n_cell_types
        coords = np.zeros((self.n_cells, 2))
        cell_types = np.zeros(self.n_cells, dtype=int)

        center = np.array([75, 75])  # Center point

        # Define radius range for each cell type
        radius_ranges = [
            (0, 15),    # CellType_0: Core
            (15, 35),   # CellType_1: Inner ring
            (35, 55),   # CellType_2: Middle ring
            (55, 75),   # CellType_3: Outer ring
        ]

        for i in range(self.n_cell_types):
            start_idx = i * n_per_type
            end_idx = start_idx + n_per_type

            r_min, r_max = radius_ranges[i]

            # Uniform sampling in ring area
            n_cells_type = n_per_type

            # Use polar coordinates to generate
            if r_min == 0:  # Core area
                # Uniform distribution in disk
                r = np.sqrt(np.random.uniform(0, r_max**2, n_cells_type))
            else:  # Ring area
                r = np.sqrt(np.random.uniform(r_min**2, r_max**2, n_cells_type))

            theta = np.random.uniform(0, 2*np.pi, n_cells_type)

            # Convert to Cartesian coordinates
            x = center[0] + r * np.cos(theta)
            y = center[1] + r * np.sin(theta)

            coords[start_idx:end_idx, 0] = x
            coords[start_idx:end_idx, 1] = y
            cell_types[start_idx:end_idx] = i

        # Ensure within boundaries
        coords = np.clip(coords, [0, 0], self.spatial_dim)

        return coords, cell_types

# Create radial dataset
simulator = RadialSimulator(
    n_cells=1000,
    n_cell_types=4,
    spatial_dim=(150, 150),
    communication_radius=25.0,  # Sufficient to cover adjacent rings
    spatial_pattern='radial',
    expression_strength=15.0,
    heterogeneity=0.12,
    random_seed=42
)

adata = simulator.generate_dataset(
    expression_pattern='simple',
    add_noise=True,
    save_path='Radial.h5ad'
)

print(f"    Generated: {adata.n_obs} cells")

# Verify spatial distribution
print("\n[2] Verifying radial structure...")
coords = adata.obsm['spatial']
center = np.array([75, 75])

for ct in sorted(adata.obs['cell_type'].unique()):
    mask = adata.obs['cell_type'] == ct
    ct_coords = coords[mask]

    # Calculate distance to center
    distances = np.sqrt(((ct_coords - center)**2).sum(axis=1))

    print(f"    {ct}:")
    print(f"      Radius range: [{distances.min():.2f}, {distances.max():.2f}]")
    print(f"      Mean radius: {distances.mean():.2f}")
    print(f"      Cell count: {len(ct_coords)}")

print("\n[3] Running COMMOT analysis...")
import commot as ct

df_ligrec = pd.DataFrame({
    'ligand': ['Ligand_A', 'Ligand_B', 'Ligand_C', 'Ligand_D'],
    'receptor': ['Receptor_A', 'Receptor_B', 'Receptor_C', 'Receptor_D'],
    'pathway': ['Pathway_1', 'Pathway_2', 'Pathway_3', 'Pathway_4']
})

ct.tl.spatial_communication(
    adata,
    database_name='simulated',
    df_ligrec=df_ligrec,
    dis_thr=25,
    heteromeric=False,
    pathway_sum=False
)
print("    COMMOT SUCCESS!")

# Extracting features
sender = adata.obsm['commot-simulated-sum-sender'].values
receiver = adata.obsm['commot-simulated-sum-receiver'].values
features = np.hstack([sender, receiver])
print(f"    Features: {features.shape}")

# Clustering
print("\n[4] Clustering...")
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans

# Leiden
adata_temp = sc.AnnData(features)
adata_temp.obs_names = adata.obs_names
sc.pp.neighbors(adata_temp, n_neighbors=20, use_rep='X')
sc.tl.leiden(adata_temp, resolution=0.8)
adata.obs['leiden'] = adata_temp.obs['leiden'].values

# Louvain
adata_temp = sc.AnnData(features)
adata_temp.obs_names = adata.obs_names
sc.pp.neighbors(adata_temp, n_neighbors=20, use_rep='X')
sc.tl.louvain(adata_temp, resolution=0.8)
adata.obs['louvain'] = adata_temp.obs['louvain'].values

# K-means
kmeans = KMeans(n_clusters=4, random_state=42, n_init=10)
adata.obs['kmeans'] = kmeans.fit_predict(features).astype(str)

# Evaluation
print("\n[5] Evaluation:")
ground_truth = adata.obs['cell_type'].values
results = {}

print(f"    {'Method':<12} {'ARI':<10} {'NMI':<10} {'N_Clusters':<12}")
print("    " + "-"*44)

for method in ['leiden', 'louvain', 'kmeans']:
    pred = adata.obs[method].values
    ari = adjusted_rand_score(ground_truth, pred)
    nmi = normalized_mutual_info_score(ground_truth, pred)
    n_clusters = len(np.unique(pred))
    results[method] = {'ARI': ari, 'NMI': nmi, 'N_Clusters': n_clusters}
    print(f"    {method.capitalize():<12} {ari:<10.4f} {nmi:<10.4f} {n_clusters:<12}")

best_method = max(results.items(), key=lambda x: x[1]['ARI'])[0]
print(f"\n    Best: {best_method.capitalize()} (ARI={results[best_method]['ARI']:.4f})")

# Visualization
print("\n[6] Creating visualization...")
fig = plt.figure(figsize=(20, 10))

# Ground truth - polar coordinate style
ax = plt.subplot(2, 4, 1)
for ct in sorted(adata.obs['cell_type'].unique()):
    mask = adata.obs['cell_type'] == ct
    ax.scatter(coords[mask, 0], coords[mask, 1], label=ct, s=15, alpha=0.7)

# Draw ring boundaries
center_point = [75, 75]
for radius in [15, 35, 55, 75]:
    circle = plt.Circle(center_point, radius, fill=False, color='gray',
                       linestyle='--', linewidth=1.5, alpha=0.5)
    ax.add_patch(circle)

ax.scatter(center_point[0], center_point[1], s=200, marker='*',
          color='red', edgecolors='black', linewidths=2, zorder=10)
ax.set_title('Ground Truth\n(Radial Structure)', fontweight='bold')
ax.legend(fontsize=8)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_aspect('equal')
ax.grid(alpha=0.3)

# Clustering Results
titles = ['Leiden', 'Louvain', 'K-means']
columns = ['leiden', 'louvain', 'kmeans']

for idx, (title, col) in enumerate(zip(titles, columns)):
    ax = plt.subplot(2, 4, idx+2)
    labels = adata.obs[col]
    for label in sorted(labels.unique()):
        mask = labels == label
        ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=15, alpha=0.7)

    # Draw ring reference lines
    for radius in [15, 35, 55, 75]:
        circle = plt.Circle(center_point, radius, fill=False, color='gray',
                           linestyle='--', linewidth=1, alpha=0.3)
        ax.add_patch(circle)

    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, ncol=2)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    ax.grid(alpha=0.3)

# CCC Signal
ax = plt.subplot(2, 4, 5)
sender_total = sender.sum(axis=1)
scatter = ax.scatter(coords[:, 0], coords[:, 1], c=sender_total, s=15, cmap='Reds', alpha=0.7)
ax.set_title('Sender Signals', fontweight='bold')
plt.colorbar(scatter, ax=ax)
ax.set_aspect('equal')

ax = plt.subplot(2, 4, 6)
receiver_total = receiver.sum(axis=1)
scatter = ax.scatter(coords[:, 0], coords[:, 1], c=receiver_total, s=15, cmap='Blues', alpha=0.7)
ax.set_title('Receiver Signals', fontweight='bold')
plt.colorbar(scatter, ax=ax)
ax.set_aspect('equal')

# Performance Compare
ax = plt.subplot(2, 4, 7)
methods = list(results.keys())
ari_scores = [results[m]['ARI'] for m in methods]
nmi_scores = [results[m]['NMI'] for m in methods]
x = np.arange(len(methods))
width = 0.35

bars1 = ax.bar(x - width/2, ari_scores, width, label='ARI', alpha=0.8, color='steelblue')
bars2 = ax.bar(x + width/2, nmi_scores, width, label='NMI', alpha=0.8, color='coral')

for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.3f}', ha='center', va='bottom', fontsize=9)

ax.set_xlabel('Method')
ax.set_ylabel('Score')
ax.set_title('Performance', fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([m.capitalize() for m in methods])
ax.legend()
ax.set_ylim([0, 1])
ax.grid(axis='y', alpha=0.3)

# Design description
ax = plt.subplot(2, 4, 8)
ax.axis('off')

description = "RADIAL DESIGN\n" + "="*30 + "\n\n"
description += "Structure:\n"
description += "  Core (0-15)\n"
description += "  Inner ring (15-35)\n"
description += "  Middle ring (35-55)\n"
description += "  Outer ring (55-75)\n\n"
description += "CCC Pattern:\n"
description += "  Center -> Periphery\n"
description += "  Directional signaling\n\n"
description += f"Performance:\n"
description += f"  Best: {best_method.capitalize()}\n"
description += f"  ARI: {results[best_method]['ARI']:.3f}\n"
description += f"  NMI: {results[best_method]['NMI']:.3f}\n"

ax.text(0.1, 0.9, description, transform=ax.transAxes,
        fontsize=10, verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

plt.suptitle('Radial Dataset Analysis', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('Radial_Analysis.png', dpi=300, bbox_inches='tight')
print("    Saved: Radial_Analysis.png")

# Save Results
adata.write_h5ad('Radial.h5ad')
print("    Saved: Radial.h5ad")

print("\n" + "="*80)
print("RADIAL DATASET COMPLETE!")
print("="*80)
print(f"\nResults:")
print(f"  Best method: {best_method.capitalize()}")
print(f"  ARI: {results[best_method]['ARI']:.4f}")
print(f"  NMI: {results[best_method]['NMI']:.4f}")
print("\nRadial structure provides directional CCC patterns")
print("from center to periphery, ideal for testing CCC-based clustering.")
print("="*80)
