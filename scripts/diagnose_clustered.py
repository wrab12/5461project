"""
Analyze Clustered Dataset - Why Performance is Poor
"""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cdist

# Loading data
print("="*80)
print("Analyzing why Clustered dataset performs poorly")
print("="*80)

adata = sc.read_h5ad('Optimized_Clustered.h5ad')

print(f"\n[1] Basic Info:")
print(f"  Cells: {adata.n_obs}")
print(f"  Genes: {adata.n_vars}")
print(f"  Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")

# Analysis Spatial distribution
print(f"\n[2] Spatial Distribution Analysis:")
coords = adata.obsm['spatial']

for ct in sorted(adata.obs['cell_type'].unique()):
    mask = adata.obs['cell_type'] == ct
    ct_coords = coords[mask]
    center = ct_coords.mean(axis=0)
    std = ct_coords.std(axis=0)
    print(f"  {ct}:")
    print(f"    Center: ({center[0]:.2f}, {center[1]:.2f})")
    print(f"    Std: ({std[0]:.2f}, {std[1]:.2f})")
    print(f"    Range X: [{ct_coords[:,0].min():.2f}, {ct_coords[:,0].max():.2f}]")
    print(f"    Range Y: [{ct_coords[:,1].min():.2f}, {ct_coords[:,1].max():.2f}]")

# Calculate inter-cluster Distance
print(f"\n[3] Inter-cluster Distances:")
centers = []
for ct in sorted(adata.obs['cell_type'].unique()):
    mask = adata.obs['cell_type'] == ct
    center = coords[mask].mean(axis=0)
    centers.append(center)

centers = np.array(centers)
dist_matrix = cdist(centers, centers)

print("  Distance matrix between cluster centers:")
for i in range(len(centers)):
    for j in range(i+1, len(centers)):
        print(f"    CellType_{i} <-> CellType_{j}: {dist_matrix[i,j]:.2f}")

# Analysis expression pattern
print(f"\n[4] Expression Pattern Analysis:")
ligands = ['Ligand_A', 'Ligand_B', 'Ligand_C', 'Ligand_D']
receptors = ['Receptor_A', 'Receptor_B', 'Receptor_C', 'Receptor_D']

for i, (lig, rec) in enumerate(zip(ligands, receptors)):
    if lig in adata.var_names and rec in adata.var_names:
        lig_idx = adata.var_names.tolist().index(lig)
        rec_idx = adata.var_names.tolist().index(rec)

        print(f"\n  {lig} -> {rec}:")
        for ct in sorted(adata.obs['cell_type'].unique()):
            mask = adata.obs['cell_type'] == ct
            lig_expr = adata.X[mask, lig_idx].mean()
            rec_expr = adata.X[mask, rec_idx].mean()
            print(f"    {ct}: Ligand={lig_expr:.2f}, Receptor={rec_expr:.2f}")

# Analysis communication Radius coverage
print(f"\n[5] Communication Radius Coverage:")
comm_radius = adata.uns['simulation_params']['communication_radius']
print(f"  Communication radius: {comm_radius}")

# Calculate average distance within each cell type
for ct in sorted(adata.obs['cell_type'].unique()):
    mask = adata.obs['cell_type'] == ct
    ct_coords = coords[mask]

    # Calculate intra-type distance
    intra_dist = cdist(ct_coords, ct_coords)
    avg_intra_dist = intra_dist[np.triu_indices_from(intra_dist, k=1)].mean()

    # Calculate what proportion of cell pairs are within communication radius
    within_radius = (intra_dist < comm_radius).sum() / 2  # Divide by 2 because symmetric
    total_pairs = len(ct_coords) * (len(ct_coords) - 1) / 2
    coverage = within_radius / total_pairs * 100

    print(f"  {ct}:")
    print(f"    Avg intra-cluster distance: {avg_intra_dist:.2f}")
    print(f"    Pairs within comm radius: {coverage:.1f}%")

# Analysis CCC Features
print(f"\n[6] CCC Feature Analysis:")
if 'commot-simulated-sum-sender' in adata.obsm.keys():
    sender = adata.obsm['commot-simulated-sum-sender'].values
    receiver = adata.obsm['commot-simulated-sum-receiver'].values

    print(f"  Sender signal range: [{sender.min():.2f}, {sender.max():.2f}]")
    print(f"  Receiver signal range: [{receiver.min():.2f}, {receiver.max():.2f}]")

    print(f"\n  Average CCC signals by cell type:")
    for ct in sorted(adata.obs['cell_type'].unique()):
        mask = adata.obs['cell_type'] == ct
        avg_sender = sender[mask].mean(axis=0)
        avg_receiver = receiver[mask].mean(axis=0)
        print(f"    {ct}:")
        print(f"      Sender: {avg_sender}")
        print(f"      Receiver: {avg_receiver}")
        print(f"      Total: {avg_sender.sum() + avg_receiver.sum():.2f}")

# Visualization Diagnosis
print(f"\n[7] Creating diagnostic visualization...")
fig = plt.figure(figsize=(20, 10))

# Spatial distribution
ax = plt.subplot(2, 4, 1)
for ct in sorted(adata.obs['cell_type'].unique()):
    mask = adata.obs['cell_type'] == ct
    ax.scatter(coords[mask, 0], coords[mask, 1], label=ct, s=20, alpha=0.7)
    # Draw center point
    center = coords[mask].mean(axis=0)
    ax.scatter(center[0], center[1], s=200, marker='x', color='black', linewidths=3)
    # Draw communication radius
    circle = plt.Circle(center, comm_radius, fill=False, color='red', linestyle='--', linewidth=2)
    ax.add_patch(circle)

ax.set_title('Spatial Distribution\n(Red circles = comm radius)', fontweight='bold')
ax.legend()
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.grid(alpha=0.3)

# Inter-cluster Distance heatmap
ax = plt.subplot(2, 4, 2)
sns.heatmap(dist_matrix, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
            xticklabels=[f'CT{i}' for i in range(4)],
            yticklabels=[f'CT{i}' for i in range(4)])
ax.set_title('Inter-cluster Distances', fontweight='bold')

# Expression heatmap
ax = plt.subplot(2, 4, 3)
expr_matrix = np.zeros((4, 8))  # 4 cell types x 8 genes (4 ligands + 4 receptors)
for i, ct in enumerate(sorted(adata.obs['cell_type'].unique())):
    mask = adata.obs['cell_type'] == ct
    for j, gene in enumerate(ligands + receptors):
        if gene in adata.var_names:
            gene_idx = adata.var_names.tolist().index(gene)
            expr_matrix[i, j] = adata.X[mask, gene_idx].mean()

sns.heatmap(expr_matrix, annot=True, fmt='.1f', cmap='viridis', ax=ax,
            xticklabels=ligands + receptors,
            yticklabels=[f'CT{i}' for i in range(4)])
ax.set_title('Expression Pattern', fontweight='bold')
plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

# CCC Signal distribution
if 'commot-simulated-sum-sender' in adata.obsm.keys():
    ax = plt.subplot(2, 4, 4)
    total_signal = sender.sum(axis=1) + receiver.sum(axis=1)
    scatter = ax.scatter(coords[:, 0], coords[:, 1], c=total_signal, s=20, cmap='viridis', alpha=0.7)
    ax.set_title('Total CCC Signal', fontweight='bold')
    plt.colorbar(scatter, ax=ax)

    # CCC Signal boxplot for each cell type
    ax = plt.subplot(2, 4, 5)
    data_for_box = []
    labels_for_box = []
    for ct in sorted(adata.obs['cell_type'].unique()):
        mask = adata.obs['cell_type'] == ct
        data_for_box.append(total_signal[mask])
        labels_for_box.append(ct)

    ax.boxplot(data_for_box, labels=labels_for_box)
    ax.set_title('CCC Signal Distribution', fontweight='bold')
    ax.set_ylabel('Total CCC Signal')
    ax.grid(axis='y', alpha=0.3)

    # Sender and Receiver Signal Comparison
    ax = plt.subplot(2, 4, 6)
    for ct in sorted(adata.obs['cell_type'].unique()):
        mask = adata.obs['cell_type'] == ct
        avg_sender = sender[mask].sum(axis=1).mean()
        avg_receiver = receiver[mask].sum(axis=1).mean()
        ax.scatter(avg_sender, avg_receiver, s=200, label=ct, alpha=0.7)

    ax.set_xlabel('Avg Sender Signal')
    ax.set_ylabel('Avg Receiver Signal')
    ax.set_title('Sender vs Receiver', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

# Problem Diagnosis summary
ax = plt.subplot(2, 4, 7)
ax.axis('off')

# Calculate key metrics
avg_inter_dist = dist_matrix[np.triu_indices_from(dist_matrix, k=1)].mean()
min_inter_dist = dist_matrix[np.triu_indices_from(dist_matrix, k=1)].min()

diagnosis = "DIAGNOSIS\n" + "="*40 + "\n\n"
diagnosis += "Potential Issues:\n\n"

# Problem 1: Inter-cluster Distance
if avg_inter_dist < comm_radius * 2:
    diagnosis += "1. CLUSTERS TOO CLOSE\n"
    diagnosis += f"   Avg inter-dist: {avg_inter_dist:.1f}\n"
    diagnosis += f"   Comm radius: {comm_radius}\n"
    diagnosis += f"   Ratio: {avg_inter_dist/comm_radius:.2f}x\n"
    diagnosis += "   -> Clusters overlap in\n"
    diagnosis += "      communication range!\n\n"

# Problem 2: Intra-cluster Distance
avg_intra_dists = []
for ct in sorted(adata.obs['cell_type'].unique()):
    mask = adata.obs['cell_type'] == ct
    ct_coords = coords[mask]
    intra_dist = cdist(ct_coords, ct_coords)
    avg_intra_dists.append(intra_dist[np.triu_indices_from(intra_dist, k=1)].mean())

avg_intra = np.mean(avg_intra_dists)
if avg_intra < comm_radius:
    diagnosis += "2. HIGH INTRA-CLUSTER\n"
    diagnosis += "   COMMUNICATION\n"
    diagnosis += f"   Avg intra-dist: {avg_intra:.1f}\n"
    diagnosis += f"   Comm radius: {comm_radius}\n"
    diagnosis += "   -> Most cells within\n"
    diagnosis += "      cluster communicate\n\n"

# Problem 3: CCC Signal variability
if 'commot-simulated-sum-sender' in adata.obsm.keys():
    ccc_by_type = []
    for ct in sorted(adata.obs['cell_type'].unique()):
        mask = adata.obs['cell_type'] == ct
        total = sender[mask].sum() + receiver[mask].sum()
        ccc_by_type.append(total / mask.sum())

    ccc_std = np.std(ccc_by_type)
    ccc_mean = np.mean(ccc_by_type)
    cv = ccc_std / ccc_mean if ccc_mean > 0 else 0

    if cv < 0.3:
        diagnosis += "3. LOW CCC VARIABILITY\n"
        diagnosis += f"   CV: {cv:.3f}\n"
        diagnosis += "   -> CCC signals too\n"
        diagnosis += "      similar across types\n\n"

diagnosis += "\nRECOMMENDATIONS:\n"
diagnosis += "- Increase cluster spacing\n"
diagnosis += "- Reduce cluster std\n"
diagnosis += "- Increase expr strength\n"
diagnosis += "- Or reduce comm radius\n"

ax.text(0.05, 0.95, diagnosis, transform=ax.transAxes,
        fontsize=9, verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.suptitle('Clustered Dataset - Diagnostic Analysis', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('Clustered_Diagnosis.png', dpi=300, bbox_inches='tight')
print("  Saved: Clustered_Diagnosis.png")

print("\n" + "="*80)
print("CONCLUSION:")
print("="*80)
print(f"Average inter-cluster distance: {avg_inter_dist:.2f}")
print(f"Communication radius: {comm_radius}")
print(f"Ratio: {avg_inter_dist/comm_radius:.2f}x")
print(f"\nProblem: Clusters are too close relative to communication radius!")
print(f"This causes:")
print(f"  1. Overlapping communication ranges")
print(f"  2. Similar CCC patterns across cell types")
print(f"  3. Poor clustering performance")
print(f"\nSolution: Increase cluster spacing or reduce communication radius")
print("="*80)
