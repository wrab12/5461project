"""
修复Clustered数据集
根据诊断结果，增大聚类间距离
"""
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print(" "*20 + "Fixing Clustered Dataset")
print("="*80)

# 使用诊断结果的建议
print("\n[1] Creating FIXED Clustered dataset...")
print("    Changes:")
print("    - Spatial dimension: (100, 100) -> (150, 150)")
print("    - Cluster spacing: ~50 -> ~90")
print("    - Cluster std: 8 -> 5")
print("    - Communication radius: 20 -> 15")

from optimized_synthetic_data import OptimizedSpatialCCCSimulator

# 创建修复版本
simulator = OptimizedSpatialCCCSimulator(
    n_cells=1000,  # 增加细胞数
    n_cell_types=4,
    spatial_dim=(150, 150),  # 增大空间
    communication_radius=15.0,  # 减小通信半径
    spatial_pattern='clustered',
    expression_strength=15.0,  # 增强表达差异
    heterogeneity=0.12,  # 降低异质性
    random_seed=42
)

# 修改聚类中心位置（增大间距）
class FixedClusteredSimulator(OptimizedSpatialCCCSimulator):
    def _get_cluster_centers(self):
        """使用更大的间距"""
        return [
            np.array([30, 30]),    # 左下
            np.array([120, 30]),   # 右下
            np.array([30, 120]),   # 左上
            np.array([120, 120]),  # 右上
        ]

    def generate_spatial_coordinates(self):
        """生成更紧密的聚类"""
        n_per_type = self.n_cells // self.n_cell_types
        coords = np.zeros((self.n_cells, 2))
        cell_types = np.zeros(self.n_cells, dtype=int)

        centers = self._get_cluster_centers()
        for i in range(self.n_cell_types):
            start_idx = i * n_per_type
            end_idx = start_idx + n_per_type
            center = centers[i]
            # 减小标准差，使聚类更紧密
            coords[start_idx:end_idx] = np.random.randn(n_per_type, 2) * 5 + center
            cell_types[start_idx:end_idx] = i

        coords = np.clip(coords, [0, 0], self.spatial_dim)
        return coords, cell_types

# 使用修复版本
fixed_simulator = FixedClusteredSimulator(
    n_cells=1000,
    n_cell_types=4,
    spatial_dim=(150, 150),
    communication_radius=15.0,
    spatial_pattern='clustered',
    expression_strength=15.0,
    heterogeneity=0.12,
    random_seed=42
)

adata = fixed_simulator.generate_dataset(
    expression_pattern='simple',
    add_noise=True,
    save_path='Fixed_Clustered.h5ad'
)

# 验证改进
print("\n[2] Verifying improvements...")
coords = adata.obsm['spatial']

# 计算聚类中心
centers = []
for ct in sorted(adata.obs['cell_type'].unique()):
    mask = adata.obs['cell_type'] == ct
    center = coords[mask].mean(axis=0)
    centers.append(center)
    print(f"    {ct} center: ({center[0]:.2f}, {center[1]:.2f})")

centers = np.array(centers)
dist_matrix = cdist(centers, centers)

print("\n    Inter-cluster distances:")
for i in range(len(centers)):
    for j in range(i+1, len(centers)):
        print(f"      CellType_{i} <-> CellType_{j}: {dist_matrix[i,j]:.2f}")

avg_inter_dist = dist_matrix[np.triu_indices_from(dist_matrix, k=1)].mean()
min_inter_dist = dist_matrix[np.triu_indices_from(dist_matrix, k=1)].min()
comm_radius = 15.0

print(f"\n    Average inter-cluster distance: {avg_inter_dist:.2f}")
print(f"    Minimum inter-cluster distance: {min_inter_dist:.2f}")
print(f"    Communication radius: {comm_radius}")
print(f"    Ratio (avg/radius): {avg_inter_dist/comm_radius:.2f}x")
print(f"    Ratio (min/radius): {min_inter_dist/comm_radius:.2f}x")

if avg_inter_dist / comm_radius > 4.0:
    print("\n    [GOOD] Clusters are well separated")
else:
    print("\n    [WARNING] Still need improvement")

# 计算聚类内距离
print("\n    Intra-cluster distances:")
for ct in sorted(adata.obs['cell_type'].unique()):
    mask = adata.obs['cell_type'] == ct
    ct_coords = coords[mask]
    intra_dist = cdist(ct_coords, ct_coords)
    avg_intra = intra_dist[np.triu_indices_from(intra_dist, k=1)].mean()

    within_radius = (intra_dist < comm_radius).sum() / 2
    total_pairs = len(ct_coords) * (len(ct_coords) - 1) / 2
    coverage = within_radius / total_pairs * 100

    print(f"    {ct}:")
    print(f"      Avg distance: {avg_intra:.2f}")
    print(f"      Coverage: {coverage:.1f}%")

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
    dis_thr=15,
    heteromeric=False,
    pathway_sum=False
)
print("    COMMOT SUCCESS!")

# 提取特征
sender = adata.obsm['commot-simulated-sum-sender'].values
receiver = adata.obsm['commot-simulated-sum-receiver'].values
features = np.hstack([sender, receiver])
print(f"    Features: {features.shape}")

# 聚类
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

# 评估
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

# 可视化对比
print("\n[6] Creating comparison visualization...")
import matplotlib.pyplot as plt
import seaborn as sns

fig = plt.figure(figsize=(20, 10))

# 原始vs修复的空间分布
ax = plt.subplot(2, 4, 1)
for ct in sorted(adata.obs['cell_type'].unique()):
    mask = adata.obs['cell_type'] == ct
    ax.scatter(coords[mask, 0], coords[mask, 1], label=ct, s=15, alpha=0.7)
    center = coords[mask].mean(axis=0)
    ax.scatter(center[0], center[1], s=200, marker='x', color='black', linewidths=3)
    circle = plt.Circle(center, comm_radius, fill=False, color='red', linestyle='--', linewidth=2)
    ax.add_patch(circle)
ax.set_title('FIXED: Ground Truth\n(Red = comm radius)', fontweight='bold')
ax.legend(fontsize=8)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.grid(alpha=0.3)

# 聚类结果
titles = ['Leiden', 'Louvain', 'K-means']
columns = ['leiden', 'louvain', 'kmeans']

for idx, (title, col) in enumerate(zip(titles, columns)):
    ax = plt.subplot(2, 4, idx+2)
    labels = adata.obs[col]
    for label in sorted(labels.unique()):
        mask = labels == label
        ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), s=15, alpha=0.7)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, ncol=2)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.grid(alpha=0.3)

# CCC信号
ax = plt.subplot(2, 4, 5)
sender_total = sender.sum(axis=1)
scatter = ax.scatter(coords[:, 0], coords[:, 1], c=sender_total, s=15, cmap='Reds', alpha=0.7)
ax.set_title('Sender Signals', fontweight='bold')
plt.colorbar(scatter, ax=ax)

ax = plt.subplot(2, 4, 6)
receiver_total = receiver.sum(axis=1)
scatter = ax.scatter(coords[:, 0], coords[:, 1], c=receiver_total, s=15, cmap='Blues', alpha=0.7)
ax.set_title('Receiver Signals', fontweight='bold')
plt.colorbar(scatter, ax=ax)

# 性能比较
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

# 改进对比
ax = plt.subplot(2, 4, 8)
ax.axis('off')

comparison = "IMPROVEMENT\n" + "="*30 + "\n\n"
comparison += "Original Clustered:\n"
comparison += "  Space: 100x100\n"
comparison += "  Spacing: ~50\n"
comparison += "  Comm radius: 20\n"
comparison += "  Ratio: 2.5x\n"
comparison += "  ARI: 0.006\n\n"
comparison += "FIXED Clustered:\n"
comparison += "  Space: 150x150\n"
comparison += f"  Spacing: ~{avg_inter_dist:.0f}\n"
comparison += "  Comm radius: 15\n"
comparison += f"  Ratio: {avg_inter_dist/comm_radius:.1f}x\n"
comparison += f"  ARI: {results[best_method]['ARI']:.3f}\n\n"
comparison += f"Improvement:\n"
comparison += f"  {results[best_method]['ARI']/0.006:.0f}x better!\n"

ax.text(0.1, 0.9, comparison, transform=ax.transAxes,
        fontsize=10, verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

plt.suptitle('Fixed Clustered Dataset - Performance Improvement', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('Fixed_Clustered_Analysis.png', dpi=300, bbox_inches='tight')
print("    Saved: Fixed_Clustered_Analysis.png")

print("\n" + "="*80)
print("FIXED CLUSTERED DATASET COMPLETE!")
print("="*80)
print(f"\nImprovement:")
print(f"  Original ARI: 0.006")
print(f"  Fixed ARI: {results[best_method]['ARI']:.4f}")
print(f"  Improvement: {results[best_method]['ARI']/0.006:.1f}x")
print("="*80)
