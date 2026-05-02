"""
Analyze Tuned UMN Letter-shaped Datasets
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

def analyze_tuned_dataset(h5ad_path, level):
    """Analyze single letter dataset"""
    print("\n" + "="*80)
    print(f"Analyzing: {level} dataset")
    print("="*80)

    # Loading data
    adata = sc.read_h5ad(h5ad_path)
    print(f"[1] Data: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"    Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")

    # Running COMMOT
    print("\n[2] Running COMMOT...")
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

    # Extract CCC Features
    sender = adata.obsm['commot-simulated-sum-sender'].values
    receiver = adata.obsm['commot-simulated-sum-receiver'].values
    features = np.hstack([sender, receiver])
    print(f"    Features: {features.shape}")

    # Clustering
    print("\n[3] Clustering...")

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
    n_true_types = len(adata.obs['cell_type'].unique())
    kmeans = KMeans(n_clusters=n_true_types, random_state=42, n_init=10)
    adata.obs['kmeans'] = kmeans.fit_predict(features).astype(str)

    # Evaluation
    print("\n[4] Evaluation:")
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
    print("\n[5] Visualization...")
    fig = plt.figure(figsize=(20, 10))
    coords = adata.obsm['spatial']

    # Spatial distribution
    titles = ['Ground Truth', 'Leiden', 'Louvain', 'K-means']
    columns = ['cell_type', 'leiden', 'louvain', 'kmeans']

    for idx, (title, col) in enumerate(zip(titles, columns)):
        ax = plt.subplot(2, 4, idx+1)
        labels = adata.obs[col]

        for label in sorted(labels.unique()):
            mask = labels == label
            ax.scatter(coords[mask, 0], coords[mask, 1],
                      label=str(label), s=15, alpha=0.7)

        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, ncol=2)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_aspect('equal')
        ax.grid(alpha=0.3)

    # CCC Signal
    ax = plt.subplot(2, 4, 5)
    sender_total = sender.sum(axis=1)
    scatter = ax.scatter(coords[:, 0], coords[:, 1],
                        c=sender_total, s=15, cmap='Reds', alpha=0.7)
    ax.set_title('Sender Signals', fontsize=12, fontweight='bold')
    plt.colorbar(scatter, ax=ax)
    ax.set_aspect('equal')

    ax = plt.subplot(2, 4, 6)
    receiver_total = receiver.sum(axis=1)
    scatter = ax.scatter(coords[:, 0], coords[:, 1],
                        c=receiver_total, s=15, cmap='Blues', alpha=0.7)
    ax.set_title('Receiver Signals', fontsize=12, fontweight='bold')
    plt.colorbar(scatter, ax=ax)
    ax.set_aspect('equal')

    # Performance Compare
    ax = plt.subplot(2, 4, 7)
    methods = list(results.keys())
    ari_scores = [results[m]['ARI'] for m in methods]
    nmi_scores = [results[m]['NMI'] for m in methods]
    x = np.arange(len(methods))
    width = 0.35

    ax.bar(x - width/2, ari_scores, width, label='ARI', alpha=0.8, color='steelblue')
    ax.bar(x + width/2, nmi_scores, width, label='NMI', alpha=0.8, color='coral')
    ax.set_xlabel('Method')
    ax.set_ylabel('Score')
    ax.set_title('Performance', fontsize=12, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([m.capitalize() for m in methods])
    ax.legend()
    ax.set_ylim([0, 1])
    ax.grid(axis='y', alpha=0.3)

    plt.suptitle(f'{level} Dataset Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'UMN_{level}_analysis.png', dpi=300, bbox_inches='tight')
    print(f"    Saved: UMN_{level}_analysis.png")
    plt.close()

    return results


def main():
    print("="*80)
    print(" "*20 + "Tuned UMN Letter Datasets Analysis")
    print("="*80)

    levels = ['easy', 'medium', 'hard']
    all_results = {}

    for level in levels:
        try:
            results = analyze_tuned_dataset(f'umn_data/{level}.h5ad', level)
            all_results[level] = results
        except Exception as e:
            print(f"\n[ERROR] {level} analysis failed: {e}")
            import traceback
            traceback.print_exc()

    # Comprehensive report
    print("\n" + "="*80)
    print(" "*25 + "COMPREHENSIVE REPORT")
    print("="*80)

    print(f"\n{'Letter':<10} {'Best Method':<15} {'ARI':<10} {'NMI':<10}")
    print("-"*45)

    for level, results in all_results.items():
        best_method = max(results.items(), key=lambda x: x[1]['ARI'])[0]
        ari = results[best_method]['ARI']
        nmi = results[best_method]['NMI']
        print(f"{letter:<10} {best_method.capitalize():<15} {ari:<10.4f} {nmi:<10.4f}")

if __name__ == '__main__':
    main()
