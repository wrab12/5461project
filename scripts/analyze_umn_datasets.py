"""
Analyze UMN Letter-shaped Datasets
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

def analyze_letter_dataset(h5ad_path, letter_name):
    """Analyze single letter dataset"""
    print("\n" + "="*80)
    print(f"Analyzing: {letter_name}-shaped dataset")
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

    plt.suptitle(f'{letter_name}-shaped Dataset Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'Letter_{letter_name}_analysis.png', dpi=300, bbox_inches='tight')
    print(f"    Saved: Letter_{letter_name}_analysis.png")
    plt.close()

    return results


def main():
    print("="*80)
    print(" "*20 + "UMN Letter Datasets Analysis")
    print("="*80)

    letters = ['U', 'M', 'N', 'UMN']
    all_results = {}

    for letter in letters:
        try:
            results = analyze_letter_dataset(f'Letter_{letter}.h5ad', letter)
            all_results[letter] = results
        except Exception as e:
            print(f"\n[ERROR] {letter} analysis failed: {e}")
            import traceback
            traceback.print_exc()

    # Comprehensive report
    print("\n" + "="*80)
    print(" "*25 + "COMPREHENSIVE REPORT")
    print("="*80)

    print(f"\n{'Letter':<10} {'Best Method':<15} {'ARI':<10} {'NMI':<10}")
    print("-"*45)

    for letter, results in all_results.items():
        best_method = max(results.items(), key=lambda x: x[1]['ARI'])[0]
        ari = results[best_method]['ARI']
        nmi = results[best_method]['NMI']
        print(f"{letter:<10} {best_method.capitalize():<15} {ari:<10.4f} {nmi:<10.4f}")

    # Creating comprehensive comparison chart
    print("\nCreating comprehensive comparison...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    letter_names = list(all_results.keys())
    methods = ['leiden', 'louvain', 'kmeans']

    # ARI Compare
    ax = axes[0, 0]
    for method in methods:
        ari_scores = [all_results[l][method]['ARI'] for l in letter_names]
        ax.plot(range(len(letter_names)), ari_scores, marker='o', label=method.capitalize(), linewidth=2)
    ax.set_xticks(range(len(letter_names)))
    ax.set_xticklabels(letter_names)
    ax.set_ylabel('ARI Score', fontsize=11)
    ax.set_title('ARI Comparison', fontweight='bold', fontsize=12)
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_ylim([0, 1])

    # NMI Compare
    ax = axes[0, 1]
    for method in methods:
        nmi_scores = [all_results[l][method]['NMI'] for l in letter_names]
        ax.plot(range(len(letter_names)), nmi_scores, marker='s', label=method.capitalize(), linewidth=2)
    ax.set_xticks(range(len(letter_names)))
    ax.set_xticklabels(letter_names)
    ax.set_ylabel('NMI Score', fontsize=11)
    ax.set_title('NMI Comparison', fontweight='bold', fontsize=12)
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_ylim([0, 1])

    # Average Performance
    ax = axes[1, 0]
    avg_ari = {method: np.mean([all_results[l][method]['ARI'] for l in letter_names])
               for method in methods}
    avg_nmi = {method: np.mean([all_results[l][method]['NMI'] for l in letter_names])
               for method in methods}

    x = np.arange(len(methods))
    width = 0.35
    bars1 = ax.bar(x - width/2, [avg_ari[m] for m in methods], width, label='Avg ARI', alpha=0.8, color='steelblue')
    bars2 = ax.bar(x + width/2, [avg_nmi[m] for m in methods], width, label='Avg NMI', alpha=0.8, color='coral')

    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.3f}', ha='center', va='bottom', fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([m.capitalize() for m in methods])
    ax.set_ylabel('Score', fontsize=11)
    ax.set_title('Average Performance', fontweight='bold', fontsize=12)
    ax.legend()
    ax.set_ylim([0, 1])
    ax.grid(axis='y', alpha=0.3)

    # Summary
    ax = axes[1, 1]
    ax.axis('off')

    summary = "UMN LETTER DATASETS\n" + "="*30 + "\n\n"
    summary += f"Total Datasets: {len(letter_names)}\n\n"
    summary += "Average Performance:\n"
    for method in methods:
        summary += f"  {method.capitalize()}:\n"
        summary += f"    ARI: {avg_ari[method]:.4f}\n"
        summary += f"    NMI: {avg_nmi[method]:.4f}\n"

    summary += f"\nBest Overall:\n"
    overall_best = max(methods, key=lambda m: avg_ari[m])
    summary += f"  {overall_best.capitalize()}\n"
    summary += f"  Avg ARI: {avg_ari[overall_best]:.4f}\n"

    ax.text(0.1, 0.9, summary, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    plt.suptitle('UMN Letter Datasets - Comprehensive Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('UMN_Comprehensive_Analysis.png', dpi=300, bbox_inches='tight')
    print("  Saved: UMN_Comprehensive_Analysis.png")

    print("\n" + "="*80)
    print("UMN Letter Datasets Analysis Complete!")
    print("="*80)


if __name__ == '__main__':
    main()
