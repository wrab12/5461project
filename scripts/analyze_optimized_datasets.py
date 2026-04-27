"""
Analyze Optimized Datasets
"""
import numpy as np
import pandas as pd
import scanpy as sc
import commot as ct
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, confusion_matrix
from sklearn.cluster import KMeans
import warnings
warnings.filterwarnings('ignore')

def analyze_optimized_dataset(h5ad_path, dataset_name):
    """Analyze single optimized dataset"""
    print("\n" + "="*80)
    print(f"Analyzing: {dataset_name}")
    print("="*80)

    # Loading data
    print("[1] Loading data...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")
    print(f"  Spatial pattern: {adata.uns['simulation_params']['spatial_pattern']}")
    print(f"  Cell type distribution: {adata.obs['cell_type'].value_counts().to_dict()}")

    # Running COMMOT
    print("\n[2] Running COMMOT analysis...")

    # Creating Ligand-Receptor DataFrame
    df_ligrec = pd.DataFrame({
        'ligand': ['Ligand_A', 'Ligand_B', 'Ligand_C', 'Ligand_D'],
        'receptor': ['Receptor_A', 'Receptor_B', 'Receptor_C', 'Receptor_D'],
        'pathway': ['Pathway_1', 'Pathway_2', 'Pathway_3', 'Pathway_4']
    })

    # Based on spatial pattern adjust communication radius
    dis_thr = adata.uns['simulation_params']['communication_radius']
    print(f"  Communication radius: {dis_thr}")

    try:
        ct.tl.spatial_communication(
            adata,
            database_name='simulated',
            df_ligrec=df_ligrec,
            dis_thr=dis_thr,
            heteromeric=False,
            pathway_sum=False
        )
        print("  COMMOT analysis SUCCESS!")
    except Exception as e:
        print(f"  COMMOT failed: {e}")
        return None, None

    # Extract CCC Features
    print("\n[3] Extracting CCC features...")
    commot_keys = [k for k in adata.obsm.keys() if 'commot' in k]
    sender_key = [k for k in commot_keys if 'sender' in k][0]
    receiver_key = [k for k in commot_keys if 'receiver' in k][0]

    sender = adata.obsm[sender_key].values
    receiver = adata.obsm[receiver_key].values
    features = np.hstack([sender, receiver])
    print(f"  Feature dimension: {features.shape}")
    print(f"  Sender signal range: [{sender.min():.2f}, {sender.max():.2f}]")
    print(f"  Receiver signal range: [{receiver.min():.2f}, {receiver.max():.2f}]")

    # Clustering
    print("\n[4] Running clustering algorithms...")

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
    print("\n[5] Evaluating clustering performance...")
    ground_truth = adata.obs['cell_type'].values

    results = {}
    print(f"\n  {'Method':<12} {'ARI':<10} {'NMI':<10} {'N_Clusters':<12}")
    print("  " + "-"*44)

    for method in ['leiden', 'louvain', 'kmeans']:
        pred = adata.obs[method].values
        ari = adjusted_rand_score(ground_truth, pred)
        nmi = normalized_mutual_info_score(ground_truth, pred)
        n_clusters = len(np.unique(pred))
        results[method] = {'ARI': ari, 'NMI': nmi, 'N_Clusters': n_clusters}
        print(f"  {method.capitalize():<12} {ari:<10.4f} {nmi:<10.4f} {n_clusters:<12}")

    best_method = max(results.items(), key=lambda x: x[1]['ARI'])[0]
    print(f"\n  Best method: {best_method.capitalize()} (ARI={results[best_method]['ARI']:.4f})")

    # Visualization
    print("\n[6] Creating visualization...")
    fig = plt.figure(figsize=(20, 10))
    coords = adata.obsm['spatial']

    # Spatial distribution
    titles = ['Ground Truth', 'Leiden', 'Louvain', 'K-means']
    columns = ['cell_type', 'leiden', 'louvain', 'kmeans']

    for idx, (title, col) in enumerate(zip(titles, columns)):
        ax = plt.subplot(2, 4, idx+1)
        labels = adata.obs[col]
        unique_labels = sorted(labels.unique())

        for i, label in enumerate(unique_labels):
            mask = labels == label
            ax.scatter(coords[mask, 0], coords[mask, 1],
                      label=str(label), s=15, alpha=0.7)

        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.legend(fontsize=7, loc='best', ncol=2)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.grid(alpha=0.2)

    # CCC Signal
    ax = plt.subplot(2, 4, 5)
    sender_total = sender.sum(axis=1)
    scatter = ax.scatter(coords[:, 0], coords[:, 1],
                        c=sender_total, s=15, cmap='Reds', alpha=0.7)
    ax.set_title('Sender Signals', fontsize=12, fontweight='bold')
    plt.colorbar(scatter, ax=ax)

    ax = plt.subplot(2, 4, 6)
    receiver_total = receiver.sum(axis=1)
    scatter = ax.scatter(coords[:, 0], coords[:, 1],
                        c=receiver_total, s=15, cmap='Blues', alpha=0.7)
    ax.set_title('Receiver Signals', fontsize=12, fontweight='bold')
    plt.colorbar(scatter, ax=ax)

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

    # Confusion matrix
    ax = plt.subplot(2, 4, 8)
    cm = confusion_matrix(adata.obs['cell_type'], adata.obs[best_method])
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax, cbar=False)
    ax.set_title(f'Confusion Matrix\n({best_method.capitalize()})', fontsize=12, fontweight='bold')
    ax.set_xlabel('Predicted')
    ax.set_ylabel('True')

    plt.suptitle(f'{dataset_name} - Optimized Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_name = f'{dataset_name}_analysis.png'
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_name}")
    plt.close()

    return results, adata


def main():
    """Main function"""
    print("="*80)
    print(" "*20 + "Optimized Dataset Analysis")
    print("="*80)

    datasets = [
        ('Optimized_Clustered.h5ad', 'Optimized_Clustered'),
        ('Optimized_Tissue.h5ad', 'Optimized_Tissue'),
        ('Optimized_Striped.h5ad', 'Optimized_Striped'),
        ('Optimized_Gradient.h5ad', 'Optimized_Gradient'),
        ('Optimized_Mixed.h5ad', 'Optimized_Mixed'),
    ]

    all_results = {}

    for h5ad_path, dataset_name in datasets:
        try:
            results, adata = analyze_optimized_dataset(h5ad_path, dataset_name)
            if results is not None:
                all_results[dataset_name] = results
        except Exception as e:
            print(f"\n[ERROR] {dataset_name} analysis failed: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Comprehensive report
    print("\n" + "="*80)
    print(" "*25 + "COMPREHENSIVE REPORT")
    print("="*80)

    print(f"\n{'Dataset':<25} {'Best Method':<15} {'ARI':<10} {'NMI':<10}")
    print("-"*60)

    for dataset_name, results in all_results.items():
        best_method = max(results.items(), key=lambda x: x[1]['ARI'])[0]
        ari = results[best_method]['ARI']
        nmi = results[best_method]['NMI']
        print(f"{dataset_name:<25} {best_method.capitalize():<15} {ari:<10.4f} {nmi:<10.4f}")

    # Creating comprehensive comparison chart
    print("\nCreating comprehensive comparison...")
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    dataset_names = list(all_results.keys())
    methods = ['leiden', 'louvain', 'kmeans']

    # ARI Compare
    ax = axes[0, 0]
    for method in methods:
        ari_scores = [all_results[ds][method]['ARI'] for ds in dataset_names]
        ax.plot(range(len(dataset_names)), ari_scores, marker='o', label=method.capitalize(), linewidth=2)
    ax.set_xticks(range(len(dataset_names)))
    ax.set_xticklabels([ds.replace('Optimized_', '') for ds in dataset_names],
                       fontsize=9, rotation=45, ha='right')
    ax.set_ylabel('ARI Score', fontsize=11)
    ax.set_title('ARI Comparison Across Datasets', fontweight='bold', fontsize=12)
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_ylim([0, 1])

    # NMI Compare
    ax = axes[0, 1]
    for method in methods:
        nmi_scores = [all_results[ds][method]['NMI'] for ds in dataset_names]
        ax.plot(range(len(dataset_names)), nmi_scores, marker='s', label=method.capitalize(), linewidth=2)
    ax.set_xticks(range(len(dataset_names)))
    ax.set_xticklabels([ds.replace('Optimized_', '') for ds in dataset_names],
                       fontsize=9, rotation=45, ha='right')
    ax.set_ylabel('NMI Score', fontsize=11)
    ax.set_title('NMI Comparison Across Datasets', fontweight='bold', fontsize=12)
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_ylim([0, 1])

    # Average Performance
    ax = axes[0, 2]
    avg_ari = {method: np.mean([all_results[ds][method]['ARI'] for ds in dataset_names])
               for method in methods}
    avg_nmi = {method: np.mean([all_results[ds][method]['NMI'] for ds in dataset_names])
               for method in methods}

    x = np.arange(len(methods))
    width = 0.35
    bars1 = ax.bar(x - width/2, [avg_ari[m] for m in methods], width, label='Avg ARI', alpha=0.8, color='steelblue')
    bars2 = ax.bar(x + width/2, [avg_nmi[m] for m in methods], width, label='Avg NMI', alpha=0.8, color='coral')

    # Add value labels
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

    # Best Method frequency
    ax = axes[1, 0]
    best_methods = [max(all_results[ds].items(), key=lambda x: x[1]['ARI'])[0]
                    for ds in dataset_names]
    method_counts = {m: best_methods.count(m) for m in methods}
    colors = ['steelblue', 'coral', 'lightgreen']
    bars = ax.bar(method_counts.keys(), method_counts.values(), alpha=0.8, color=colors)
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}', ha='center', va='bottom', fontsize=10)
    ax.set_xlabel('Method', fontsize=11)
    ax.set_ylabel('Number of Datasets', fontsize=11)
    ax.set_title('Best Method Frequency', fontweight='bold', fontsize=12)
    ax.grid(axis='y', alpha=0.3)

    # Clustering number Compare
    ax = axes[1, 1]
    for method in methods:
        n_clusters = [all_results[ds][method]['N_Clusters'] for ds in dataset_names]
        ax.plot(range(len(dataset_names)), n_clusters, marker='o', label=method.capitalize(), linewidth=2)
    ax.axhline(y=4, color='red', linestyle='--', label='True (4)', linewidth=2)
    ax.set_xticks(range(len(dataset_names)))
    ax.set_xticklabels([ds.replace('Optimized_', '') for ds in dataset_names],
                       fontsize=9, rotation=45, ha='right')
    ax.set_ylabel('Number of Clusters', fontsize=11)
    ax.set_title('Detected Clusters', fontweight='bold', fontsize=12)
    ax.legend()
    ax.grid(alpha=0.3)

    # Summary statistics
    ax = axes[1, 2]
    ax.axis('off')

    summary = "OPTIMIZATION RESULTS\n" + "="*35 + "\n\n"
    summary += f"Total Datasets: {len(dataset_names)}\n\n"
    summary += "Average Performance:\n"
    for method in methods:
        summary += f"  {method.capitalize()}:\n"
        summary += f"    ARI: {avg_ari[method]:.4f}\n"
        summary += f"    NMI: {avg_nmi[method]:.4f}\n"

    summary += f"\nBest Method Overall:\n"
    overall_best = max(methods, key=lambda m: avg_ari[m])
    summary += f"  {overall_best.capitalize()}\n"
    summary += f"  Avg ARI: {avg_ari[overall_best]:.4f}\n"
    summary += f"  Avg NMI: {avg_nmi[overall_best]:.4f}\n"

    ax.text(0.1, 0.9, summary, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

    plt.suptitle('Optimized Datasets - Comprehensive Performance', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('Optimized_Comprehensive_Comparison.png', dpi=300, bbox_inches='tight')
    print("  Saved: Optimized_Comprehensive_Comparison.png")

    print("\n" + "="*80)
    print(" "*30 + "ANALYSIS COMPLETE!")
    print("="*80)

    return all_results


if __name__ == '__main__':
    results = main()
