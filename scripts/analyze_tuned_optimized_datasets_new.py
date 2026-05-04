"""
Analyze Tuned Optimized Datasets
"""
import numpy as np
import pandas as pd
import scanpy as sc
import commot as ct
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, confusion_matrix, silhouette_score
from sklearn.cluster import KMeans
import warnings
warnings.filterwarnings('ignore')

def analyze_tuned_dataset(h5ad_path, dataset_name):
    """Analyze single tuned dataset"""
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
    print(f"\n  {'Method':<12} {'ARI':<10} {'NMI':<10} {'Silhouette':<12} {'N_Clusters':<12}")
    print("  " + "-"*58)
    
    for method in ['leiden', 'louvain', 'kmeans']:
        pred = adata.obs[method].values
        ari = adjusted_rand_score(ground_truth, pred)
        nmi = normalized_mutual_info_score(ground_truth, pred)
        n_clusters = len(np.unique(pred))
        
        if n_clusters > 1 and n_clusters < features.shape[0]:
            sil = silhouette_score(features, pred)
            
        else:
            sil = np.nan
            
        results[method] = {
                'ARI': ari,
                'NMI': nmi,
                'Silhouette': sil,
                'N_Clusters': n_clusters
                }
        print(f"  {method.capitalize():<12} {ari:<10.4f} {nmi:<10.4f} {sil:<12.4f} {n_clusters:<12}")

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

    plt.suptitle(f'{dataset_name} - Tuned Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_name = f'{dataset_name}_analysis.png'
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_name}")
    plt.close()

    return results, adata


def main():
    """Main function"""
    print("="*80)
    print(" "*20 + "Tuned Dataset Analysis")
    print("="*80)

    datasets = [
        ('optimized_synthetic_data/easy.h5ad', 'Optimized_Easy'),
        ('optimized_synthetic_data/medium.h5ad', 'Optimized_Medium'),
        ('optimized_synthetic_data/hard.h5ad', 'Optimized_Hard'),
    ]

    all_results = {}

    for h5ad_path, dataset_name in datasets:
        try:
            results, adata = analyze_tuned_dataset(h5ad_path, dataset_name)
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

    print(f"\n{'Dataset':<25} {'Best Method':<15} {'ARI':<10} {'NMI':<10} {'Silhouette':<12}")
    print("-"*75)
    
    
    for dataset_name, results in all_results.items():
        best_method = max(results.items(), key=lambda x: x[1]['ARI'])[0]
        ari = results[best_method]['ARI']
        nmi = results[best_method]['NMI']
        sil = results[best_method]['Silhouette']
        print(f"{dataset_name:<25} {best_method.capitalize():<15} {ari:<10.4f} {nmi:<10.4f} {sil:<12.4f}")
    
    # Save full evaluation results
    rows = []
    
    for dataset_name, results in all_results.items():
        for method, metrics in results.items():
            rows.append({
                "Dataset": dataset_name,
                "Method": method,
                "ARI": metrics["ARI"],
                "NMI": metrics["NMI"],
                "Silhouette": metrics["Silhouette"],
                "N_Clusters": metrics["N_Clusters"]
            })

    results_df = pd.DataFrame(rows)
    results_df.to_csv("tuned_optimized_clustering_evaluation_results.csv", index=False)

    print("\nSaved evaluation results to: tuned_optimized_clustering_evaluation_results.csv")

    return all_results


if __name__ == '__main__':
    results = main()
