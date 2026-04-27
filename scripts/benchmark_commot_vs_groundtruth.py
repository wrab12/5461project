"""
True COMMOT Benchmark: COMMOT Inference vs Ground-truth
Evaluate COMMOT method accuracy
"""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print(" "*15 + "COMMOT vs Ground-truth Benchmark")
print("="*80)

# Select datasets to test
datasets = [
    'Optimized_Tissue.h5ad',
    'Optimized_Striped.h5ad',
    'Optimized_Gradient.h5ad',
    'Radial.h5ad',
]

all_results = {}

for dataset_name in datasets:
    print(f"\n{'='*80}")
    print(f"Dataset: {dataset_name}")
    print(f"{'='*80}")

    try:
        adata = sc.read_h5ad(dataset_name)
    except:
        print(f"  [SKIP] File not found")
        continue

    # Check if has ground-truth
    if 'ground_truth_communication' not in adata.uns:
        print(f"  [SKIP] No ground-truth communication matrix")
        continue

    print(f"\n[1] Loading data...")
    print(f"    Cells: {adata.n_obs}")
    print(f"    Cell types: {adata.obs['cell_type'].nunique()}")

    coords = adata.obsm['spatial']
    ground_truth = adata.obs['cell_type'].values

    # Extract ground-truth communication matrix
    gt_comm = adata.uns['ground_truth_communication']
    print(f"    Ground-truth communication matrices: {len(gt_comm)} pairs")

    # ========================================================================
    # Method 1: Extract CCC features from Ground-truth communication matrix
    # ========================================================================
    print(f"\n[2] Extracting CCC features from GROUND-TRUTH...")

    # Calculate sender and receiver signals for each cell (based on ground-truth)
    n_cells = adata.n_obs
    n_pairs = len(gt_comm)

    gt_sender = np.zeros((n_cells, n_pairs))
    gt_receiver = np.zeros((n_cells, n_pairs))

    for idx, (pair_name, comm_matrix) in enumerate(gt_comm.items()):
        # Sender signal: total signal this cell sends to all other cells
        gt_sender[:, idx] = comm_matrix.sum(axis=1)
        # Receiver signal: total signal this cell receives from all other cells
        gt_receiver[:, idx] = comm_matrix.sum(axis=0)

    gt_features = np.hstack([gt_sender, gt_receiver])
    print(f"    Ground-truth features shape: {gt_features.shape}")
    print(f"    Feature range: [{gt_features.min():.4f}, {gt_features.max():.4f}]")

    # ========================================================================
    # Method 2: Extract CCC features from COMMOT inferred communication matrix
    # ========================================================================
    print(f"\n[3] Extracting CCC features from COMMOT...")

    # Check if already has COMMOT results
    commot_keys = [k for k in adata.obsm.keys() if 'commot' in k.lower() and 'sender' in k.lower()]

    if len(commot_keys) == 0:
        print("    Running COMMOT analysis...")
        import commot as ct

        df_ligrec = pd.DataFrame({
            'ligand': ['Ligand_A', 'Ligand_B', 'Ligand_C', 'Ligand_D'],
            'receptor': ['Receptor_A', 'Receptor_B', 'Receptor_C', 'Receptor_D'],
            'pathway': ['Pathway_1', 'Pathway_2', 'Pathway_3', 'Pathway_4']
        })

        comm_radius = adata.uns.get('simulation_params', {}).get('communication_radius', 20)

        ct.tl.spatial_communication(
            adata,
            database_name='simulated',
            df_ligrec=df_ligrec,
            dis_thr=comm_radius,
            heteromeric=False,
            pathway_sum=False
        )
        print("    COMMOT completed!")

    # Extract COMMOT features
    sender_key = [k for k in adata.obsm.keys() if 'commot' in k.lower() and 'sender' in k.lower()][0]
    receiver_key = [k for k in adata.obsm.keys() if 'commot' in k.lower() and 'receiver' in k.lower()][0]

    commot_sender = adata.obsm[sender_key]
    commot_receiver = adata.obsm[receiver_key]

    if hasattr(commot_sender, 'values'):
        commot_sender = commot_sender.values
    if hasattr(commot_receiver, 'values'):
        commot_receiver = commot_receiver.values

    commot_features = np.hstack([commot_sender, commot_receiver])
    print(f"    COMMOT features shape: {commot_features.shape}")
    print(f"    Feature range: [{commot_features.min():.4f}, {commot_features.max():.4f}]")

    # ========================================================================
    # Comparison 1: Feature similarity
    # ========================================================================
    print(f"\n[4] Comparing feature similarity...")

    # Calculate correlation for each feature dimension
    feature_correlations = []
    for i in range(min(gt_features.shape[1], commot_features.shape[1])):
        corr, _ = pearsonr(gt_features[:, i], commot_features[:, i])
        feature_correlations.append(corr)

    avg_corr = np.mean(feature_correlations)
    print(f"    Average feature correlation: {avg_corr:.4f}")
    print(f"    Min correlation: {min(feature_correlations):.4f}")
    print(f"    Max correlation: {max(feature_correlations):.4f}")

    # ========================================================================
    # Comparison 2: Clustering performance
    # ========================================================================
    print(f"\n[5] Clustering comparison...")

    results = {
        'ground_truth': {},
        'commot': {}
    }

    # Use Ground-truth features for clustering
    print(f"    [A] Using GROUND-TRUTH features:")
    for method_name, features in [('ground_truth', gt_features), ('commot', commot_features)]:

        # Leiden
        adata_temp = sc.AnnData(features)
        adata_temp.obs_names = adata.obs_names
        sc.pp.neighbors(adata_temp, n_neighbors=20, use_rep='X')
        sc.tl.leiden(adata_temp, resolution=0.8)
        pred_leiden = adata_temp.obs['leiden'].values

        # Louvain
        adata_temp = sc.AnnData(features)
        adata_temp.obs_names = adata.obs_names
        sc.pp.neighbors(adata_temp, n_neighbors=20, use_rep='X')
        sc.tl.louvain(adata_temp, resolution=0.8)
        pred_louvain = adata_temp.obs['louvain'].values

        # K-means
        n_types = len(np.unique(ground_truth))
        kmeans = KMeans(n_clusters=n_types, random_state=42, n_init=10)
        pred_kmeans = kmeans.fit_predict(features).astype(str)

        # Evaluation
        for method, pred in [('Leiden', pred_leiden), ('Louvain', pred_louvain), ('K-means', pred_kmeans)]:
            ari = adjusted_rand_score(ground_truth, pred)
            nmi = normalized_mutual_info_score(ground_truth, pred)
            results[method_name][method] = {'ARI': ari, 'NMI': nmi}

            if method_name == 'ground_truth':
                print(f"        {method:<10} ARI={ari:.4f}, NMI={nmi:.4f}")

    print(f"    [B] Using COMMOT features:")
    for method in ['Leiden', 'Louvain', 'K-means']:
        ari = results['commot'][method]['ARI']
        nmi = results['commot'][method]['NMI']
        print(f"        {method:<10} ARI={ari:.4f}, NMI={nmi:.4f}")

    # Calculate performance gap
    print(f"\n[6] Performance gap (COMMOT vs Ground-truth):")
    for method in ['Leiden', 'Louvain', 'K-means']:
        gt_ari = results['ground_truth'][method]['ARI']
        cm_ari = results['commot'][method]['ARI']
        gap = cm_ari - gt_ari
        retention = (cm_ari / gt_ari * 100) if gt_ari > 0 else 0
        print(f"    {method:<10} Gap={gap:+.4f}, Retention={retention:.1f}%")

    # Save results
    all_results[dataset_name] = {
        'feature_correlation': avg_corr,
        'ground_truth_performance': results['ground_truth'],
        'commot_performance': results['commot']
    }

# ========================================================================
# Comprehensive report
# ========================================================================
print(f"\n{'='*80}")
print(" "*20 + "COMPREHENSIVE REPORT")
print(f"{'='*80}")

if len(all_results) == 0:
    print("\nNo results to report.")
else:
    print(f"\n{'Dataset':<25} {'Feature Corr':<15} {'GT ARI':<12} {'COMMOT ARI':<12} {'Gap':<10}")
    print("-"*80)

    for dataset_name, result in all_results.items():
        # Use K-means as representative
        gt_ari = result['ground_truth_performance']['K-means']['ARI']
        cm_ari = result['commot_performance']['K-means']['ARI']
        gap = cm_ari - gt_ari
        corr = result['feature_correlation']

        short_name = dataset_name.replace('Optimized_', '').replace('.h5ad', '')
        print(f"{short_name:<25} {corr:<15.4f} {gt_ari:<12.4f} {cm_ari:<12.4f} {gap:+.4f}")

    # Calculate averages
    avg_corr = np.mean([r['feature_correlation'] for r in all_results.values()])
    avg_gt_ari = np.mean([r['ground_truth_performance']['K-means']['ARI'] for r in all_results.values()])
    avg_cm_ari = np.mean([r['commot_performance']['K-means']['ARI'] for r in all_results.values()])
    avg_gap = avg_cm_ari - avg_gt_ari

    print("-"*80)
    print(f"{'AVERAGE':<25} {avg_corr:<15.4f} {avg_gt_ari:<12.4f} {avg_cm_ari:<12.4f} {avg_gap:+.4f}")

    print(f"\n{'='*80}")
    print("CONCLUSIONS")
    print(f"{'='*80}")

    print(f"\n1. Feature Correlation: {avg_corr:.4f}")
    if avg_corr > 0.9:
        print("   [EXCELLENT] COMMOT features highly correlated with ground-truth")
    elif avg_corr > 0.7:
        print("   [GOOD] COMMOT features well correlated with ground-truth")
    elif avg_corr > 0.5:
        print("   [MODERATE] COMMOT features moderately correlated with ground-truth")
    else:
        print("   [POOR] COMMOT features poorly correlated with ground-truth")

    print(f"\n2. Performance Retention: {(avg_cm_ari/avg_gt_ari*100):.1f}%")
    if avg_cm_ari / avg_gt_ari > 0.95:
        print("   [EXCELLENT] COMMOT retains >95% of ground-truth performance")
    elif avg_cm_ari / avg_gt_ari > 0.85:
        print("   [GOOD] COMMOT retains >85% of ground-truth performance")
    elif avg_cm_ari / avg_gt_ari > 0.7:
        print("   [MODERATE] COMMOT retains >70% of ground-truth performance")
    else:
        print("   [POOR] COMMOT loses significant performance vs ground-truth")

    print(f"\n3. Average Performance Gap: {avg_gap:+.4f}")
    if abs(avg_gap) < 0.05:
        print("   [EXCELLENT] Negligible difference between COMMOT and ground-truth")
    elif abs(avg_gap) < 0.1:
        print("   [GOOD] Small difference between COMMOT and ground-truth")
    elif abs(avg_gap) < 0.2:
        print("   [MODERATE] Moderate difference between COMMOT and ground-truth")
    else:
        print("   [SIGNIFICANT] Large difference between COMMOT and ground-truth")

    print(f"\n{'='*80}")
    print("This benchmark demonstrates:")
    print("  1. How accurately COMMOT infers cell-cell communication")
    print("  2. Whether COMMOT-inferred features match ground-truth features")
    print("  3. How much clustering performance is lost when using COMMOT vs ground-truth")
    print(f"{'='*80}")

print("\nBenchmark complete!")
