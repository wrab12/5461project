# Spatial Transcriptomics Cell-Cell Communication Simulated Dataset Design Documentation

## 1. Project Background and Objectives

### 1.1 Project Motivation
According to the CSCI 5461 project proposal, our goal is to create a standardized simulated dataset for the COMMOT (COMMunication analysis by Optimal Transport) method, used for:
- Benchmark testing of the COMMOT method
- Evaluating clustering algorithm performance based on cell-cell communication (CCC) features
- Providing a testing platform for future CCC algorithm research

### 1.2 Core Objectives
1. Generate spatial transcriptomics data with ground-truth labels
2. Simulate realistic ligand-receptor interaction patterns
3. Input data into COMMOT to obtain communication features
4. Perform cell clustering using communication features
5. Evaluate consistency between clustering results and ground-truth

## 2. Dataset Design Principles

### 2.1 Spatial Structure Design

**Why this design?**
- **Spatial aggregation**: In real tissues, cells of the same type tend to cluster together to form functional regions
- **Implementation method**: Define a center point for each cell type, generate cell coordinates around the center using Gaussian distribution
- **Parameter selection**:
  - 4 cell types, evenly distributed on a circle
  - Standard deviation of 8 units for each cluster, ensuring some overlap but mainly aggregated

```python
# Centers distributed in a circle
angle = 2 * π * i / n_cell_types
center = (radius * cos(angle), radius * sin(angle))

# Generate cells around center
coords = N(center, σ=8)
```

**Advantages**:
1. Simulates spatial organization patterns of real tissues
2. Creates spatial proximity relationships between cell types
3. Facilitates verification of whether COMMOT can correctly identify spatial communication patterns

### 2.2 Ligand-Receptor Pair Design

**Designed 4 ligand-receptor pairs, divided into two patterns:**

#### Pattern 1: One-to-One Signal Transmission
- `Ligand_A → Receptor_A` (Pathway_1)
- `Ligand_B → Receptor_B` (Pathway_2)

**Why need one-to-one pattern?**
- This is the simplest, most direct communication pattern
- Facilitates verification of COMMOT's basic functionality
- Simulates specific signaling pathways (such as specific growth factor signals)

#### Pattern 2: Many-to-One Competition
- `Ligand_C → Receptor_C` (Pathway_3)
- `Ligand_D → Receptor_C` (Pathway_3)

**Why need competition pattern?**
- In real biological systems, multiple ligands often compete for the same receptor
- Tests COMMOT's ability to handle complex interactions
- Simulates signal integration and competitive inhibition phenomena

### 2.3 Expression Pattern Design

**Functional positioning of cell types:**

| Cell Type | Functional Role | Expression Characteristics |
|-----------|----------------|---------------------------|
| CellType_0 | Ligand_A sender | High Ligand_A, low Receptor_A |
| CellType_1 | Ligand_A receiver | Low Ligand_A, high Receptor_A |
| CellType_2 | Multi-ligand sender | High Ligand_B and Ligand_C |
| CellType_3 | Multi-receptor receiver | High Receptor_B and Receptor_C |

**Why this allocation?**

1. **Clear sender-receiver relationship**:
   - CellType_0 → CellType_1 forms directional communication
   - Facilitates verification of whether COMMOT can identify signal direction

2. **Competitive signal testing**:
   - CellType_2 simultaneously sends Ligand_C and Ligand_D
   - CellType_3 needs to integrate signals from two ligands
   - Tests COMMOT's signal decomposition ability

3. **Expression level differences**:
   - High expression: 5.0 ± 0.5
   - Low expression (background): 0.1 ± 0.1
   - Signal-to-noise ratio approximately 50:1, ensuring signal detectability

### 2.4 Ground-Truth Communication Matrix Generation

**Communication strength calculation formula:**
```
Communication(i→j) = Ligand_expr(i) × Receptor_expr(j) × Distance_decay(d_ij)
```

**Distance decay function:**
```
Distance_decay(d) = exp(-d² / (2σ²))
```
where σ = communication_radius / 2

**Why use Gaussian decay?**
1. **Biological plausibility**: Ligand molecules propagate through diffusion, concentration decays with distance in Gaussian manner
2. **Smoothness**: Avoids hard cutoff, closer to reality
3. **Controllability**: Communication range can be controlled by adjusting σ

**Reason for setting communication radius to 15 units:**
- Relative to 100×100 space, 15 units is approximately 15% of spatial scale
- Large enough to allow cross-cluster communication
- Small enough to maintain spatial locality

### 2.5 Noise and Dropout Design

**Why add noise?**
Real spatial transcriptomics data has multiple types of technical noise:
1. **Sequencing noise**: Random fluctuations in read counts
2. **Technical variation**: Batch effects from sample processing
3. **Biological variation**: Natural differences between cells

**Noise parameters:**
- **Gaussian noise level**: 30% (noise_level=0.3)
  - Standard deviation relative to signal strength
  - Maintains signal-to-noise ratio in acceptable range
  
- **Dropout rate**: 10% (dropout_rate=0.1)
  - Simulates sparsity of gene expression
  - Reflects technical limitations of single-cell sequencing

**Why choose these parameters?**
- 30% noise: Sufficient to introduce challenge, but won't completely mask signal
- 10% dropout: Close to dropout level of real Visium data
- Can test algorithm robustness by adjusting these parameters

### 2.6 Background Gene Design

**Reasons for adding 100 background genes:**
1. **Dimensional realism**: Real data contains thousands of genes, most unrelated to CCC
2. **Noise dimensions**: Tests whether algorithm can extract useful signals from noise
3. **Computational complexity**: Simulates computational burden of real data

**Background gene expression distribution:**
- Uses exponential distribution: Exp(λ=0.5)
- Simulates long-tail distribution of gene expression
- Most genes low expression, few genes high expression

## 3. Dataset Workflow for COMMOT Analysis

### 3.1 Input to COMMOT

```python
import commot as ct

# 1. Prepare ligand-receptor pair list
lr_pairs = [
    ['Ligand_A', 'Receptor_A', 'Pathway_1'],
    ['Ligand_B', 'Receptor_B', 'Pathway_2'],
    ['Ligand_C', 'Receptor_C', 'Pathway_3'],
    ['Ligand_D', 'Receptor_C', 'Pathway_3'],
]
df_ligrec = pd.DataFrame(lr_pairs)

# 2. Run COMMOT analysis
ct.tl.spatial_communication(
    adata,
    database_name='simulated',
    df_ligrec=df_ligrec,
    dis_thr=15,  # Consistent with simulated communication radius
    heteromeric=True,
    pathway_sum=True
)
```

### 3.2 COMMOT Output Features

COMMOT generates the following features that can be used for clustering:

1. **Sender signal features** (`adata.obsm['commot-simulated-sum-sender']`)
   - Total amount of signal each cell sends through each LR pair
   - Dimensions: n_cells × n_lr_pairs
   - Reflects cell's "sender" identity

2. **Receiver signal features** (`adata.obsm['commot-simulated-sum-receiver']`)
   - Total amount of signal each cell receives through each LR pair
   - Dimensions: n_cells × n_lr_pairs
   - Reflects cell's "receiver" identity

3. **Pathway-level features**
   - Sender/receiver signals aggregated by signaling pathway
   - Higher-level functional features

### 3.3 Clustering Based on CCC Features

**Why use CCC features for clustering?**
1. **Functional similarity**: Cells with similar communication patterns may belong to the same functional population
2. **Spatial context**: Integrates spatial neighborhood information
3. **Biological significance**: Reflects functional role of cells in tissue

**Clustering method selection:**
- **Leiden algorithm**: Graph-based community detection, suitable for identifying hierarchical structure
- **Louvain algorithm**: Predecessor of Leiden, classic community detection method
- **K-means**: Distance-based clustering, as baseline

**Feature construction strategies:**
```python
# Strategy 1: Use only receiver signals
features = adata.obsm['commot-simulated-sum-receiver']

# Strategy 2: Use only sender signals
features = adata.obsm['commot-simulated-sum-sender']

# Strategy 3: Combine sender and receiver signals
features = np.hstack([
    adata.obsm['commot-simulated-sum-sender'],
    adata.obsm['commot-simulated-sum-receiver']
])

# Strategy 4: Use pathway-level features
# Extract pathway-related columns
```

## 4. Evaluation Metrics

### 4.1 Clustering Performance Evaluation

**ARI (Adjusted Rand Index)**
- Range: [-1, 1], 1 indicates perfect match
- Adjusted for random clustering influence
- Suitable for clusters of different sizes

**NMI (Normalized Mutual Information)**
- Range: [0, 1], 1 indicates perfect match
- Based on information theory
- Insensitive to number of clusters

### 4.2 Communication Reconstruction Evaluation

**Signal strength correlation**
- Compare COMMOT-inferred signal strength with ground-truth
- Use Pearson or Spearman correlation coefficient

**Spatial communication region identification**
- Compare identified communication-active regions with ground-truth
- Use spatial overlap metrics

## 5. Dataset Scalability

### 5.1 Parameter Adjustability

All key parameters can be adjusted to create benchmarks of different difficulty:

```python
simulator = SpatialCCCSimulator(
    n_cells=500,              # Cell count
    n_cell_types=4,           # Number of cell types
    spatial_dim=(100, 100),   # Spatial size
    communication_radius=15,  # Communication range
    random_seed=42            # Reproducibility
)
```

### 5.2 Complexity Variation

**Simple scenarios**:
- Few cell types (2-3 types)
- Large communication radius
- Low noise level
- Clear spatial separation

**Difficult scenarios**:
- Many cell types (6-8 types)
- Small communication radius
- High noise and dropout
- Cell types spatially mixed

### 5.3 Future Extension Directions

1. **3D space**: Extend to three-dimensional tissue structures
2. **Time series**: Add temporal dimension, simulate dynamic communication
3. **More complex communication patterns**:
   - Feedback loops
   - Cascade signaling
   - Paracrine vs autocrine
4. **Heterogeneity**: Expression heterogeneity within same cell type

## 6. Summary

### 6.1 Design Advantages

1. **Known ground-truth**: Can precisely evaluate algorithm performance
2. **Biologically plausible**: Simulates realistic spatial and expression patterns
3. **Controllable complexity**: Adjust difficulty through parameters
4. **Standardized**: Provides fair comparison platform for different algorithms

### 6.2 Correspondence with Project Objectives

| Project Requirement | Implementation |
|---------------------|----------------|
| Spatial coordinates | Gaussian clustering generation |
| Cell type labels | 4 types, clearly defined |
| Ligand expression | Type-specific high expression |
| Receptor expression | Type-specific high expression |
| Communication radius | 15 units, adjustable |
| One-to-one signaling | Ligand_A-Receptor_A, Ligand_B-Receptor_B |
| Many-to-one competition | Ligand_C/D → Receptor_C |
| Ground-truth | Complete communication matrix |
| Noise and dropout | 30% noise, 10% dropout |
| Clustering evaluation | ARI and NMI metrics |

### 6.3 Usage Recommendations

1. **First use**: Generate data with default parameters, familiarize with workflow
2. **Algorithm testing**: Gradually increase noise and complexity
3. **Robustness testing**: Change communication radius and dropout rate
4. **Visualization verification**: Always check spatial distribution of generated data

## 7. Code Usage Example

```python
# 1. Generate data
from generate_synthetic_data import SpatialCCCSimulator

simulator = SpatialCCCSimulator(
    n_cells=500,
    n_cell_types=4,
    spatial_dim=(100, 100),
    communication_radius=15,
    random_seed=42
)

adata = simulator.generate_dataset(
    add_noise=True,
    save_path='synthetic_data.h5ad'
)

# 2. Run COMMOT
import commot as ct
import pandas as pd

lr_pairs = [
    ['Ligand_A', 'Receptor_A', 'Pathway_1'],
    ['Ligand_B', 'Receptor_B', 'Pathway_2'],
    ['Ligand_C', 'Receptor_C', 'Pathway_3'],
    ['Ligand_D', 'Receptor_C', 'Pathway_3'],
]
df_ligrec = pd.DataFrame(lr_pairs)

ct.tl.spatial_communication(
    adata,
    database_name='simulated',
    df_ligrec=df_ligrec,
    dis_thr=15,
    heteromeric=True,
    pathway_sum=True
)

# 3. Extract CCC features for clustering
import scanpy as sc

# Use receiver signal features
receiver_features = adata.obsm['commot-simulated-sum-receiver']

# Build neighbor graph
sc.pp.neighbors(adata, use_rep='commot-simulated-sum-receiver')

# Leiden clustering
sc.tl.leiden(adata, resolution=0.5)

# 4. Evaluate clustering results
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

ari = adjusted_rand_score(adata.obs['cell_type'], adata.obs['leiden'])
nmi = normalized_mutual_info_score(adata.obs['cell_type'], adata.obs['leiden'])

print(f"ARI: {ari:.3f}")
print(f"NMI: {nmi:.3f}")
```

## 8. References

1. Cang, Z., et al. (2023). Screening cell–cell communication in spatial transcriptomics via collective optimal transport. Nature methods, 20(2), 218-228.
2. Traag, V. A., et al. (2019). From Louvain to Leiden: guaranteeing well-connected communities. Scientific reports, 9(1), 5233.
