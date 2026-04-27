# COMMOT Benchmark Project Report

## Project Overview

This project is for the CSCI 5461 course, aiming to create simulated datasets for the COMMOT (COMMunication analysis by Optimal Transport) method and conduct clustering algorithm benchmark testing.

**Core Workflow**: Simulated data generation → COMMOT analysis → Extract CCC features → Clustering → Performance evaluation

## Team Division

- **Rui**: Data generation module
- **Evan**: Leiden clustering algorithm implementation
- **Grace**: Louvain clustering algorithm implementation
- **Lechen**: K-means clustering algorithm implementation

## I. Simulated Dataset Design

### 1.1 Design Objectives

Create a spatial transcriptomics dataset with clear ground-truth to test whether COMMOT can correctly identify cell-cell communication patterns and evaluate clustering algorithm performance based on these communication features.

### 1.2 Dataset Parameters

```python
n_cells = 200              # Cell count
n_cell_types = 4           # Number of cell types
spatial_dim = (100, 100)   # Spatial range
communication_radius = 15  # Communication radius
```

### 1.3 Design Principles

#### 1.3.1 Spatial Structure Design

**Why this design**:
- Each cell type exhibits circular aggregation distribution in space
- Simulates spatial organization patterns of cells in real tissues
- Clear spatial boundaries between different cell types, facilitating clustering validation

**Implementation**:
```python
# Each cell type occupies one quadrant
centers = [(25, 25), (75, 25), (25, 75), (75, 75)]
# Gaussian sampling around center point
coords = np.random.randn(n_per_type, 2) * 10 + center
```

**Biological Significance**:
- In real tissues, cells of the same type tend to cluster together
- Spatial proximity is a prerequisite for cell-cell communication

#### 1.3.2 Ligand-Receptor Pair Design

**Designed 4 ligand-receptor interactions**:

1. **Ligand_A → Receptor_A (Pathway_1)**
   - One-to-one signaling
   - CellType_0 highly expresses Ligand_A (sender)
   - CellType_1 highly expresses Receptor_A (receiver)

2. **Ligand_B → Receptor_B (Pathway_2)**
   - One-to-one signaling
   - CellType_1 highly expresses Ligand_B (sender)
   - CellType_2 highly expresses Receptor_B (receiver)

3. **Ligand_C → Receptor_C (Pathway_3)**
   - Many-to-one competitive signaling (first ligand)
   - CellType_2 highly expresses Ligand_C (sender)
   - CellType_3 highly expresses Receptor_C (receiver)

4. **Ligand_D → Receptor_C (Pathway_3)**
   - Many-to-one competitive signaling (second ligand)
   - CellType_3 highly expresses Ligand_D (sender)
   - CellType_3 highly expresses Receptor_C (receiver, autocrine)

**Why this design**:

1. **One-to-one signaling**:
   - Simplest communication pattern
   - Validates whether COMMOT can identify basic ligand-receptor interactions
   - Forms clear sender-receiver relationships

2. **Many-to-one competition**:
   - Simulates complexity in real biological systems
   - Multiple ligands compete for the same receptor
   - Tests COMMOT's ability to handle complex interactions

3. **Autocrine signaling**:
   - Both Ligand_D and Receptor_C highly expressed in CellType_3
   - Simulates cellular self-regulation mechanisms
   - Increases biological realism of dataset

#### 1.3.3 Expression Pattern Design

**Cell type-specific expression**:

```python
# CellType_0: High expression of Ligand_A
expression[cell_type_mask, ligand_idx] = base_expression * 5.0

# CellType_1: High expression of Receptor_A and Ligand_B
expression[cell_type_mask, receptor_idx] = base_expression * 5.0
```

**Why this design**:
- 5-fold expression difference is sufficiently obvious for COMMOT identification
- Base expression level (base_expression = 1.0) ensures all genes have some expression
- Avoids complete 0/1 expression, more consistent with real transcriptome data

#### 1.3.4 Ground-Truth Communication Matrix

**Distance decay model**:

```python
# Gaussian distance decay
distance_factor = np.exp(-distances**2 / (2 * communication_radius**2))

# Communication strength = sender expression × receiver expression × distance decay
comm_matrix[i, j] = sender_expr * receiver_expr * distance_factor
```

**Why this design**:

1. **Gaussian decay function**:
   - Simulates physical process of ligand diffusion
   - Greater distance, weaker communication strength
   - Consistent with diffusion laws in biology

2. **Expression level dependence**:
   - Sender must express ligand
   - Receiver must express receptor
   - Both are indispensable

3. **Communication radius**:
   - Set to 15 units
   - Ensures cells within same cell type can communicate
   - Communication between different cell types is limited

### 1.4 Noise Simulation

**Why noise is needed**:
- Real transcriptome data always has technical noise
- Tests algorithm robustness
- More realistically simulates experimental data

**Noise types**:
1. **Gaussian noise** (30%): Simulates sequencing errors
2. **Dropout** (10%): Simulates gene detection failures

**Note**: Current version temporarily disables noise for stability (`add_noise=False`)

## II. COMMOT Analysis Workflow

### 2.1 Role of COMMOT

COMMOT is a spatial cell-cell communication inference tool based on optimal transport theory. Its core functions are:

1. **Input**:
   - Spatial transcriptomics data (gene expression + spatial coordinates)
   - Ligand-receptor pair list

2. **Output**:
   - Sender signal strength for each cell
   - Receiver signal strength for each cell

3. **Principle**:
   - Uses optimal transport theory to calculate communication flow between cells
   - Considers spatial distance and expression levels
   - Infers which cells are sending signals and which are receiving

### 2.2 Why Use COMMOT

1. **Spatial information integration**:
   - Considers not only gene expression but also spatial position of cells
   - Only spatially proximate cells can communicate

2. **Quantitative analysis**:
   - Provides quantified communication strength values for each cell
   - Can compare communication activity levels of different cells

3. **Feature extraction**:
   - Converts high-dimensional gene expression data to low-dimensional communication features
   - These features more directly reflect functional states of cells

### 2.3 COMMOT Parameter Settings

```python
ct.tl.spatial_communication(
    adata,
    database_name='simulated',
    df_ligrec=df_ligrec,        # Ligand-receptor pairs
    dis_thr=15,                  # Communication distance threshold
    heteromeric=False,           # Don't consider heteromeric complexes
    pathway_sum=False            # Don't perform pathway-level aggregation
)
```

**Parameter explanation**:
- `dis_thr=15`: Consistent with communication radius during data generation
- `heteromeric=False`: Simplifies analysis, only considers individual ligand-receptor pairs
- `pathway_sum=False`: Retains independent signals for each ligand-receptor pair

### 2.4 COMMOT Output Features

**Generated feature matrices**:
- `commot-simulated-sum-sender`: Sender signal matrix (200 cells × 5 LR pairs)
- `commot-simulated-sum-receiver`: Receiver signal matrix (200 cells × 5 LR pairs)
- **Combined features**: (200 cells × 10 features)

**Feature meaning**:
- Each cell has 10 CCC features
- First 5 features: Signal strength as sender
- Last 5 features: Signal strength as receiver
- These features capture communication roles of cells

## III. Clustering Algorithm Implementation

### 3.1 Why Cluster on CCC Features

**Core hypothesis**:
- Cells of the same type have similar communication patterns
- Communication features can reflect functional states of cells
- Clustering based on communication features should recover cell types

**Difference from traditional clustering**:
- Traditional: Clustering based on gene expression
- This project: Clustering based on COMMOT-extracted communication features
- Advantage: Communication features integrate spatial and expression information

### 3.2 Leiden Clustering (Evan)

**Algorithm principle**:
- Graph-based community detection algorithm
- Improved version of Louvain algorithm
- Guarantees generation of connected communities

**Implementation**:
```python
sc.pp.neighbors(adata_temp, n_neighbors=15, use_rep='X')
sc.tl.leiden(adata_temp, resolution=0.5)
```

**Parameters**:
- `n_neighbors=15`: Construct k-nearest neighbor graph
- `resolution=0.5`: Controls clustering granularity

**Results**:
- Identified 6 clusters
- ARI = 0.4569
- NMI = 0.5924

### 3.3 Louvain Clustering (Grace)

**Algorithm principle**:
- Classic community detection algorithm
- Optimizes modularity
- Fast and widely used

**Implementation**:
```python
sc.pp.neighbors(adata_temp, n_neighbors=15, use_rep='X')
sc.tl.louvain(adata_temp, resolution=0.5)
```

**Results**:
- Identified 5 clusters
- ARI = 0.4732 (best)
- NMI = 0.5675

**Why best performance**:
- Cluster count (5) close to true cell type count (4)
- Modularity optimization suits this dataset's structure

### 3.4 K-means Clustering (Lechen)

**Algorithm principle**:
- Distance-based clustering method
- Minimizes within-cluster sum of squares
- Requires pre-specifying cluster count

**Implementation**:
```python
kmeans = KMeans(n_clusters=4, random_state=42, n_init=10)
adata.obs['kmeans'] = kmeans.fit_predict(features)
```

**Parameters**:
- `n_clusters=4`: Consistent with true cell type count
- `n_init=10`: Run 10 times and take best result

**Results**:
- Identified 4 clusters
- ARI = 0.0028 (worst)
- NMI = 0.0854

**Why poorer performance**:
- K-means assumes spherical clusters
- CCC feature space may not satisfy this assumption
- Sensitive to initialization

## IV. Performance Evaluation

### 4.1 Evaluation Metrics

#### ARI (Adjusted Rand Index)
- **Range**: [-1, 1], 1 indicates perfect match
- **Meaning**: Measures consistency between clustering results and true labels
- **Advantage**: Adjusted for random clustering influence

#### NMI (Normalized Mutual Information)
- **Range**: [0, 1], 1 indicates perfect match
- **Meaning**: Measures information overlap between clustering results and true labels
- **Advantage**: Insensitive to cluster count

### 4.2 Results Analysis

| Method | ARI | NMI | Cluster Count | Team Member |
|--------|-----|-----|---------------|-------------|
| Leiden | 0.4569 | 0.5924 | 6 | Evan |
| **Louvain** | **0.4732** | **0.5675** | **5** | **Grace** |
| K-means | 0.0028 | 0.0854 | 4 | Lechen |

**Best method**: Louvain (Grace)

**Result interpretation**:

1. **Louvain performs best**:
   - ARI = 0.4732 indicates moderate clustering quality
   - Cluster count (5) close to true value (4)
   - Successfully captured most cell type structure

2. **Leiden performs second**:
   - Performance close to Louvain
   - Generated more clusters (6)
   - May have over-segmented some cell types

3. **K-means performs poorly**:
   - ARI close to 0, almost equivalent to random clustering
   - Indicates CCC feature space not suitable for K-means
   - May need feature transformation or dimensionality reduction

### 4.3 Reasons for Suboptimal Performance

**Why is ARI only around 0.47?**

1. **Small dataset scale**:
   - Only 200 cells
   - Only 50 cells per cell type
   - Weak statistical signal

2. **Simplified COMMOT settings**:
   - Noise disabled (`add_noise=False`)
   - Used simplified parameters (`heteromeric=False`, `pathway_sum=False`)
   - Sacrificed complexity for stability

3. **Complexity of CCC features**:
   - Communication features may not be as direct as gene expression features
   - Requires more cells and stronger signals

4. **Limitations of clustering algorithms**:
   - Graph-based methods depend on neighborhood structure
   - K-means not suitable for non-spherical distributions

## V. Visualization Results

Generated visualization charts include:

1. **Spatial distribution plots** (Row 1):
   - Ground Truth: Spatial distribution of true cell types
   - Leiden, Louvain, K-means: Spatial distribution of three clustering results
   - Can visually compare clustering results with true labels

2. **CCC signal plots** (Row 2):
   - Sender signal heatmap: Shows total sender signal strength for each cell
   - Receiver signal heatmap: Shows total receiver signal strength for each cell
   - Can see communication activity levels of different cell types

3. **Performance comparison plots** (Row 2):
   - Bar chart comparison of ARI and NMI
   - Confusion matrix: Shows correspondence between clustering results and true labels

4. **Statistical summary** (Row 3):
   - Dataset information
   - COMMOT analysis parameters
   - Best method and its performance

## VI. Project Summary

### 6.1 Completed Work

✅ **Data generation module**:
- Created spatial transcriptomics dataset with clear ground-truth
- Designed reasonable ligand-receptor interaction patterns
- Implemented spatial aggregation and distance decay models

✅ **COMMOT analysis**:
- Successfully ran COMMOT spatial communication analysis
- Extracted sender and receiver signal features
- Validated that COMMOT can identify communication patterns

✅ **Clustering algorithm implementation**:
- Evan: Leiden clustering
- Grace: Louvain clustering (best performance)
- Lechen: K-means clustering

✅ **Performance evaluation**:
- Used ARI and NMI metrics
- Generated confusion matrices
- Created comprehensive visualization charts

✅ **Documentation**:
- Detailed design documentation
- Complete analysis report
- Clear code comments

### 6.2 Key Findings

1. **COMMOT effectiveness**:
   - COMMOT successfully extracted cell-cell communication features
   - CCC features can partially recover cell type information

2. **Algorithm comparison**:
   - Graph-based methods (Leiden, Louvain) outperform K-means
   - Louvain performs best on this dataset
   - K-means not suitable for CCC feature space

3. **Dataset design**:
   - Spatial aggregation pattern effective
   - Ligand-receptor pair design reasonable
   - Larger dataset needed for better performance

### 6.3 Improvement Directions

**Short-term improvements**:
1. Increase cell count (500-1000)
2. Add moderate noise
3. Test different COMMOT parameters
4. Try feature dimensionality reduction (PCA, UMAP)

**Long-term extensions**:
1. 3D spatial data
2. Time series data
3. More complex communication networks
4. Real dataset validation

### 6.4 Biological Significance

This project demonstrates:
1. Spatial cell-cell communication is an important feature of cell types
2. COMMOT can extract communication signals from spatial transcriptomics data
3. Clustering based on communication features is feasible
4. Provides benchmark platform for future CCC analysis methods

## VII. File List

### Core Code
1. `generate_synthetic_data.py` - Data generator
2. `FINAL_COMMOT_BENCHMARK.py` - Complete analysis pipeline (final version)

### Output Files
3. `FINAL_COMMOT_BENCHMARK.png` - Comprehensive visualization chart
4. `FINAL_COMMOT_RESULTS.h5ad` - Complete results data
5. `Project_Report.md` - This report

### Documentation
6. `Design_Documentation.md` - Dataset design documentation
7. `Project_Completion_Summary.md` - Project completion summary

## VIII. Usage

### Environment Requirements
```bash
conda activate commot
```

### Run Analysis
```bash
python FINAL_COMMOT_BENCHMARK.py
```

### Output
- `FINAL_COMMOT_BENCHMARK.png`: Visualization results
- `FINAL_COMMOT_RESULTS.h5ad`: Complete data object

## IX. References

1. Cang, Z., et al. (2023). Screening cell–cell communication in spatial transcriptomics via collective optimal transport. *Nature methods*, 20(2), 218-228.

2. Traag, V. A., et al. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific reports*, 9(1), 5233.

3. Blondel, V. D., et al. (2008). Fast unfolding of communities in large networks. *Journal of statistical mechanics: theory and experiment*, 2008(10), P10008.

---

**Project Status**: ✅ Completed  
**Last Updated**: 2026-04-26  
**Version**: Final 1.0
