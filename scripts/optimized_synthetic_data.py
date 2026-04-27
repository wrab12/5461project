"""
Optimized Spatial Transcriptomics Data Generator
Performance improvements:
1. Increase cell count to 800-1000
2. Enhance expression difference (10x instead of 5x)
3. Fix gradient pattern bug
4. Reduce heterogeneity
5. Optimize COMMOT parameters
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import cdist
from typing import Literal, List, Tuple
import warnings
warnings.filterwarnings('ignore')


class OptimizedSpatialCCCSimulator:
    """Optimized Spatial Cell-Cell Communication Simulator"""

    def __init__(
        self,
        n_cells: int = 800,
        n_cell_types: int = 4,
        spatial_dim: Tuple[int, int] = (100, 100),
        communication_radius: float = 20.0,  # Increase communication radius
        spatial_pattern: Literal['clustered', 'striped', 'gradient', 'mixed', 'tissue'] = 'clustered',
        heterogeneity: float = 0.15,  # Reduce heterogeneity
        expression_strength: float = 10.0,  # Enhance expression difference
        random_seed: int = 42
    ):
        self.n_cells = n_cells
        self.n_cell_types = n_cell_types
        self.spatial_dim = spatial_dim
        self.communication_radius = communication_radius
        self.spatial_pattern = spatial_pattern
        self.heterogeneity = heterogeneity
        self.expression_strength = expression_strength
        self.random_seed = random_seed

        np.random.seed(random_seed)

        # Simplified Ligand-Receptor pairs (clearer signals)
        self.lr_pairs = [
            ('Ligand_A', 'Receptor_A', 'Pathway_1'),
            ('Ligand_B', 'Receptor_B', 'Pathway_2'),
            ('Ligand_C', 'Receptor_C', 'Pathway_3'),
            ('Ligand_D', 'Receptor_D', 'Pathway_4'),
        ]

        self.n_lr_pairs = len(self.lr_pairs)
        self.n_background_genes = 100
        self.n_genes = len(set([lr[0] for lr in self.lr_pairs] +
                               [lr[1] for lr in self.lr_pairs])) + self.n_background_genes

    def generate_spatial_coordinates(self) -> Tuple[np.ndarray, np.ndarray]:
        """Generate spatial coordinates"""
        n_per_type = self.n_cells // self.n_cell_types
        coords = np.zeros((self.n_cells, 2))
        cell_types = np.zeros(self.n_cells, dtype=int)

        if self.spatial_pattern == 'clustered':
            # Clustered distribution, increase clustering
            centers = self._get_cluster_centers()
            for i in range(self.n_cell_types):
                start_idx = i * n_per_type
                end_idx = start_idx + n_per_type
                center = centers[i]
                # Reduce standard deviation, increase clustering
                coords[start_idx:end_idx] = np.random.randn(n_per_type, 2) * 8 + center
                cell_types[start_idx:end_idx] = i

        elif self.spatial_pattern == 'striped':
            # Striped distribution
            stripe_width = self.spatial_dim[1] / self.n_cell_types
            for i in range(self.n_cell_types):
                start_idx = i * n_per_type
                end_idx = start_idx + n_per_type
                x = np.random.uniform(0, self.spatial_dim[0], n_per_type)
                y_center = (i + 0.5) * stripe_width
                y = np.random.randn(n_per_type) * 4 + y_center  # Reduce Y direction variation
                coords[start_idx:end_idx] = np.column_stack([x, y])
                cell_types[start_idx:end_idx] = i

        elif self.spatial_pattern == 'gradient':
            # FixGradient distribution
            all_x = np.random.uniform(0, self.spatial_dim[0], self.n_cells)
            all_y = np.random.uniform(0, self.spatial_dim[1], self.n_cells)
            coords = np.column_stack([all_x, all_y])

            # Fix: Based on X coordinate assign cell types, ensure each type has cells
            x_normalized = all_x / self.spatial_dim[0]

            # Use piecewise function instead of sigmoid
            for i in range(self.n_cells):
                # Divide X axis into n_cell_types segments, each segment corresponds to one cell type
                # But has 20% overlap at boundaries
                x_val = x_normalized[i]
                segment_width = 1.0 / self.n_cell_types

                # Find main cell type
                main_type = int(x_val / segment_width)
                main_type = min(main_type, self.n_cell_types - 1)

                # At boundaries has probability to assign to adjacent type
                boundary_dist = (x_val % segment_width) / segment_width
                if boundary_dist < 0.2 and main_type > 0:
                    # Near left boundary, may be previous type
                    if np.random.random() < 0.3:
                        main_type -= 1
                elif boundary_dist > 0.8 and main_type < self.n_cell_types - 1:
                    # Near right boundary, may be next type
                    if np.random.random() < 0.3:
                        main_type += 1

                cell_types[i] = main_type

        elif self.spatial_pattern == 'mixed':
            # Mixed distribution: 70% clustered, 30% random
            n_clustered = int(self.n_cells * 0.7)
            n_random = self.n_cells - n_clustered

            centers = self._get_cluster_centers()
            n_per_cluster = n_clustered // self.n_cell_types
            for i in range(self.n_cell_types):
                start_idx = i * n_per_cluster
                end_idx = start_idx + n_per_cluster
                center = centers[i]
                coords[start_idx:end_idx] = np.random.randn(n_per_cluster, 2) * 7 + center
                cell_types[start_idx:end_idx] = i

            coords[n_clustered:] = np.random.uniform([0, 0], self.spatial_dim, (n_random, 2))
            cell_types[n_clustered:] = np.random.randint(0, self.n_cell_types, n_random)

        elif self.spatial_pattern == 'tissue':
            # Tissue-like structure
            center = np.array(self.spatial_dim) / 2

            for i in range(self.n_cells):
                ct = i % self.n_cell_types
                cell_types[i] = ct

                if ct == 0:  # Core
                    r = np.random.rayleigh(8)
                    theta = np.random.uniform(0, 2*np.pi)
                    coords[i] = center + r * np.array([np.cos(theta), np.sin(theta)])

                elif ct == 1:  # Surrounding
                    r = np.random.uniform(12, 25)
                    theta = np.random.uniform(0, 2*np.pi)
                    coords[i] = center + r * np.array([np.cos(theta), np.sin(theta)])

                elif ct == 2:  # Edge
                    r = np.random.uniform(25, 40)
                    theta = np.random.uniform(0, 2*np.pi)
                    coords[i] = center + r * np.array([np.cos(theta), np.sin(theta)])

                else:  # Vessel-like
                    vessel_id = i % 3
                    t = np.random.uniform(0, 1)
                    if vessel_id == 0:
                        coords[i] = [t * self.spatial_dim[0], self.spatial_dim[1] * 0.3]
                    elif vessel_id == 1:
                        coords[i] = [t * self.spatial_dim[0], self.spatial_dim[1] * 0.7]
                    else:
                        coords[i] = [self.spatial_dim[0] * 0.5, t * self.spatial_dim[1]]
                    coords[i] += np.random.randn(2) * 2

        coords = np.clip(coords, [0, 0], self.spatial_dim)
        return coords, cell_types

    def _get_cluster_centers(self) -> List[np.ndarray]:
        """Get cluster centers"""
        if self.n_cell_types == 4:
            return [
                np.array([self.spatial_dim[0]*0.25, self.spatial_dim[1]*0.25]),
                np.array([self.spatial_dim[0]*0.75, self.spatial_dim[1]*0.25]),
                np.array([self.spatial_dim[0]*0.25, self.spatial_dim[1]*0.75]),
                np.array([self.spatial_dim[0]*0.75, self.spatial_dim[1]*0.75]),
            ]
        else:
            centers = []
            n_cols = int(np.ceil(np.sqrt(self.n_cell_types)))
            for i in range(self.n_cell_types):
                row = i // n_cols
                col = i % n_cols
                x = (col + 0.5) * self.spatial_dim[0] / n_cols
                y = (row + 0.5) * self.spatial_dim[1] / n_cols
                centers.append(np.array([x, y]))
            return centers

    def generate_expression_patterns(
        self,
        cell_types: np.ndarray,
        pattern: Literal['simple', 'cascade'] = 'simple'
    ) -> np.ndarray:
        """Generate gene expression pattern"""
        expression = np.ones((self.n_cells, self.n_genes))

        ligands = [lr[0] for lr in self.lr_pairs]
        receptors = [lr[1] for lr in self.lr_pairs]
        all_lr_genes = list(set(ligands + receptors))

        if pattern == 'simple':
            # Simple pattern: clear Sender-Receiver relationship
            for i, (lig, rec, pathway) in enumerate(self.lr_pairs):
                lig_idx = all_lr_genes.index(lig)
                rec_idx = all_lr_genes.index(rec)

                sender_ct = i % self.n_cell_types
                receiver_ct = (i + 1) % self.n_cell_types

                sender_mask = cell_types == sender_ct
                receiver_mask = cell_types == receiver_ct

                # Use stronger expression difference
                expression[sender_mask, lig_idx] = self.expression_strength
                expression[receiver_mask, rec_idx] = self.expression_strength

                # Other cell types maintain low expression
                other_mask = ~(sender_mask | receiver_mask)
                expression[other_mask, lig_idx] = 1.0
                expression[other_mask, rec_idx] = 1.0

        elif pattern == 'cascade':
            # cascadeSignal
            cascade_pairs = [
                ('Ligand_A', 'Receptor_A', 0, 1),
                ('Ligand_B', 'Receptor_B', 1, 2),
                ('Ligand_C', 'Receptor_C', 2, 3),
                ('Ligand_D', 'Receptor_D', 3, 0),
            ]

            for lig, rec, sender_ct, receiver_ct in cascade_pairs:
                if lig in all_lr_genes and rec in all_lr_genes:
                    lig_idx = all_lr_genes.index(lig)
                    rec_idx = all_lr_genes.index(rec)

                    sender_mask = cell_types == sender_ct
                    receiver_mask = cell_types == receiver_ct

                    expression[sender_mask, lig_idx] = self.expression_strength
                    expression[receiver_mask, rec_idx] = self.expression_strength

        # Add moderate heterogeneity
        if self.heterogeneity > 0:
            noise = np.random.randn(self.n_cells, len(all_lr_genes)) * self.heterogeneity
            expression[:, :len(all_lr_genes)] *= (1 + noise)
            expression = np.maximum(expression, 0.1)

        # Background genes
        expression[:, len(all_lr_genes):] = np.random.lognormal(0, 0.5,
                                                                 (self.n_cells, self.n_background_genes))

        return expression

    def generate_ground_truth_communication(
        self,
        coords: np.ndarray,
        expression: np.ndarray,
        cell_types: np.ndarray
    ) -> dict:
        """Generate ground-truth communication matrix"""
        distances = cdist(coords, coords, metric='euclidean')

        ligands = [lr[0] for lr in self.lr_pairs]
        receptors = [lr[1] for lr in self.lr_pairs]
        all_lr_genes = list(set(ligands + receptors))

        comm_matrices = {}

        for lig, rec, pathway in self.lr_pairs:
            lig_idx = all_lr_genes.index(lig)
            rec_idx = all_lr_genes.index(rec)

            comm_matrix = np.zeros((self.n_cells, self.n_cells))

            for i in range(self.n_cells):
                for j in range(self.n_cells):
                    if i != j:
                        distance = distances[i, j]

                        if distance < self.communication_radius:
                            distance_factor = np.exp(-distance**2 / (2 * self.communication_radius**2))
                            sender_expr = expression[i, lig_idx]
                            receiver_expr = expression[j, rec_idx]
                            comm_matrix[i, j] = sender_expr * receiver_expr * distance_factor

            comm_matrices[f'{lig}-{rec}'] = comm_matrix

        return comm_matrices

    def add_noise_and_dropout(
        self,
        expression: np.ndarray,
        noise_level: float = 0.2,
        dropout_rate: float = 0.05
    ) -> np.ndarray:
        """Add moderate noise"""
        noise = np.random.randn(*expression.shape) * noise_level
        expression_noisy = expression * (1 + noise)
        expression_noisy = np.maximum(expression_noisy, 0)

        dropout_mask = np.random.random(expression.shape) > dropout_rate
        expression_noisy = expression_noisy * dropout_mask

        return expression_noisy

    def create_anndata(
        self,
        coords: np.ndarray,
        expression: np.ndarray,
        cell_types: np.ndarray,
        comm_matrices: dict
    ) -> sc.AnnData:
        """Create AnnData object"""
        ligands = [lr[0] for lr in self.lr_pairs]
        receptors = [lr[1] for lr in self.lr_pairs]
        all_lr_genes = list(set(ligands + receptors))

        gene_names = all_lr_genes + [f'Gene_{i}' for i in range(self.n_background_genes)]

        adata = sc.AnnData(X=expression)
        adata.var_names = gene_names
        adata.obs_names = [f'Cell_{i}' for i in range(self.n_cells)]

        adata.obs['cell_type'] = [f'CellType_{ct}' for ct in cell_types]
        adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

        adata.obsm['spatial'] = coords
        adata.uns['ground_truth_communication'] = comm_matrices

        adata.uns['simulation_params'] = {
            'n_cells': self.n_cells,
            'n_cell_types': self.n_cell_types,
            'spatial_dim': list(self.spatial_dim),
            'communication_radius': self.communication_radius,
            'spatial_pattern': self.spatial_pattern,
            'heterogeneity': self.heterogeneity,
            'expression_strength': self.expression_strength,
            'n_lr_pairs': self.n_lr_pairs,
            'random_seed': self.random_seed
        }

        return adata

    def generate_dataset(
        self,
        expression_pattern: Literal['simple', 'cascade'] = 'simple',
        add_noise: bool = True,
        noise_level: float = 0.2,
        dropout_rate: float = 0.05,
        save_path: str = None
    ) -> sc.AnnData:
        """Generate complete dataset"""
        print(f"Generate spatial coordinates (pattern: {self.spatial_pattern})...")
        coords, cell_types = self.generate_spatial_coordinates()

        print(f"Generate gene expression (pattern: {expression_pattern})...")
        expression = self.generate_expression_patterns(cell_types, pattern=expression_pattern)

        print("Generate ground-truth communication matrix...")
        comm_matrices = self.generate_ground_truth_communication(coords, expression, cell_types)

        if add_noise:
            print(f"Add noise (level: {noise_level}, dropout: {dropout_rate})...")
            expression = self.add_noise_and_dropout(expression, noise_level, dropout_rate)

        print("Create AnnData object...")
        adata = self.create_anndata(coords, expression, cell_types, comm_matrices)

        if save_path:
            print(f"Save to: {save_path}")
            adata.write_h5ad(save_path)

        print("Dataset generation complete!")
        print(f"Number of cells: {adata.n_obs}")
        print(f"Number of genes: {adata.n_vars}")
        print(f"Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")

        return adata


if __name__ == '__main__':
    print("="*80)
    print(" "*20 + "Optimized Data Generator Test")
    print("="*80)

    # GenerateOptimizedDataset
    examples = [
        {
            'name': 'Optimized_Clustered',
            'params': {
                'n_cells': 800,
                'spatial_pattern': 'clustered',
                'expression_strength': 10.0,
                'heterogeneity': 0.15,
            },
            'expression_pattern': 'simple'
        },
        {
            'name': 'Optimized_Tissue',
            'params': {
                'n_cells': 1000,
                'spatial_pattern': 'tissue',
                'expression_strength': 12.0,
                'heterogeneity': 0.15,
            },
            'expression_pattern': 'cascade'
        },
        {
            'name': 'Optimized_Striped',
            'params': {
                'n_cells': 800,
                'spatial_pattern': 'striped',
                'communication_radius': 30.0,  # Increase communication radius
                'expression_strength': 10.0,
                'heterogeneity': 0.12,
            },
            'expression_pattern': 'simple'
        },
        {
            'name': 'Optimized_Gradient',
            'params': {
                'n_cells': 1000,
                'spatial_pattern': 'gradient',
                'expression_strength': 10.0,
                'heterogeneity': 0.18,
            },
            'expression_pattern': 'simple'
        },
        {
            'name': 'Optimized_Mixed',
            'params': {
                'n_cells': 1000,
                'spatial_pattern': 'mixed',
                'expression_strength': 10.0,
                'heterogeneity': 0.2,
            },
            'expression_pattern': 'simple'
        },
    ]

    for example in examples:
        print("\n" + "="*80)
        print(f"Generate: {example['name']}")
        print("="*80)

        params = example['params']
        expression_pattern = example['expression_pattern']

        simulator = OptimizedSpatialCCCSimulator(**params)
        adata = simulator.generate_dataset(
            expression_pattern=expression_pattern,
            add_noise=True,
            save_path=f"{example['name']}.h5ad"
        )

    print("\n" + "="*80)
    print("All optimized datasets generated successfully!")
    print("="*80)
