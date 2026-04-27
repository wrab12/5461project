"""
UMN Letter-shaped Spatial Transcriptomics Data Generator
Creative Test: Use letter shapes as spatial distribution patterns
"""
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')


class LetterShapeSimulator:
    """Letter shape data generator"""

    def __init__(
        self,
        letter: str,
        n_cells: int = 800,
        n_cell_types: int = 4,
        scale: float = 80.0,
        expression_strength: float = 12.0,
        heterogeneity: float = 0.15,
        random_seed: int = 42
    ):
        self.letter = letter.upper()
        self.n_cells = n_cells
        self.n_cell_types = n_cell_types
        self.scale = scale
        self.expression_strength = expression_strength
        self.heterogeneity = heterogeneity
        self.random_seed = random_seed

        np.random.seed(random_seed)

        # Ligand-Receptor pairs
        self.lr_pairs = [
            ('Ligand_A', 'Receptor_A', 'Pathway_1'),
            ('Ligand_B', 'Receptor_B', 'Pathway_2'),
            ('Ligand_C', 'Receptor_C', 'Pathway_3'),
            ('Ligand_D', 'Receptor_D', 'Pathway_4'),
        ]

        self.n_background_genes = 100
        self.n_genes = 8 + self.n_background_genes

    def generate_letter_shape(self, letter: str, n_points: int) -> np.ndarray:
        """Generate letter shape coordinate points"""

        if letter == 'U':
            # U shape: left vertical + bottom arc + right vertical
            coords = []

            # Left vertical (30% of points)
            n_left = int(n_points * 0.3)
            x_left = np.ones(n_left) * 0.2
            y_left = np.linspace(0.3, 1.0, n_left)
            coords.extend(zip(x_left, y_left))

            # Bottom arc (40% of points)
            n_bottom = int(n_points * 0.4)
            theta = np.linspace(np.pi, 0, n_bottom)
            x_bottom = 0.5 + 0.3 * np.cos(theta)
            y_bottom = 0.3 + 0.15 * np.sin(theta)
            coords.extend(zip(x_bottom, y_bottom))

            # Right vertical (30% of points)
            n_right = n_points - n_left - n_bottom
            x_right = np.ones(n_right) * 0.8
            y_right = np.linspace(0.3, 1.0, n_right)
            coords.extend(zip(x_right, y_right))

        elif letter == 'M':
            # M shape: left vertical + left diagonal + right diagonal + right vertical
            coords = []

            # Left vertical (25% of points)
            n_left = int(n_points * 0.25)
            x_left = np.ones(n_left) * 0.1
            y_left = np.linspace(0.0, 1.0, n_left)
            coords.extend(zip(x_left, y_left))

            # Left diagonal (25% of points)
            n_left_diag = int(n_points * 0.25)
            x_left_diag = np.linspace(0.1, 0.5, n_left_diag)
            y_left_diag = np.linspace(1.0, 0.4, n_left_diag)
            coords.extend(zip(x_left_diag, y_left_diag))

            # Right diagonal (25% of points)
            n_right_diag = int(n_points * 0.25)
            x_right_diag = np.linspace(0.5, 0.9, n_right_diag)
            y_right_diag = np.linspace(0.4, 1.0, n_right_diag)
            coords.extend(zip(x_right_diag, y_right_diag))

            # Right vertical (25% of points)
            n_right = n_points - n_left - n_left_diag - n_right_diag
            x_right = np.ones(n_right) * 0.9
            y_right = np.linspace(1.0, 0.0, n_right)
            coords.extend(zip(x_right, y_right))

        elif letter == 'N':
            # N shape: left vertical + diagonal + right vertical
            coords = []

            # Left vertical (33% of points)
            n_left = int(n_points * 0.33)
            x_left = np.ones(n_left) * 0.2
            y_left = np.linspace(0.0, 1.0, n_left)
            coords.extend(zip(x_left, y_left))

            # Diagonal (34% of points)
            n_diag = int(n_points * 0.34)
            x_diag = np.linspace(0.2, 0.8, n_diag)
            y_diag = np.linspace(1.0, 0.0, n_diag)
            coords.extend(zip(x_diag, y_diag))

            # Right vertical (33% of points)
            n_right = n_points - n_left - n_diag
            x_right = np.ones(n_right) * 0.8
            y_right = np.linspace(0.0, 1.0, n_right)
            coords.extend(zip(x_right, y_right))

        elif letter == 'UMN':
            # UMN combination: three letters side by side
            coords = []
            n_per_letter = n_points // 3

            # U (left)
            u_coords = self.generate_letter_shape('U', n_per_letter)
            u_coords = u_coords * np.array([0.3, 1.0]) + np.array([0.0, 0.0])
            coords.extend(u_coords)

            # M (middle)
            m_coords = self.generate_letter_shape('M', n_per_letter)
            m_coords = m_coords * np.array([0.3, 1.0]) + np.array([0.35, 0.0])
            coords.extend(m_coords)

            # N (right)
            n_coords = self.generate_letter_shape('N', n_points - 2*n_per_letter)
            n_coords = n_coords * np.array([0.3, 1.0]) + np.array([0.70, 0.0])
            coords.extend(n_coords)

        else:
            raise ValueError(f"Unsupported letter: {letter}")

        coords = np.array(coords)
        # Scale to specified size
        coords = coords * self.scale

        return coords

    def assign_cell_types_to_shape(self, coords: np.ndarray) -> np.ndarray:
        """Assign cell types based on position"""
        n_points = len(coords)
        cell_types = np.zeros(n_points, dtype=int)

        if self.letter == 'UMN':
            # UMN: each letter is different cell type
            # U: CellType_0, M: CellType_1,2, N: CellType_3
            n_per_letter = n_points // 3

            # U part
            cell_types[:n_per_letter] = 0

            # M part (split into two types)
            m_start = n_per_letter
            m_end = 2 * n_per_letter
            m_mid = (m_start + m_end) // 2
            cell_types[m_start:m_mid] = 1
            cell_types[m_mid:m_end] = 2

            # N part
            cell_types[2*n_per_letter:] = 3

        else:
            # Single letter: assign cell types along shape segments
            segment_size = n_points // self.n_cell_types
            for i in range(self.n_cell_types):
                start = i * segment_size
                end = start + segment_size if i < self.n_cell_types - 1 else n_points
                cell_types[start:end] = i

        return cell_types

    def add_noise_to_coords(self, coords: np.ndarray, noise_level: float = 2.0) -> np.ndarray:
        """Add noise to coordinates to make shape more natural"""
        noise = np.random.randn(*coords.shape) * noise_level
        return coords + noise

    def generate_expression(self, cell_types: np.ndarray) -> np.ndarray:
        """Generate gene expression"""
        n_cells = len(cell_types)
        expression = np.ones((n_cells, self.n_genes))

        # LigandandReceptor
        ligands = ['Ligand_A', 'Ligand_B', 'Ligand_C', 'Ligand_D']
        receptors = ['Receptor_A', 'Receptor_B', 'Receptor_C', 'Receptor_D']

        for i in range(4):
            sender_ct = i % self.n_cell_types
            receiver_ct = (i + 1) % self.n_cell_types

            sender_mask = cell_types == sender_ct
            receiver_mask = cell_types == receiver_ct

            # Ligand expression
            expression[sender_mask, i] = self.expression_strength
            # Receptor expression
            expression[receiver_mask, i + 4] = self.expression_strength

        # Add heterogeneity
        if self.heterogeneity > 0:
            noise = np.random.randn(n_cells, 8) * self.heterogeneity
            expression[:, :8] *= (1 + noise)
            expression = np.maximum(expression, 0.1)

        # Background genes
        expression[:, 8:] = np.random.lognormal(0, 0.5, (n_cells, self.n_background_genes))

        return expression

    def generate_ground_truth_communication(
        self,
        coords: np.ndarray,
        expression: np.ndarray,
        communication_radius: float = 15.0
    ) -> dict:
        """Generate ground-truth communication matrix"""
        distances = cdist(coords, coords, metric='euclidean')
        comm_matrices = {}

        for i, (lig, rec, pathway) in enumerate(self.lr_pairs):
            comm_matrix = np.zeros((len(coords), len(coords)))

            for cell_i in range(len(coords)):
                for cell_j in range(len(coords)):
                    if cell_i != cell_j:
                        distance = distances[cell_i, cell_j]

                        if distance < communication_radius:
                            distance_factor = np.exp(-distance**2 / (2 * communication_radius**2))
                            sender_expr = expression[cell_i, i]
                            receiver_expr = expression[cell_j, i + 4]
                            comm_matrix[cell_i, cell_j] = sender_expr * receiver_expr * distance_factor

            comm_matrices[f'{lig}-{rec}'] = comm_matrix

        return comm_matrices

    def create_anndata(
        self,
        coords: np.ndarray,
        expression: np.ndarray,
        cell_types: np.ndarray,
        comm_matrices: dict
    ) -> sc.AnnData:
        """Create AnnData object"""
        gene_names = ['Ligand_A', 'Ligand_B', 'Ligand_C', 'Ligand_D',
                      'Receptor_A', 'Receptor_B', 'Receptor_C', 'Receptor_D']
        gene_names += [f'Gene_{i}' for i in range(self.n_background_genes)]

        adata = sc.AnnData(X=expression)
        adata.var_names = gene_names
        adata.obs_names = [f'Cell_{i}' for i in range(len(coords))]

        adata.obs['cell_type'] = [f'CellType_{ct}' for ct in cell_types]
        adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

        adata.obsm['spatial'] = coords

        adata.uns['ground_truth_communication'] = comm_matrices
        adata.uns['simulation_params'] = {
            'letter': self.letter,
            'n_cells': self.n_cells,
            'n_cell_types': self.n_cell_types,
            'scale': self.scale,
            'expression_strength': self.expression_strength,
            'heterogeneity': self.heterogeneity,
            'random_seed': self.random_seed
        }

        return adata

    def generate_dataset(self, save_path: str = None) -> sc.AnnData:
        """Generate complete dataset"""
        print(f"Generating {self.letter}-shaped dataset...")

        # Generate letter shape
        print("  Creating letter shape...")
        coords = self.generate_letter_shape(self.letter, self.n_cells)

        # Add noise
        print("  Adding spatial noise...")
        coords = self.add_noise_to_coords(coords, noise_level=2.0)

        # Assign cell types
        print("  Assigning cell types...")
        cell_types = self.assign_cell_types_to_shape(coords)

        # Generate expression
        print("  Generating expression...")
        expression = self.generate_expression(cell_types)

        # Generate communication matrix
        print("  Generating communication matrices...")
        comm_matrices = self.generate_ground_truth_communication(coords, expression)

        # CreatingAnnData
        print("  Creating AnnData...")
        adata = self.create_anndata(coords, expression, cell_types, comm_matrices)

        if save_path:
            print(f"  Saving to {save_path}...")
            adata.write_h5ad(save_path)

        print(f"  Done! {adata.n_obs} cells, {adata.n_vars} genes")
        print(f"  Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")

        return adata


if __name__ == '__main__':
    print("="*80)
    print(" "*25 + "UMN Letter Shape Datasets")
    print("="*80)

    letters = ['U', 'M', 'N', 'UMN']

    for letter in letters:
        print(f"\n{'='*80}")
        print(f"Creating {letter}-shaped dataset")
        print('='*80)

        simulator = LetterShapeSimulator(
            letter=letter,
            n_cells=800,
            n_cell_types=4,
            scale=80.0,
            expression_strength=12.0,
            heterogeneity=0.15,
            random_seed=42
        )

        adata = simulator.generate_dataset(save_path=f'Letter_{letter}.h5ad')

        # Quick visualization
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        coords = adata.obsm['spatial']

        for ct in sorted(adata.obs['cell_type'].unique()):
            mask = adata.obs['cell_type'] == ct
            ax.scatter(coords[mask, 0], coords[mask, 1], label=ct, s=20, alpha=0.7)

        ax.set_title(f'{letter}-shaped Dataset', fontsize=16, fontweight='bold')
        ax.legend()
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_aspect('equal')
        ax.grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(f'Letter_{letter}_shape.png', dpi=300, bbox_inches='tight')
        plt.close()

        print(f"  Saved visualization: Letter_{letter}_shape.png")

    print("\n" + "="*80)
    print("All UMN letter datasets created!")
    print("="*80)
