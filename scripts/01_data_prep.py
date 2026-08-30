import scanpy as sc
import squidpy as sq
import os

print("Step 1: Data preparation and QC")

# Anchor all paths to the project root (parent of this script's scripts/ folder),
# so this runs the same whether launched from root, scripts/, or anywhere else.
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
OUTPUTS_DIR = os.path.join(PROJECT_ROOT, 'outputs')

# 1. Ensure that the directories exist
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUTPUTS_DIR, exist_ok=True)

# 2. Download the spatial dataset
print("Downloading 10x Visium Mouse Brain (Sagittal) dataset")
adata= sq.datasets.visium_hne_adata()
adata.var_names_make_unique()

# 3. Quality control
print("Running QC")
sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_genes(adata, min_cells=10)

# 4. Normalization
print("Normalizing and log-transforming")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# 5. Save the stat to the data folder
output_path = os.path.join(DATA_DIR, "01_processed_adata.h5ad")
adata.write(output_path)

print(f"Sucess, preprocessed data safely stored at {output_path}")


