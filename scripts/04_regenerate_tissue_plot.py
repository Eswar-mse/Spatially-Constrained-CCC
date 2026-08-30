import scanpy as sc
import os

print("Regenerating tissue-overlay cluster plot...")

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
OUTPUTS_DIR = os.path.join(PROJECT_ROOT, 'outputs')
os.makedirs(OUTPUTS_DIR, exist_ok=True)

# Load the clustered adata saved by 02_spatial_clustering.py
input_path = os.path.join(DATA_DIR, '02_clustered_adata.h5ad')
print(f"Loading {input_path}...")
adata = sc.read_h5ad(input_path)

# Plot clusters overlaid on the H&E tissue image
out_path = os.path.join(OUTPUTS_DIR, 'tissue_clusters_overlay.png')
sc.pl.spatial(
    adata,
    color='clusters',
    title='Mouse Brain (Sagittal) — Spatial Domains',
    save=None,
    show=False,
)
import matplotlib.pyplot as plt
plt.savefig(out_path, dpi=300, bbox_inches='tight')
print(f"Saved: {out_path}")
