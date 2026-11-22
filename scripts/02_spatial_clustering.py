import scanpy as sc
import squidpy as sq
import pandas as pd
import plotly.express as px
import os

print("--- STEP 2: CLUSTERING & INTERACTIVE VISUALIZATION ---")

# Anchor all paths to the project root (parent of this script's scripts/ folder),
# so this runs the same whether launched from root, scripts/, or anywhere else.
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
OUTPUTS_DIR = os.path.join(PROJECT_ROOT, 'outputs')
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUTPUTS_DIR, exist_ok=True)

# 1. Load the preprocessed data
input_path = os.path.join(DATA_DIR, '01_processed_adata.h5ad')
print(f"Loading data from {input_path}...")
adata = sc.read_h5ad(input_path)

# 2. Dimensionality Reduction & Clustering
print("Running PCA and Leiden clustering...")
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added="clusters", resolution=0.5) 

# 3. Spatial Graph Construction (The Physical Mesh)
print("Calculating physical spatial neighborhoods...")
sq.gr.spatial_neighbors(adata)

# 4. Extract Data for Plotly
print("Generating interactive Plotly tissue map...")
spatial_coords = adata.obsm['spatial']
plot_df = pd.DataFrame({
    'X': spatial_coords[:, 0],
    'Y': -spatial_coords[:, 1], # Invert Y for correct orientation
    'Cluster': adata.obs['clusters']
})

# 5. Build the Interactive Visual
fig = px.scatter(
    plot_df, 
    x='X', 
    y='Y', 
    color='Cluster',
    title='Interactive Spatial Domains of Mouse Brain Tissue',
    color_discrete_sequence=px.colors.qualitative.Bold,
    template='plotly_dark'
)

fig.update_traces(marker=dict(size=7, opacity=0.9, line=dict(width=0.5, color='DarkSlateGrey')))
fig.update_layout(xaxis_title="Spatial X", yaxis_title="Spatial Y")

# Save as an interactive HTML file
html_out = os.path.join(OUTPUTS_DIR, 'interactive_clusters.html')
fig.write_html(html_out)
print(f"Stunning interactive visual saved to: {html_out}")

# 6. Save the state for Script 3
output_path = os.path.join(DATA_DIR, '02_clustered_adata.h5ad')
adata.write(output_path)
print(f"Clustered data safely stored at: {output_path}")
