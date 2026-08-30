import pandas as pd
import numpy as np
import liana as li
import scanpy as sc
import anndata as ad
import squidpy as sq
import os

def main():
    print("--- STEP 3: SPATIAL CELL-CELL COMMUNICATION ---")

    PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
    OUTPUTS_DIR = os.path.join(PROJECT_ROOT, 'outputs')
    os.makedirs(OUTPUTS_DIR, exist_ok=True)

    # 1. Load Data
    input_path = os.path.join(DATA_DIR, '02_clustered_adata.h5ad')
    print(f"Loading data from {input_path}...")
    adata = sc.read_h5ad(input_path)

    # 2. Run Global CCC (The "Soup" Model)
    print("Running standard cell-cell communication (ignoring space)...")
    li.mt.rank_aggregate(adata, groupby='clusters', resource_name='mouseconsensus', expr_prop=0.1, use_raw=False, verbose=False)
    bulk_res = adata.uns['liana_res'].copy()

    # Filter for the most robust, significant interactions
    sig_bulk = bulk_res[bulk_res['magnitude_rank'] <= 0.05].copy()
    print(f"\nTotal significant interactions predicted globally: {len(sig_bulk)}")

    # 3. Calculate Spatial Adjacency (The Reality Check)
    print("Calculating physical spatial adjacency between clusters...")
    sq.gr.nhood_enrichment(adata, cluster_key="clusters")

    cluster_names = adata.obs['clusters'].cat.categories
    enrichment_df = pd.DataFrame(
        adata.uns['clusters_nhood_enrichment']['zscore'], 
        index=cluster_names, 
        columns=cluster_names
    )

    # 4. Apply the Spatial Constraint
    # Squidpy's nhood_enrichment z-score comes from a permutation test against a null
    # distribution of shuffled cluster labels. Under the normal approximation, z > 1.96
    # corresponds to p < 0.05 (two-tailed) -- i.e. "significantly more adjacent than
    # expected by chance," not merely "positive." A bare >0 check lets weakly/noisily
    # positive pairs slip through and inflates how many interactions look validated.
    ZSCORE_SIG_THRESHOLD = 1.96  # p < 0.05, two-tailed

    def check_adjacency(row):
        source = row['source']
        target = row['target']
        return enrichment_df.loc[source, target] > ZSCORE_SIG_THRESHOLD

    sig_bulk['physically_adjacent'] = sig_bulk.apply(check_adjacency, axis=1)
    spatial_filtered_res = sig_bulk[sig_bulk['physically_adjacent'] == True]

    false_positives = len(sig_bulk) - len(spatial_filtered_res)
    fp_rate = (false_positives / len(sig_bulk)) * 100 if len(sig_bulk) > 0 else 0

    print(f"\nAdjacency significance threshold: z-score > {ZSCORE_SIG_THRESHOLD} (p < 0.05, two-tailed)")
    print(f"Interactions remaining after enforcing physical proximity: {len(spatial_filtered_res)}")
    print(f"False-Positive Interactions eliminated: {false_positives}")
    print(f"False-Positive Rate: {fp_rate:.2f}%")

    # 5. Generate and Save Final Visualization
    print("\nGenerating publication-ready dot plot...")
    border_interactions = spatial_filtered_res[spatial_filtered_res['source'] != spatial_filtered_res['target']].copy()
    top_hits = border_interactions.sort_values(by='magnitude_rank').head(15)

    g = li.pl.dotplot(
        liana_res=top_hits,
        colour='magnitude_rank',
        size='magnitude_rank',
        inverse_colour=True,
        inverse_size=True
    )

    g.save(os.path.join(OUTPUTS_DIR, "spatial_ccc_border_dotplot.pdf"), width=10, height=8, dpi=300, verbose=False)
    g.save(os.path.join(OUTPUTS_DIR, "spatial_ccc_border_dotplot.png"), width=10, height=8, dpi=300, verbose=False)

    print("Success! Final pipeline complete. Plots saved to 'outputs' folder.")

# THIS is the magic line that makes it safe for terminal execution!
if __name__ == '__main__':
    main()
