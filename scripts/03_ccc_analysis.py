import pandas as pd
import numpy as np
import liana as li
import scanpy as sc
import anndata as ad
import squidpy as sq
import os
import json

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

    # magnitude_rank is a rank aggregation score, not a p value.
    sig_bulk = bulk_res[bulk_res['magnitude_rank'] <= 0.05].copy()
    print(f"\nCandidates below magnitude-rank cutoff: {len(sig_bulk)}")

    # 3. Calculate Spatial Adjacency (The Reality Check)
    print("Calculating physical spatial adjacency between clusters...")
    sq.gr.nhood_enrichment(adata, cluster_key="clusters", n_perms=1000, seed=42)

    cluster_names = adata.obs['clusters'].cat.categories
    enrichment_df = pd.DataFrame(
        adata.uns['clusters_nhood_enrichment']['zscore'], 
        index=cluster_names, 
        columns=cluster_names
    )

    # 4. Apply the Spatial Constraint
    # This is an exploratory cluster-pair enrichment screen. It is not a
    # ligand-receptor-specific test and does not include multiple-test correction.
    ZSCORE_SIG_THRESHOLD = 1.96

    def check_adjacency(row):
        source = row['source']
        target = row['target']
        return enrichment_df.loc[source, target] > ZSCORE_SIG_THRESHOLD

    sig_bulk['physically_adjacent'] = sig_bulk.apply(check_adjacency, axis=1)
    spatial_filtered_res = sig_bulk[sig_bulk['physically_adjacent']].copy()

    excluded = len(sig_bulk) - len(spatial_filtered_res)
    retention_rate = len(spatial_filtered_res) / len(sig_bulk) * 100 if len(sig_bulk) else 0

    print(f"\nCluster adjacency screen: z-score > {ZSCORE_SIG_THRESHOLD}")
    print(f"Interactions remaining after enforcing physical proximity: {len(spatial_filtered_res)}")
    print(f"Candidates excluded by adjacency screen: {excluded}")
    print(f"Candidate retention rate: {retention_rate:.2f}%")

    # Preserve the complete result for inspection and later visualizations.
    sig_bulk.to_csv(os.path.join(OUTPUTS_DIR, "all_ranked_candidates.csv"), index=False)
    spatial_filtered_res.to_csv(os.path.join(OUTPUTS_DIR, "retained_candidates.csv"), index=False)
    enrichment_df.to_csv(os.path.join(OUTPUTS_DIR, "cluster_adjacency_zscores.csv"))
    summary = {
        "dataset": "10x Visium mouse brain sagittal section (Squidpy visium_hne_adata)",
        "candidate_count": len(sig_bulk),
        "retained_count": len(spatial_filtered_res),
        "excluded_count": excluded,
        "retention_percent": round(retention_rate, 2),
        "magnitude_rank_cutoff": 0.05,
        "adjacency_zscore_cutoff": ZSCORE_SIG_THRESHOLD,
        "neighborhood_permutations": 1000,
        "random_seed": 42,
        "interpretation": "Cluster adjacency screen; not a biological false-positive estimate",
    }
    with open(os.path.join(OUTPUTS_DIR, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    # 5. Generate and Save Final Visualization
    print("\nGenerating exploratory dot plot...")
    border_interactions = spatial_filtered_res[spatial_filtered_res['source'] != spatial_filtered_res['target']].copy()
    top_hits = border_interactions.sort_values(by='magnitude_rank').head(15)

    if top_hits.empty:
        print("No cross-cluster candidates passed the screen; skipping dot plot.")
        return

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

if __name__ == '__main__':
    main()
