# Spatially Constrained Cell Communication

An exploratory analysis of ligand–receptor candidates in a **10x Visium mouse brain sagittal section**. The pipeline compares expression based candidates with the spatial neighborhood enrichment of their source and target clusters.

**[Explore the interactive case study](https://spatial-ccc-eswar.archeddegu7.chatgpt.site)** · [View the interactive tissue map](assets/interactive_clusters.html)

![Mouse brain tissue with Leiden clusters](assets/tissue_clusters_overlay.png)

## Question and result

Can tissue position help prioritize cell communication candidates inferred from expression alone? In the original exploratory run, LIANA+ produced **7,659** candidates below a magnitude rank cutoff of 0.05. **876** were retained after requiring source and target clusters to have Squidpy neighborhood enrichment **z > 1.96**; **6,783** were excluded. The retention rate was **11.44%**.

These figures are from the original run and have not been independently reproduced as part of the portfolio update. An excluded candidate is **not a proven false positive**. Neighborhood enrichment is a cluster pair statistic, not evidence for or against an individual ligand–receptor event. In particular, this screen can exclude biologically real longer range signaling.

![Top retained ligand receptor candidates](assets/spatial_ccc_border_dotplot.png)

## Pipeline

1. `scripts/01_data_prep.py` downloads Squidpy's `visium_hne_adata` mouse brain example, filters spots with fewer than 500 counts and genes seen in fewer than 10 spots, then normalizes and log transforms expression.
2. `scripts/02_spatial_clustering.py` runs PCA, a transcriptomic neighbor graph and Leiden clustering at resolution 0.5. It separately builds a **physical** spatial neighbor graph and exports an interactive tissue map.
3. `scripts/03_ccc_analysis.py` uses LIANA+ rank aggregation with the `mouseconsensus` resource, selects candidates with `magnitude_rank <= 0.05`, computes neighborhood enrichment with 1,000 permutations and seed 42, and applies the cluster pair z score screen. It writes full and retained candidate tables, an adjacency matrix, a JSON summary, and a dot plot.
4. `scripts/04_regenerate_tissue_plot.py` renders clusters on the tissue image.

The interactive case study source is in [`docs/`](docs/), ready for static hosting from the repository. Its assets are independent of the compute pipeline.

## Run locally

Use **Python 3.12**. The pinned analysis packages are in [`requirements.txt`](requirements.txt).

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/01_data_prep.py
python scripts/02_spatial_clustering.py
python scripts/03_ccc_analysis.py
python scripts/04_regenerate_tissue_plot.py
```

The first step downloads a large `.h5ad` dataset. Generated data remains in `data/`; analysis exports are written to `outputs/`. The repo includes selected example figures in `assets/`. If you rerun with a different library version, dataset revision, or clustering outcome, expect the counts to change.

## Interpretation limits

- This is one tissue section with no biological replicate.
- Visium spots may contain multiple cells. Numeric clusters are transcriptomic domains, not annotated cell types.
- Ligand and receptor expression plus proximity does not prove an active signaling event.
- The z score threshold is an exploratory screen based on a normal approximation. It does not correct for multiple comparisons across cluster pairs.
- A cluster pair enrichment score applies equally to every candidate assigned to that pair. It cannot distinguish which individual ligand–receptor interactions occur.
- Some signaling acts beyond immediate neighboring spots.

## Next research steps

Replicate across tissue sections, annotate cell populations, apply multiple testing correction, and compare spatial screening with ligand–receptor specific spatial methods and experimental validation.

This is a portfolio case study and a reproducible starting point for further analysis, not a clinical or publication grade finding.
