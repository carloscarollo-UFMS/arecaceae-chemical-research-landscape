# Arecaceae chemical research landscape

Reproducible data and code supporting a curated evidence synthesis of the taxonomic, phylogenetic, geographic, functional, socioeconomic and analytical structure of phytochemical research in Arecaceae.

## Reproducible archive

The final manuscript pipeline is `scripts/Arecaceae_single_pipeline_end.R`. A direct audit of that archived script identifies **15 named normalized inputs**. Their repository paths, file sizes and checksums are recorded in `results/manifests/`.

The upstream selection-model workflow is deposited as `scripts/upstream/06_model_research_selection.R`. Its normalized analytical matrix, fitted model object, fit statistics, convergence checks, sensitivity analyses, execution log and session information are retained under `diagnostics/stage05B/`. Regenerable figures, duplicate tables, manuscripts, article files and raw external datasets are excluded.

## Run

Restore the archived environment and execute the final pipeline from the repository root:

```r
renv::restore()
source("scripts/Arecaceae_single_pipeline_end.R")
```

The exact R version and platform from the final archived run are recorded in `results/logs/final_pipeline_sessionInfo.txt`.

## Release and citation

- Stable repository: <https://github.com/carloscarollo-UFMS/arecaceae-chemical-research-landscape>
- Immutable release: `v1.0.0`
- Code licence: MIT
- Data licence: CC BY 4.0

The repository DOI will be added to `CITATION.cff` after the GitHub release is archived in Zenodo. Until then, cite the immutable GitHub release URL.

Copyrighted article full texts, raw publisher files and credentials are not included.
