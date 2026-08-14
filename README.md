# Arecaceae chemical research landscape

Reproducible data and code supporting a curated evidence synthesis of the taxonomic, phylogenetic, geographic, functional, socioeconomic and analytical structure of phytochemical research in Arecaceae.

## Reproducible archive

The final manuscript pipeline is `scripts/Arecaceae_single_pipeline_end.R`. A direct audit of that archived script identifies **15 named normalized inputs**. Their repository paths, file sizes and checksums are recorded in `results/manifests/`.

The archive retains the normalized inputs and selected diagnostics required to audit the final analyses. Superseded exploratory and upstream-generation scripts are intentionally excluded. Regenerable figures, duplicate tables, manuscripts, article files and raw external datasets are also excluded.

## Run

Restore the archived environment and execute the final pipeline from the repository root:

```r
renv::restore()
source("scripts/Arecaceae_single_pipeline_end.R")
```

The exact R version and platform from the final archived run are recorded in `results/logs/final_pipeline_sessionInfo.txt`.

## Version and citation

- Stable repository: <https://github.com/carloscarollo-UFMS/arecaceae-chemical-research-landscape>
- Versioned archive: <https://github.com/carloscarollo-UFMS/arecaceae-chemical-research-landscape/tree/v1.0.0>
- Code licence: MIT
- Data licence: CC BY 4.0

To cite this repository, use the metadata in `CITATION.cff` and report the version tag or exact commit used in the analysis.

Copyrighted article full texts, raw publisher files and credentials are not included.
