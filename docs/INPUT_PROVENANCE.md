# Input and analysis provenance

The archived final pipeline reads 15 normalized inputs. The earlier editorial note stated 16, but direct inspection of the executable `input_paths` object and its generated manifest confirmed 15 named runtime inputs; no sixteenth input is referenced by the archived script.

| Input group | Archived repository inputs |
|---|---|
| Species universe and pruned WCVP tree | `data/reference/01_species_tree_WCVP_pruned.tre` and `data/interim/02_species_analysis_with_GBIF_GHS.csv` |
| Chemical incidence and study corpus | `data/interim/04C_incidence_article_taxon_organ_entity.csv` and `data/interim/04C_articles_250.csv` |
| Chemical architecture | `data/derived/Table_04C_*.csv` |
| Research-selection analyses | `data/derived/Table_05B_*.csv` and selected files under `diagnostics/stage05B/` |
| Phylogenetic reference | `data/reference/Species_tree_all_genes_rooted_original.tre` |

The normalized files in `data/` are the authoritative inputs for the archived final pipeline. Superseded exploratory and upstream-generation scripts are intentionally excluded. The selected Stage 05B diagnostic package includes the normalized 2,582-taxon model matrix, fitted model object, model-fit tables, convergence checks, sensitivity analyses, predictor correlations, provenance table, execution log and session information. Regenerable Stage 05B figures and duplicate output tables are intentionally omitted.

The exact final pipeline archived on 14 August 2026 has MD5 `00bfb4ae76ada6383e04ed0493368548` and SHA-256 `4ee5b9acdc2687574f48ec8535ecef4f3d0379effc52f16c708af857a9b5edb7`.
