# Input and script provenance

The archived final pipeline reads 15 normalized inputs. The earlier editorial note stated 16, but direct inspection of the executable `input_paths` object and its generated manifest confirmed 15 named runtime inputs; no sixteenth input is referenced by the archived script.

| Input group | Upstream workflow |
|---|---|
| Species universe and pruned WCVP tree | `scripts/upstream/01_prepare_species_universe.R` |
| GBIF and GHS-POP species summaries | `scripts/upstream/02_process_gbif_ghs.R` |
| Chemical incidence and architecture tables | `scripts/upstream/03_build_chemical_architecture.R` and `04_regenerate_chemical_architecture_figures.R` |
| Research-selection model tables | `scripts/upstream/06_model_research_selection.R` |
| Genus-level phylogenetic input | `scripts/upstream/07_prepare_genus_phylogeny.R` |

The exact Stage 05B script has MD5 `9f38d16e009b848dcad793ce250d8ed0`. The deposited diagnostic package includes the normalized 2,582-taxon model matrix, the fitted model object, model-fit tables, convergence checks, sensitivity analyses, predictor correlations, provenance table, execution log and session information. Regenerable Stage 05B figures and duplicate output tables are intentionally omitted.

The exact final pipeline archived on 14 August 2026 has MD5 `00bfb4ae76ada6383e04ed0493368548`. The publication script refuses to proceed if the local final script differs from this checksum.
