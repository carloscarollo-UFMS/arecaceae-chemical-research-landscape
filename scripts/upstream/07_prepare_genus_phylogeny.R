# =============================================================================
# Arecaceae phylogeny with 184 genera and annotation datasets for iTOL
# Includes genus richness to make the denominator of study coverage explicit.
#
# Phylogenetic source:
# Bellot et al. (2026). Phylogenomics and a New Fossil Synthesis Illuminate
# the Early Evolution of Palms (Arecaceae). Systematic Biology, syag022.
# DOI: 10.1093/sysbio/syag022
# Data repository: https://doi.org/10.5061/dryad.pzgmsbcwg
#
# Instructions:
# 1. Open this file in RStudio.
# 2. Install the packages listed below, if necessary.
# 3. Run Source > Source or execute: Rscript generate_itol_files.R
# 4. Results will be written to the R_results_with_genus_richness directory.
# =============================================================================

required_packages <- c("ape", "readr", "dplyr", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the missing packages before running this script:\n",
    "install.packages(c(",
    paste(sprintf('\"%s\"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(ape)
  library(readr)
  library(dplyr)
  library(tidyr)
})

get_script_directory <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  script_file <- grep("^--file=", arguments, value = TRUE)
  if (length(script_file) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", script_file))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(path)) return(dirname(normalizePath(path)))
  }
  normalizePath(getwd())
}

project_directory <- get_script_directory()
input_directory <- file.path(project_directory, "input")
output_directory <- file.path(project_directory, "R_results_with_genus_richness")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

resolve_input_file <- function(candidates) {
  search_directories <- unique(c(input_directory, project_directory))
  possible_paths <- unlist(
    lapply(search_directories, function(directory) file.path(directory, candidates)),
    use.names = FALSE
  )
  existing_paths <- possible_paths[file.exists(possible_paths)]
  if (length(existing_paths) == 0) return(NA_character_)
  normalizePath(existing_paths[[1]])
}

tree_file <- resolve_input_file(c(
  "Species_tree_all_genes_rooted_original.tre"
))
powo_file <- resolve_input_file(c(
  "powo_native_distribution.csv",
  "10_distribuicao_nativa_atual_nao_duvidosa.csv",
  "10_distribuicao_nativa_atual_nao_duvidosa(2).csv"
))
study_file <- resolve_input_file(c(
  "studied_arecaceae_species.csv",
  "especies_Arecaceae_com_estudos_no_corpus.csv"
))

missing_inputs <- c(
  "phylogenetic tree" = tree_file,
  "POWO distribution table" = powo_file,
  "studied-species table" = study_file
)
missing_inputs <- names(missing_inputs)[is.na(missing_inputs)]
if (length(missing_inputs) > 0) {
  stop(
    "Missing input file(s): ", paste(missing_inputs, collapse = ", "), ".\n",
    "Extract the complete ZIP file and preserve the input subdirectory.\n",
    "Expected project directory: ", project_directory, "\n",
    "Expected input directory: ", input_directory
  )
}

# -----------------------------------------------------------------------------
# 1. Published tree: clean, prune, and reduce to genus level
# -----------------------------------------------------------------------------

tree <- ape::read.tree(tree_file)

# Three samples are masked in the deposited tree. Their identities were
# recovered from Supplementary Table S1 and confirmed in Figure 1.
sample_replacements <- c(
  Genus_species4 = "Truongsonia_lecongkietii",
  Genus_species7 = "Areca_chaiana",
  Genus_species8 = "Oncosperma_fasciculatum"
)
replacement_index <- match(names(sample_replacements), tree$tip.label)
tree$tip.label[replacement_index[!is.na(replacement_index)]] <-
  sample_replacements[!is.na(replacement_index)]

# Remove the outgroup and one additional sample from each duplicated genus.
# Areca catechu and Oncosperma tigillarium are retained as representatives.
tips_to_remove <- intersect(
  c("Dasypogon_bromeliifolius", "Areca_chaiana", "Oncosperma_fasciculatum"),
  tree$tip.label
)
tree <- ape::drop.tip(tree, tips_to_remove, collapse.singles = TRUE)

# Convert representative species labels to genus labels and align spelling
# with the supplied POWO extract.
tree$tip.label <- sub("_.*$", "", tree$tip.label)
tree$tip.label[tree$tip.label == "Acoelorrhaphe"] <- "Acoelorraphe"
tree$node.label <- NULL

if (ape::Ntip(tree) != 184 || anyDuplicated(tree$tip.label)) {
  stop("The final tree should contain exactly 184 unique genera.")
}

tree_with_astral_branch_lengths <- tree
ape::write.tree(
  tree_with_astral_branch_lengths,
  file.path(output_directory, "S3_arecaceae_184_genera_astral_branch_lengths.nwk")
)
tree$edge.length <- NULL
ape::write.tree(tree, file.path(output_directory, "01_arecaceae_184_genera_cladogram.nwk"))

# -----------------------------------------------------------------------------
# 2. POWO: species richness per genus and native distribution
# -----------------------------------------------------------------------------

powo <- readr::read_delim(
  powo_file,
  delim = ";",
  locale = locale(encoding = "UTF-8"),
  show_col_types = FALSE,
  trim_ws = TRUE
)

# Backward compatibility with the column names in the original supplied file.
if ("taxon_hibrido" %in% names(powo) && !("hybrid_taxon" %in% names(powo))) {
  powo <- powo %>% rename(hybrid_taxon = taxon_hibrido)
}
if (!("hybrid_taxon" %in% names(powo))) {
  stop("The POWO table does not contain the required hybrid_taxon column.")
}

powo_nonhybrids <- powo %>%
  mutate(hybrid_taxon = toupper(as.character(hybrid_taxon))) %>%
  filter(hybrid_taxon != "TRUE")

powo_species <- powo_nonhybrids %>%
  distinct(plant_name_id, .keep_all = TRUE)

powo_genera <- sort(unique(powo_species$genus))
if (!setequal(tree$tip.label, powo_genera)) {
  stop(
    "Tree and POWO genera do not match.\n",
    "Tree only: ", paste(setdiff(tree$tip.label, powo_genera), collapse = ", "), "\n",
    "POWO only: ", paste(setdiff(powo_genera, tree$tip.label), collapse = ", ")
  )
}

richness <- powo_species %>%
  count(genus, name = "powo_species")

# -----------------------------------------------------------------------------
# 3. Literature corpus reconciliation and study coverage
# -----------------------------------------------------------------------------

studies <- readr::read_csv(
  study_file,
  locale = locale(encoding = "UTF-8"),
  show_col_types = FALSE
)

# Backward compatibility with the column name in the original corpus table.
if ("especie_consolidada" %in% names(studies) &&
    !("consolidated_species" %in% names(studies))) {
  studies <- studies %>% rename(consolidated_species = especie_consolidada)
}
if (!("consolidated_species" %in% names(studies))) {
  stop("The literature table does not contain the required consolidated_species column.")
}

# Taxonomic mappings used in the analysis. NA indicates that no defensible
# accepted-name match was found in the supplied POWO extract.
name_mapping <- tibble::tribble(
  ~corpus_name, ~accepted_name,
  "Arenga micranta", "Arenga micrantha",
  "Astrocaryum aculeatum", "Astrocaryum tucuma",
  "Calamus draco", NA_character_,
  "Calamus quiquesetinervius", "Calamus formosanus",
  "Dypsis decaryi", "Chrysalidocarpus decaryi",
  "Dypsis leptocheilos", "Chrysalidocarpus leptocheilos",
  "Dypsis lutescens", "Chrysalidocarpus lutescens",
  "Geonoma litoralis", "Geonoma pohliana",
  "Jubaeopsis caffra", "Jubaeopsis afra",
  "Livistona decipiens", "Livistona decora",
  "Phoenix hanceana", "Phoenix loureiroi",
  "Phoenix humilis", "Phoenix loureiroi",
  "Phoenix pusilla", "Phoenix sylvestris",
  "Pinanga kuhlii", "Pinanga coronata",
  "Plectocomia kerrana", NA_character_,
  "Plectocomia khasiyana", "Plectocomia elongata",
  "Ptychosperma macarthurii", "Ptychosperma propinquum",
  "Sabal blackburniana", "Sabal palmetto",
  "Trachycarpus wagnerianus", "Trachycarpus fortunei",
  "Trithrinax acanthocoma", "Trithrinax brasiliensis",
  "Washingtonia robusta", "Washingtonia filifera"
)

powo_names <- powo_species$taxon_name

reconciled_studies <- studies %>%
  transmute(corpus_name = consolidated_species) %>%
  left_join(name_mapping, by = "corpus_name") %>%
  mutate(
    direct_match = corpus_name %in% powo_names,
    accepted_name = if_else(direct_match, corpus_name, accepted_name)
  ) %>%
  filter(!is.na(accepted_name), accepted_name %in% powo_names) %>%
  distinct(accepted_name)

if (nrow(reconciled_studies) != 153) {
  stop(
    "Reconciliation should produce 153 accepted studied species; obtained: ",
    nrow(reconciled_studies)
  )
}

coverage <- reconciled_studies %>%
  left_join(
    powo_species %>% select(taxon_name, genus),
    by = c("accepted_name" = "taxon_name")
  ) %>%
  count(genus, name = "studied_species")

# -----------------------------------------------------------------------------
# 4. Native ranges: six grouped continents and eight TDWG/POWO regions
# -----------------------------------------------------------------------------

tdwg_levels <- c(
  "AFRICA", "ASIA-TEMPERATE", "ASIA-TROPICAL", "AUSTRALASIA",
  "EUROPE", "NORTHERN AMERICA", "PACIFIC", "SOUTHERN AMERICA"
)

tdwg_labels <- c(
  "Africa", "Temperate Asia", "Tropical Asia", "Australasia",
  "Europe", "Northern America", "Pacific", "Southern America"
)

continents <- c("Africa", "Asia", "Europe", "North America", "South America", "Oceania")

genus_distribution <- powo_nonhybrids %>%
  distinct(genus, continent, region) %>%
  mutate(
    grouped_continent = case_when(
      continent == "AFRICA" ~ "Africa",
      continent %in% c("ASIA-TEMPERATE", "ASIA-TROPICAL") ~ "Asia",
      continent == "EUROPE" ~ "Europe",
      continent == "NORTHERN AMERICA" ~ "North America",
      continent == "SOUTHERN AMERICA" &
        region %in% c("Caribbean", "Central America") ~ "North America",
      continent == "SOUTHERN AMERICA" &
        region %in% c(
          "Northern South America", "Western South America",
          "Southern South America", "Brazil"
        ) ~ "South America",
      continent %in% c("AUSTRALASIA", "PACIFIC") ~ "Oceania",
      TRUE ~ NA_character_
    )
  )

continent_wide <- genus_distribution %>%
  filter(!is.na(grouped_continent)) %>%
  distinct(genus, grouped_continent) %>%
  mutate(value = 1L) %>%
  complete(genus = powo_genera, grouped_continent = continents, fill = list(value = 0L)) %>%
  mutate(grouped_continent = gsub(" ", "_", grouped_continent)) %>%
  pivot_wider(names_from = grouped_continent, values_from = value, names_prefix = "continent_")

tdwg_wide <- genus_distribution %>%
  filter(continent %in% tdwg_levels) %>%
  distinct(genus, continent) %>%
  mutate(value = 1L) %>%
  complete(genus = powo_genera, continent = tdwg_levels, fill = list(value = 0L)) %>%
  pivot_wider(names_from = continent, values_from = value, names_prefix = "tdwg_")

data <- tibble(genus = powo_genera) %>%
  left_join(richness, by = "genus") %>%
  left_join(coverage, by = "genus") %>%
  mutate(
    studied_species = coalesce(studied_species, 0L),
    study_coverage_percent = 100 * studied_species / powo_species,
    has_studies = as.integer(studied_species > 0)
  ) %>%
  left_join(continent_wide, by = "genus") %>%
  left_join(tdwg_wide, by = "genus")

continent_columns <- paste0("continent_", gsub(" ", "_", continents))
data$native_continents <- apply(
  data[, continent_columns],
  1,
  function(values) paste(continents[as.logical(values)], collapse = "; ")
)
data <- data %>% relocate(native_continents, .after = has_studies)

readr::write_excel_csv(data, file.path(output_directory, "06_genus_level_data.csv"))

# -----------------------------------------------------------------------------
# 5. Write iTOL annotation datasets
# -----------------------------------------------------------------------------

write_itol <- function(lines, filename) {
  writeLines(enc2utf8(lines), file.path(output_directory, filename), useBytes = TRUE)
}

# Genus richness is displayed as a neutral gray bar on a logarithmic scale.
# The values remain the actual numbers of accepted non-hybrid POWO species;
# LOG_SCALE affects only their graphical display in iTOL.
richness_lines <- c(
  "DATASET_SIMPLEBAR",
  "SEPARATOR TAB",
  "DATASET_LABEL\tGenus richness (accepted POWO species; log scale)",
  "COLOR\t#7A7A7A",
  paste("DATASET_SCALE", 1, 10, 100, max(data$powo_species), sep = "\t"),
  "WIDTH\t80",
  "MARGIN\t5",
  "LOG_SCALE\t1",
  "HEIGHT_FACTOR\t0.65",
  "SHOW_VALUE\t0",
  "SHOW_LABELS\t0",
  "BORDER_WIDTH\t0",
  "DATA",
  sprintf("%s\t%d", data$genus, data$powo_species)
)
write_itol(richness_lines, "02_itol_genus_richness_logbar.txt")

coverage_lines <- c(
  "DATASET_GRADIENT",
  "SEPARATOR TAB",
  "DATASET_LABEL\tStudy coverage (% of accepted species)",
  "COLOR\t#2166AC",
  "STRIP_WIDTH\t35",
  "MARGIN\t5",
  "BORDER_WIDTH\t0",
  "AUTO_LEGEND\t0",
  "COLOR_MIN\t#F7FBFF",
  "USE_MID_COLOR\t1",
  "COLOR_MID\t#6BAED6",
  "COLOR_MAX\t#08306B",
  "SHOW_LABELS\t0",
  "LEGEND_TITLE\tStudy coverage (% of accepted species)",
  "LEGEND_SCALE\t0.8",
  "LEGEND_HORIZONTAL\t0",
  "LEGEND_SHAPES\t1\t1\t1\t1\t1",
  "LEGEND_COLORS\t#F7FBFF\t#C6DBEF\t#6BAED6\t#2171B5\t#08306B",
  "LEGEND_LABELS\t0\t25\t50\t75\t100",
  "DATA",
  sprintf("%s\t%.6f", data$genus, data$study_coverage_percent)
)
write_itol(coverage_lines, "03_itol_study_coverage_percent.txt")

study_presence_lines <- c(
  "DATASET_BINARY",
  "SEPARATOR TAB",
  "DATASET_LABEL\tStudy presence",
  "COLOR\t#238B45",
  "FIELD_SHAPES\t2",
  "FIELD_LABELS\tAt least one studied species",
  "FIELD_COLORS\t#238B45",
  "MARGIN\t5",
  "HEIGHT_FACTOR\t0.8",
  "SHOW_LABELS\t1",
  "DATA",
  sprintf("%s\t%d", data$genus, data$has_studies)
)
write_itol(study_presence_lines, "S1_itol_study_presence.txt")

continent_colors <- c("#D55E00", "#E69F00", "#0072B2", "#CC79A7", "#009E73", "#56B4E9")
continent_matrix <- as.matrix(data[, continent_columns])
continent_matrix[continent_matrix == 0] <- -1

continent_lines <- c(
  "DATASET_BINARY",
  "SEPARATOR TAB",
  "DATASET_LABEL\tNative continents (6 groups)",
  "COLOR\t#555555",
  paste("FIELD_SHAPES", paste(rep(1, length(continents)), collapse = "\t"), sep = "\t"),
  paste("FIELD_LABELS", paste(continents, collapse = "\t"), sep = "\t"),
  paste("FIELD_COLORS", paste(continent_colors, collapse = "\t"), sep = "\t"),
  "MARGIN\t5",
  "SYMBOL_SPACING\t2",
  "HEIGHT_FACTOR\t0.85",
  "SHOW_LABELS\t0",
  "LABEL_ROTATION\t45",
  "LEGEND_TITLE\tNative continents",
  "LEGEND_SCALE\t0.8",
  "LEGEND_HORIZONTAL\t0",
  paste("LEGEND_SHAPES", paste(rep(1, length(continents)), collapse = "\t"), sep = "\t"),
  paste("LEGEND_COLORS", paste(continent_colors, collapse = "\t"), sep = "\t"),
  paste("LEGEND_LABELS", paste(continents, collapse = "\t"), sep = "\t"),
  "DATA",
  apply(cbind(data$genus, continent_matrix), 1, paste, collapse = "\t")
)
write_itol(continent_lines, "04_itol_native_continents_6_corrected.txt")

tdwg_colors <- c("#D73027", "#FDAE61", "#FC8D59", "#FFD92F", "#4575B4", "#984EA3", "#E6AB02", "#1B9E77")
tdwg_columns <- paste0("tdwg_", tdwg_levels)
tdwg_matrix <- as.matrix(data[, tdwg_columns])
tdwg_matrix[tdwg_matrix == 0] <- -1

tdwg_lines <- c(
  "DATASET_BINARY",
  "SEPARATOR TAB",
  "DATASET_LABEL\tNative TDWG/POWO regions (8)",
  "COLOR\t#555555",
  paste("FIELD_SHAPES", paste(rep(1, length(tdwg_levels)), collapse = "\t"), sep = "\t"),
  paste("FIELD_LABELS", paste(tdwg_labels, collapse = "\t"), sep = "\t"),
  paste("FIELD_COLORS", paste(tdwg_colors, collapse = "\t"), sep = "\t"),
  "MARGIN\t5",
  "SYMBOL_SPACING\t2",
  "HEIGHT_FACTOR\t0.85",
  "SHOW_LABELS\t1",
  "LABEL_ROTATION\t45",
  "DATA",
  apply(cbind(data$genus, tdwg_matrix), 1, paste, collapse = "\t")
)
write_itol(tdwg_lines, "S2_itol_native_TDWG_regions_8.txt")

# -----------------------------------------------------------------------------
# 6. Final validation
# -----------------------------------------------------------------------------

validation <- c(
  "iTOL PACKAGE VALIDATION - ARECACEAE",
  paste("Tree tips:", ape::Ntip(tree)),
  paste("POWO genera:", length(powo_genera)),
  paste("Non-hybrid POWO species:", nrow(powo_species)),
  paste("Accepted studied species:", nrow(reconciled_studies)),
  paste("Genera with studies:", sum(data$has_studies)),
  paste("Genera without studies:", sum(data$has_studies == 0)),
  paste("Median accepted species per genus:", median(data$powo_species)),
  paste("Monotypic genera:", sum(data$powo_species == 1)),
  paste("Genera with 100% study coverage:", sum(data$study_coverage_percent == 100)),
  "Tree-data match: 184/184",
  "Geographic recoding: POWO/TDWG Level 2 regions",
  "Central America and Caribbean assigned to North America",
  "Northern, Western, and Southern South America plus Brazil assigned to South America"
)
write_itol(validation, "VALIDATION.txt")

message(paste(validation, collapse = "\n"))
message("\nFiles completed in: ", normalizePath(output_directory))
