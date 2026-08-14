get_script_path <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  script_argument <- grep("^--file=", arguments, value = TRUE)
  if (length(script_argument) == 1L) {
    return(normalizePath(sub("^--file=", "", script_argument), winslash = "/"))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_path)) {
      return(normalizePath(active_path, winslash = "/"))
    }
  }
  normalizePath(
    file.path(getwd(), "Arecaceae_single_manuscript_pipeline.R"),
    winslash = "/",
    mustWork = FALSE
  )
}

SCRIPT_PATH <- get_script_path()
SCRIPT_DIRECTORY <- dirname(SCRIPT_PATH)
PROJECT_DIRECTORY <- Sys.getenv(
  "ARECACEAE_PROJECT_DIR",
  unset = normalizePath(file.path(SCRIPT_DIRECTORY, ".."), winslash = "/", mustWork = FALSE)
)

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
set.seed(20260813)

required_packages <- c(
  "ape", "brglm2", "data.table", "ggplot2", "MASS", "patchwork",
  "rWCVP", "scales", "sf", "svglite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the missing R packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(ggplot2)
  library(MASS)
  library(patchwork)
  library(scales)
  library(sf)
})

if (!capabilities("cairo")) {
  stop("The installed R build does not provide Cairo support.", call. = FALSE)
}

FIGURE_FONT <- Sys.getenv(
  "ARECACEAE_FIGURE_FONT",
  unset = if (.Platform$OS.type == "windows") "Arial" else "sans"
)
PD_PERMUTATIONS <- suppressWarnings(as.integer(Sys.getenv(
  "ARECACEAE_PD_PERMUTATIONS",
  unset = "9999"
)))
if (!is.finite(PD_PERMUTATIONS) || PD_PERMUTATIONS < 999L) {
  stop("ARECACEAE_PD_PERMUTATIONS must be an integer of at least 999.", call. = FALSE)
}

INTERIM_DIRECTORY <- file.path(PROJECT_DIRECTORY, "data", "interim")
DERIVED_DIRECTORY <- file.path(PROJECT_DIRECTORY, "data", "derived")
REFERENCE_DIRECTORY <- file.path(PROJECT_DIRECTORY, "data", "reference")
RESULTS_DIRECTORY <- file.path(PROJECT_DIRECTORY, "results")
TABLE_DIRECTORY <- file.path(RESULTS_DIRECTORY, "tables")
MODEL_DIRECTORY <- file.path(RESULTS_DIRECTORY, "models")
MAIN_FIGURE_DIRECTORY <- file.path(RESULTS_DIRECTORY, "figures", "main")
SUPPLEMENTARY_FIGURE_DIRECTORY <- file.path(RESULTS_DIRECTORY, "figures", "supplementary")
ITOL_DIRECTORY <- file.path(RESULTS_DIRECTORY, "itol")
LOG_DIRECTORY <- file.path(RESULTS_DIRECTORY, "logs")
invisible(lapply(
  c(
    INTERIM_DIRECTORY, DERIVED_DIRECTORY, REFERENCE_DIRECTORY,
    RESULTS_DIRECTORY, TABLE_DIRECTORY, MODEL_DIRECTORY,
    MAIN_FIGURE_DIRECTORY, SUPPLEMENTARY_FIGURE_DIRECTORY,
    ITOL_DIRECTORY, LOG_DIRECTORY
  ),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

input_paths <- c(
  study_design = file.path(DERIVED_DIRECTORY, "study_design_summary.csv"),
  chemical_incidence = file.path(INTERIM_DIRECTORY, "04C_incidence_article_taxon_organ_entity.csv"),
  chemical_articles = file.path(INTERIM_DIRECTORY, "04C_articles_250.csv"),
  species_analysis = file.path(INTERIM_DIRECTORY, "02_species_analysis_with_GBIF_GHS.csv"),
  species_tree = file.path(REFERENCE_DIRECTORY, "01_species_tree_WCVP_pruned.tre"),
  genus_tree = file.path(REFERENCE_DIRECTORY, "Species_tree_all_genes_rooted_original.tre"),
  domain_architecture = file.path(DERIVED_DIRECTORY, "Table_04C_02_domain_architecture.csv"),
  class_architecture = file.path(DERIVED_DIRECTORY, "Table_04C_03_class_architecture.csv"),
  organ_architecture = file.path(DERIVED_DIRECTORY, "Table_04C_04_organ_architecture.csv"),
  organ_class = file.path(DERIVED_DIRECTORY, "Table_04C_05_organ_class_incidence.csv"),
  genus_class = file.path(DERIVED_DIRECTORY, "Table_04C_07_genus_class_incidence.csv"),
  genus_coverage = file.path(DERIVED_DIRECTORY, "Table_05B_04_genus_coverage.csv"),
  main_selection_models = file.path(DERIVED_DIRECTORY, "Table_05B_05_main_selection_models.csv"),
  trait_models = file.path(DERIVED_DIRECTORY, "Table_05B_07_trait_specific_models.csv"),
  socioeconomic_models = file.path(DERIVED_DIRECTORY, "Table_05B_08_socioeconomic_models.csv")
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop(
    "The following normalized repository inputs are missing:\n- ",
    paste(missing_inputs, collapse = "\n- "),
    call. = FALSE
  )
}

read_semicolon <- function(path) {
  data.table::fread(
    path,
    sep = ";",
    encoding = "UTF-8",
    na.strings = c("", "NA"),
    check.names = FALSE
  )
}

write_semicolon <- function(data, path) {
  data.table::fwrite(
    data.table::as.data.table(data),
    path,
    sep = ";",
    bom = TRUE,
    quote = TRUE,
    na = ""
  )
  invisible(path)
}

write_itol <- function(lines, filename) {
  path <- file.path(ITOL_DIRECTORY, filename)
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
  invisible(path)
}

trim_to_na <- function(x) {
  value <- trimws(as.character(x))
  value[is.na(value) | value == "" | value == "NA"] <- NA_character_
  value
}

to_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

to_binary <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }
  value <- tolower(trimws(as.character(x)))
  output <- rep(NA_integer_, length(value))
  output[value %in% c("1", "true", "yes")] <- 1L
  output[value %in% c("0", "false", "no")] <- 0L
  numeric_value <- suppressWarnings(as.numeric(value))
  output[is.na(output) & numeric_value %in% c(0, 1)] <- as.integer(
    numeric_value[is.na(output) & numeric_value %in% c(0, 1)]
  )
  output
}

first_nonmissing <- function(x) {
  value <- trim_to_na(x)
  value <- value[!is.na(value)]
  if (length(value) == 0L) NA_character_ else value[1L]
}

safe_unique_n <- function(x) {
  value <- trim_to_na(x)
  data.table::uniqueN(value[!is.na(value)])
}

z_score <- function(x) {
  value <- to_numeric(x)
  valid <- is.finite(value)
  output <- rep(NA_real_, length(value))
  if (sum(valid) < 2L) {
    return(output)
  }
  scale_value <- stats::sd(value[valid])
  if (!is.finite(scale_value) || scale_value == 0) {
    return(output)
  }
  output[valid] <- (value[valid] - mean(value[valid])) / scale_value
  output
}

label_wrap <- function(width) {
  function(x) vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1))
}

summary_value <- function(summary_table, indicator) {
  value <- summary_table[Indicator == indicator, Value]
  if (length(value) != 1L || !is.finite(as.numeric(value))) {
    stop("Invalid or missing study-design indicator: ", indicator, call. = FALSE)
  }
  as.numeric(value)
}

semantic_palette <- c(
  universe = "#17324D",
  documented = "#147D73",
  documented_light = "#B9DCD5",
  gap = "#D9822B",
  warning = "#B44E45",
  secondary = "#6F7F8D",
  missing = "#C9D1D6",
  zero = "#F5F6F4",
  grid = "#DDE3E6",
  text = "#233746",
  background = "#FFFFFF"
)

domain_palette <- c(
  "Phenylpropanoids and polyketides" = "#3F5F9A",
  "Terpenoids and steroids" = "#A64D8F",
  "Lipids and lipid-like molecules" = "#A68B21",
  "Other specialized metabolites" = "#2796A6",
  "Carbohydrates" = "#D96868",
  "Alkaloids and nitrogen compounds" = "#6F4C9B",
  "Other primary metabolites" = "#3A8F55",
  "Amino acids and peptides" = "#A66F2C",
  "Vitamins and cofactors" = "#5C97C7"
)

theme_manuscript <- function(base_size = 10.5) {
  ggplot2::theme_minimal(base_size = base_size, base_family = FIGURE_FONT) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = semantic_palette[["universe"]]),
      plot.subtitle = ggplot2::element_text(colour = semantic_palette[["secondary"]]),
      plot.caption = ggplot2::element_text(colour = semantic_palette[["secondary"]], hjust = 0),
      plot.tag = ggplot2::element_text(face = "bold", size = 13, colour = semantic_palette[["universe"]]),
      axis.title = ggplot2::element_text(face = "bold", colour = semantic_palette[["text"]]),
      axis.text = ggplot2::element_text(colour = semantic_palette[["text"]]),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = semantic_palette[["grid"]], linewidth = 0.35),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
}

blank_itol_panel <- function(title, subtitle) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "rect",
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1,
      fill = semantic_palette[["background"]],
      colour = semantic_palette[["grid"]],
      linewidth = 0.5
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void(base_family = FIGURE_FONT) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        colour = semantic_palette[["universe"]],
        size = 11
      ),
      plot.subtitle = ggplot2::element_text(
        colour = semantic_palette[["secondary"]],
        size = 9
      ),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
}

save_vector_figure <- function(plot, directory, filename, width, height) {
  pdf_path <- file.path(directory, paste0(filename, ".pdf"))
  svg_path <- file.path(directory, paste0(filename, ".svg"))
  png_path <- file.path(directory, paste0(filename, ".png"))
  grDevices::cairo_pdf(
    pdf_path,
    width = width,
    height = height,
    family = FIGURE_FONT,
    bg = semantic_palette[["background"]],
    onefile = FALSE
  )
  print(plot)
  grDevices::dev.off()
  svglite::svglite(
    svg_path,
    width = width,
    height = height,
    bg = semantic_palette[["background"]],
    system_fonts = list(sans = FIGURE_FONT)
  )
  print(plot)
  grDevices::dev.off()
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      png_path,
      plot = plot,
      device = ragg::agg_png,
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      bg = semantic_palette[["background"]]
    )
  } else {
    ggplot2::ggsave(
      png_path,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      bg = semantic_palette[["background"]]
    )
  }
  data.table::data.table(
    Figure = filename,
    PDF = pdf_path,
    SVG = svg_path,
    PNG = png_path,
    PDF_bytes = file.info(pdf_path)$size,
    SVG_bytes = file.info(svg_path)$size,
    PNG_bytes = file.info(png_path)$size
  )
}

study_design <- read_semicolon(input_paths[["study_design"]])
incidence <- read_semicolon(input_paths[["chemical_incidence"]])
articles <- read_semicolon(input_paths[["chemical_articles"]])
species <- read_semicolon(input_paths[["species_analysis"]])
domain_architecture <- read_semicolon(input_paths[["domain_architecture"]])
class_architecture <- read_semicolon(input_paths[["class_architecture"]])
organ_architecture <- read_semicolon(input_paths[["organ_architecture"]])
organ_class <- read_semicolon(input_paths[["organ_class"]])
genus_class <- read_semicolon(input_paths[["genus_class"]])
genus_coverage <- read_semicolon(input_paths[["genus_coverage"]])
main_selection_models <- read_semicolon(input_paths[["main_selection_models"]])
trait_models <- read_semicolon(input_paths[["trait_models"]])
socioeconomic_models <- read_semicolon(input_paths[["socioeconomic_models"]])
species_tree <- ape::read.tree(input_paths[["species_tree"]])
genus_tree <- ape::read.tree(input_paths[["genus_tree"]])

required_incidence_columns <- c(
  "Incidence_record_ID", "Article_ID", "WCVP_taxon_id",
  "Species_current_analysis", "Genus_current_analysis", "Organ_primary",
  "Chemical_entity_ID", "Domain_analysis", "Class_analysis", "Class_resolved",
  "Metabolite_scope", "Strict_sensitivity_scope", "Exact_structure_resolved",
  "Method_GC", "Method_LC", "Method_NMR_isolation", "Method_MS",
  "Method_other_assay", "Method_unspecified"
)
required_species_columns <- c(
  "WCVP_taxon_id", "Analysis_species_key", "Species_current", "Genus_current",
  "Research_studied_binary", "Tree_included", "Trait_data_available",
  "POWO_native_L3_count", "POWO_native_codes_L3", "POWO_native_continents",
  "POWO_native_regions_L2", "GBIF_GHS_data_available"
)
required_article_columns <- c("Article_ID", "Year")
required_derived_columns <- list(
  domain_architecture = c(
    "Domain_analysis", "Chemical_entity_count", "Exact_structure_count"
  ),
  class_architecture = c(
    "Domain_analysis", "Class_analysis", "Metabolite_scope", "Article_count",
    "Genus_count", "Chemical_entity_count"
  ),
  organ_architecture = c(
    "Organ_primary", "Article_count", "Chemical_entity_count"
  ),
  organ_class = c(
    "Organ_primary", "Class_analysis", "Species_prevalence_within_organ"
  ),
  genus_class = c(
    "Genus_current", "Class_analysis", "Species_prevalence_within_genus"
  ),
  genus_coverage = c(
    "Genus_current", "Studied_species", "Total_species", "Representation_rate"
  ),
  main_selection_models = c(
    "Term", "P_value", "Model", "Odds_ratio", "CI_low", "CI_high"
  ),
  trait_models = c(
    "Trait", "P_value_Holm", "Holm_significant", "Odds_ratio", "CI_low", "CI_high"
  ),
  socioeconomic_models = c(
    "Term", "P_value", "Odds_ratio", "CI_low", "CI_high"
  )
)
if (length(setdiff(required_incidence_columns, names(incidence))) > 0L) {
  stop("The normalized chemical incidence table is missing required columns.", call. = FALSE)
}
if (length(setdiff(required_species_columns, names(species))) > 0L) {
  stop("The normalized species analysis table is missing required columns.", call. = FALSE)
}
if (length(setdiff(required_article_columns, names(articles))) > 0L) {
  stop("The normalized article table is missing Article_ID or Year.", call. = FALSE)
}
derived_objects <- list(
  domain_architecture = domain_architecture,
  class_architecture = class_architecture,
  organ_architecture = organ_architecture,
  organ_class = organ_class,
  genus_class = genus_class,
  genus_coverage = genus_coverage,
  main_selection_models = main_selection_models,
  trait_models = trait_models,
  socioeconomic_models = socioeconomic_models
)
for (object_name in names(required_derived_columns)) {
  missing_columns <- setdiff(
    required_derived_columns[[object_name]],
    names(derived_objects[[object_name]])
  )
  if (length(missing_columns) > 0L) {
    stop(
      "The normalized table ", object_name,
      " is missing: ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

incidence[, `:=`(
  Article_ID = trim_to_na(Article_ID),
  WCVP_taxon_id = sub("\\.0$", "", trim_to_na(WCVP_taxon_id)),
  Species_current_analysis = trim_to_na(Species_current_analysis),
  Genus_current_analysis = trim_to_na(Genus_current_analysis),
  Organ_primary = trim_to_na(Organ_primary),
  Chemical_entity_ID = trim_to_na(Chemical_entity_ID),
  Domain_analysis = trim_to_na(Domain_analysis),
  Class_analysis = trim_to_na(Class_analysis)
)]
binary_incidence_columns <- c(
  "Class_resolved", "Metabolite_scope", "Strict_sensitivity_scope",
  "Exact_structure_resolved", "Method_GC", "Method_LC",
  "Method_NMR_isolation", "Method_MS", "Method_other_assay", "Method_unspecified"
)
incidence[, (binary_incidence_columns) := lapply(.SD, to_binary), .SDcols = binary_incidence_columns]
species[, `:=`(
  WCVP_taxon_id = sub("\\.0$", "", trim_to_na(WCVP_taxon_id)),
  Analysis_species_key = trim_to_na(Analysis_species_key),
  Species_current = trim_to_na(Species_current),
  Genus_current = trim_to_na(Genus_current),
  POWO_native_codes_L3 = trim_to_na(POWO_native_codes_L3),
  POWO_native_continents = trim_to_na(POWO_native_continents),
  POWO_native_regions_L2 = trim_to_na(POWO_native_regions_L2),
  Research_studied_binary = to_binary(Research_studied_binary),
  Tree_included = to_binary(Tree_included),
  Trait_data_available = to_binary(Trait_data_available),
  GBIF_GHS_data_available = to_binary(GBIF_GHS_data_available)
)]
articles[, `:=`(
  Article_ID = trim_to_na(Article_ID),
  Year = to_numeric(Year)
)]

expected_summary <- c(
  Screened_documents = 389,
  Included_articles = 250,
  Excluded_documents = 139,
  Original_evidence_records = 6287,
  Deduplicated_analytical_incidences = 5207,
  Curated_chemical_entities = 1677,
  Biological_scope_entities = 1645,
  Current_species = 2582,
  Studied_current_species = 173,
  Family_genera = 189,
  Studied_current_genera = 75,
  Phylogeny_species = 2233,
  Studied_phylogeny_species = 172
)
summary_observed <- setNames(to_numeric(study_design$Value), study_design$Indicator)
missing_summary <- setdiff(names(expected_summary), names(summary_observed))
if (length(missing_summary) > 0L) {
  stop("The study-design summary is missing required indicators.", call. = FALSE)
}
if (any(summary_observed[names(expected_summary)] != expected_summary)) {
  stop("The frozen study-design counts do not match the audited analysis.", call. = FALSE)
}
if (nrow(incidence) != 5207L || uniqueN(articles$Article_ID) != 250L) {
  stop("The normalized chemical corpus does not contain 5,207 incidences and 250 articles.", call. = FALSE)
}
if (nrow(species) != 2582L || uniqueN(species$WCVP_taxon_id) != 2582L) {
  stop("The species analysis table does not contain 2,582 unique WCVP species.", call. = FALSE)
}
if (!inherits(species_tree, "phylo") || ape::Ntip(species_tree) != 2233L) {
  stop("The species tree does not contain the expected 2,233 WCVP tips.", call. = FALSE)
}
if (!ape::is.rooted(species_tree) || !ape::is.binary(species_tree)) {
  stop("The species tree must be rooted and binary.", call. = FALSE)
}
if (is.null(species_tree$edge.length) || any(!is.finite(species_tree$edge.length)) || any(species_tree$edge.length <= 0)) {
  stop("The species tree must contain finite positive branch lengths.", call. = FALSE)
}

message("[01/05] Computing species-level effort and richness metrics...")

main_incidence <- incidence[
  Metabolite_scope == 1L &
    !is.na(WCVP_taxon_id) &
    !is.na(Species_current_analysis)
]
if (uniqueN(main_incidence$WCVP_taxon_id) != 173L) {
  stop("The curated chemical layer must contain 173 studied WCVP species.", call. = FALSE)
}

species_effort <- main_incidence[, .(
  Species_current = first_nonmissing(Species_current_analysis),
  Genus_current = first_nonmissing(Genus_current_analysis),
  Article_count = uniqueN(Article_ID),
  Entity_richness_observed = uniqueN(Chemical_entity_ID),
  Exact_structure_richness = uniqueN(Chemical_entity_ID[Exact_structure_resolved == 1L]),
  Strict_entity_richness = uniqueN(Chemical_entity_ID[Strict_sensitivity_scope == 1L]),
  Class_richness = safe_unique_n(Class_analysis),
  Domain_richness = safe_unique_n(Domain_analysis),
  Organ_count = safe_unique_n(Organ_primary),
  Method_GC = as.integer(any(Method_GC == 1L, na.rm = TRUE)),
  Method_LC = as.integer(any(Method_LC == 1L, na.rm = TRUE)),
  Method_NMR_isolation = as.integer(any(Method_NMR_isolation == 1L, na.rm = TRUE)),
  Method_MS = as.integer(any(Method_MS == 1L, na.rm = TRUE))
), by = WCVP_taxon_id]
species_effort[, `:=`(
  Core_method_count = Method_GC + Method_LC + Method_NMR_isolation + Method_MS,
  Article_z = z_score(log1p(Article_count)),
  Organ_z = z_score(Organ_count),
  Method_z = z_score(Method_GC + Method_LC + Method_NMR_isolation + Method_MS)
)]

effort_models <- list(
  entities = MASS::glm.nb(
    Entity_richness_observed ~ Article_z + Organ_z + Method_z,
    data = species_effort,
    control = stats::glm.control(maxit = 200)
  ),
  exact = MASS::glm.nb(
    Exact_structure_richness ~ Article_z + Organ_z + Method_z,
    data = species_effort,
    control = stats::glm.control(maxit = 200)
  ),
  strict = MASS::glm.nb(
    Strict_entity_richness ~ Article_z + Organ_z + Method_z,
    data = species_effort,
    control = stats::glm.control(maxit = 200)
  ),
  classes = MASS::glm.nb(
    Class_richness ~ Article_z + Organ_z + Method_z,
    data = species_effort,
    control = stats::glm.control(maxit = 200)
  )
)

append_effort_score <- function(data, model, response, prefix) {
  prediction <- stats::predict(model, type = "link", se.fit = TRUE)
  observed <- data[[response]]
  expected_column <- paste0(prefix, "_expected")
  expected_low_column <- paste0(prefix, "_expected_low")
  expected_high_column <- paste0(prefix, "_expected_high")
  score_column <- paste0(prefix, "_conditional_score")
  data[, (expected_column) := exp(prediction$fit)]
  data[, (expected_low_column) := exp(prediction$fit - 1.96 * prediction$se.fit)]
  data[, (expected_high_column) := exp(prediction$fit + 1.96 * prediction$se.fit)]
  data[, (score_column) := log2(
    (observed + 0.5) / (exp(prediction$fit) + 0.5)
  )]
  data
}

species_effort <- append_effort_score(species_effort, effort_models$entities, "Entity_richness_observed", "Entity")
species_effort <- append_effort_score(species_effort, effort_models$exact, "Exact_structure_richness", "Exact")
species_effort <- append_effort_score(species_effort, effort_models$strict, "Strict_entity_richness", "Strict")
species_effort <- append_effort_score(species_effort, effort_models$classes, "Class_richness", "Class")

article_entity <- unique(main_incidence[, .(WCVP_taxon_id, Article_ID, Chemical_entity_ID)])
two_article_distribution <- article_entity[, {
  article_ids <- sort(unique(Article_ID))
  if (length(article_ids) < 2L) {
    list(
      Two_article_expected = NA_real_,
      Two_article_low = NA_real_,
      Two_article_high = NA_real_,
      Two_article_pair_count = 0L
    )
  } else {
    pairs <- utils::combn(article_ids, 2L)
    richness <- apply(pairs, 2L, function(pair) {
      uniqueN(Chemical_entity_ID[Article_ID %in% pair])
    })
    list(
      Two_article_expected = mean(richness),
      Two_article_low = as.numeric(stats::quantile(richness, 0.025, names = FALSE, type = 8)),
      Two_article_high = as.numeric(stats::quantile(richness, 0.975, names = FALSE, type = 8)),
      Two_article_pair_count = length(richness)
    )
  }
}, by = WCVP_taxon_id]
species_effort <- merge(
  species_effort,
  two_article_distribution,
  by = "WCVP_taxon_id",
  all.x = TRUE,
  sort = FALSE
)

sensitivity_correlations <- data.table(
  Comparison = c(
    "Entity versus exact-structure conditional score",
    "Entity versus strict-entity conditional score",
    "Entity versus chemical-class conditional score"
  ),
  N = nrow(species_effort),
  Spearman_rho = c(
    stats::cor(species_effort$Entity_conditional_score, species_effort$Exact_conditional_score, method = "spearman"),
    stats::cor(species_effort$Entity_conditional_score, species_effort$Strict_conditional_score, method = "spearman"),
    stats::cor(species_effort$Entity_conditional_score, species_effort$Class_conditional_score, method = "spearman")
  )
)

write_semicolon(species_effort, file.path(TABLE_DIRECTORY, "Table_species_effort_conditioned_richness.csv"))
write_semicolon(sensitivity_correlations, file.path(TABLE_DIRECTORY, "Table_effort_score_sensitivity.csv"))
saveRDS(effort_models, file.path(MODEL_DIRECTORY, "effort_conditioned_richness_models.rds"), compress = "xz")

message("[02/05] Computing phylogenetic coverage...")

studied_tree_tips <- intersect(
  species_tree$tip.label,
  species[Research_studied_binary == 1L, Analysis_species_key]
)
if (length(studied_tree_tips) != 172L) {
  stop("Exactly 172 studied species must match the species tree.", call. = FALSE)
}

species_tree_postorder <- ape::reorder.phylo(species_tree, order = "postorder")
faith_pd_from_indices <- function(tip_indices) {
  included_nodes <- logical(ape::Ntip(species_tree_postorder) + ape::Nnode(species_tree_postorder))
  included_nodes[tip_indices] <- TRUE
  total <- 0
  for (edge_index in seq_len(nrow(species_tree_postorder$edge))) {
    parent <- species_tree_postorder$edge[edge_index, 1L]
    child <- species_tree_postorder$edge[edge_index, 2L]
    if (included_nodes[child]) {
      total <- total + species_tree_postorder$edge.length[edge_index]
      included_nodes[parent] <- TRUE
    }
  }
  total
}

studied_tree_indices <- match(studied_tree_tips, species_tree_postorder$tip.label)
observed_pd <- faith_pd_from_indices(studied_tree_indices)
full_pd <- sum(species_tree$edge.length)
set.seed(20260813)
null_pd <- replicate(
  PD_PERMUTATIONS,
  faith_pd_from_indices(sample.int(
    ape::Ntip(species_tree_postorder),
    length(studied_tree_indices),
    replace = FALSE
  ))
)
null_mean <- mean(null_pd)
null_sd <- stats::sd(null_pd)
pd_ses <- (observed_pd - null_mean) / null_sd
pd_two_sided_p <- (
  sum(abs(null_pd - null_mean) >= abs(observed_pd - null_mean)) + 1
) / (PD_PERMUTATIONS + 1)

phylogenetic_diversity <- data.table(
  Metric = c(
    "Species_tree_tips", "Studied_species_in_tree", "Observed_Faith_PD",
    "Full_tree_Faith_PD", "Observed_PD_fraction", "Null_PD_mean",
    "Null_PD_SD", "SES_PD", "Two_sided_P", "Randomizations"
  ),
  Value = c(
    ape::Ntip(species_tree), length(studied_tree_tips), observed_pd,
    full_pd, observed_pd / full_pd, null_mean,
    null_sd, pd_ses, pd_two_sided_p, PD_PERMUTATIONS
  )
)
write_semicolon(
  phylogenetic_diversity,
  file.path(TABLE_DIRECTORY, "Table_phylogenetic_diversity_coverage.csv")
)

message("[03/05] Fitting article-level method-by-class models...")

article_method <- incidence[, .(
  Method_GC = as.integer(any(Method_GC == 1L, na.rm = TRUE)),
  Method_LC = as.integer(any(Method_LC == 1L, na.rm = TRUE)),
  Method_NMR_isolation = as.integer(any(Method_NMR_isolation == 1L, na.rm = TRUE)),
  Method_MS = as.integer(any(Method_MS == 1L, na.rm = TRUE)),
  Method_other_assay = as.integer(any(Method_other_assay == 1L, na.rm = TRUE)),
  Method_unspecified = as.integer(any(Method_unspecified == 1L, na.rm = TRUE)),
  Organ_count = safe_unique_n(Organ_primary),
  Organ_profile = {
    organs <- sort(unique(Organ_primary[!is.na(Organ_primary)]))
    if (length(organs) == 0L) {
      NA_character_
    } else if (length(organs) == 1L) {
      organs
    } else {
      "Multiple organs"
    }
  }
), by = Article_ID]
article_method <- merge(
  article_method,
  unique(articles[, .(Article_ID, Year)]),
  by = "Article_ID",
  all.x = TRUE,
  sort = FALSE
)
article_method[, Year_z := z_score(Year)]

article_class_presence <- unique(incidence[
  Metabolite_scope == 1L & Class_resolved == 1L & !is.na(Class_analysis),
  .(Article_ID, Class_analysis)
])
leading_classes <- article_class_presence[, .(
  Class_article_count = uniqueN(Article_ID)
), by = Class_analysis][Class_article_count >= 20L][order(-Class_article_count)]
if (nrow(leading_classes) != 12L) {
  stop("The method-by-class analysis expects 12 classes documented in at least 20 articles.", call. = FALSE)
}

method_grid <- data.table::CJ(
  Article_ID = sort(unique(articles$Article_ID)),
  Class_analysis = leading_classes$Class_analysis,
  unique = TRUE
)
method_grid <- merge(
  method_grid,
  article_class_presence[, .(Article_ID, Class_analysis, Class_present = 1L)],
  by = c("Article_ID", "Class_analysis"),
  all.x = TRUE,
  sort = FALSE
)
method_grid[is.na(Class_present), Class_present := 0L]
method_grid <- merge(method_grid, article_method, by = "Article_ID", all.x = TRUE, sort = FALSE)

method_terms <- c("Method_GC", "Method_LC", "Method_NMR_isolation", "Method_MS")
method_labels <- c(
  Method_GC = "GC-based",
  Method_LC = "LC-based",
  Method_NMR_isolation = "Isolation/NMR",
  Method_MS = "Mass spectrometry"
)

fit_method_class <- function(class_name) {
  analysis_data <- method_grid[Class_analysis == class_name]
  fit <- stats::glm(
    Class_present ~ Method_GC + Method_LC + Method_NMR_isolation + Method_MS,
    data = analysis_data,
    family = stats::binomial(),
    method = brglm2::brglmFit,
    type = "AS_mean",
    control = stats::glm.control(maxit = 300)
  )
  coefficient_table <- summary(fit)$coefficients
  output <- data.table(
    Class_analysis = class_name,
    Term = rownames(coefficient_table),
    Estimate = coefficient_table[, "Estimate"],
    Standard_error = coefficient_table[, "Std. Error"],
    Statistic = coefficient_table[, "z value"],
    P_value = coefficient_table[, "Pr(>|z|)"],
    N_articles = nrow(analysis_data),
    N_class_articles = sum(analysis_data$Class_present),
    Converged = isTRUE(fit$converged)
  )
  output[Term %in% method_terms]
}

method_class_models <- rbindlist(
  lapply(leading_classes$Class_analysis, fit_method_class),
  use.names = TRUE,
  fill = TRUE
)
method_class_models[, `:=`(
  Analytical_family = unname(method_labels[Term]),
  Odds_ratio = exp(Estimate),
  CI_low = exp(Estimate - 1.96 * Standard_error),
  CI_high = exp(Estimate + 1.96 * Standard_error)
)]
method_class_models[, P_value_Holm := stats::p.adjust(P_value, method = "holm")]
method_class_models[, Supported := as.integer(P_value_Holm < 0.05)]
method_class_models <- merge(
  method_class_models,
  leading_classes,
  by = "Class_analysis",
  all.x = TRUE,
  sort = FALSE
)
if (nrow(method_class_models) != 48L || any(method_class_models$Converged != 1L)) {
  write_semicolon(
    method_class_models,
    file.path(TABLE_DIRECTORY, "Table_method_class_models.csv")
  )
  stop("The method-by-class model audit failed.", call. = FALSE)
}
write_semicolon(method_class_models, file.path(TABLE_DIRECTORY, "Table_method_class_models.csv"))
saveRDS(method_class_models, file.path(MODEL_DIRECTORY, "method_class_models.rds"), compress = "xz")

message("[04/05] Building narrative figures...")

sample_replacements <- c(
  Genus_species4 = "Truongsonia_lecongkietii",
  Genus_species7 = "Areca_chaiana",
  Genus_species8 = "Oncosperma_fasciculatum"
)
replacement_index <- match(names(sample_replacements), genus_tree$tip.label)
genus_tree$tip.label[replacement_index[!is.na(replacement_index)]] <-
  sample_replacements[!is.na(replacement_index)]
genus_tree <- ape::drop.tip(
  genus_tree,
  intersect(
    c("Dasypogon_bromeliifolius", "Areca_chaiana", "Oncosperma_fasciculatum"),
    genus_tree$tip.label
  ),
  collapse.singles = TRUE
)
genus_tree$tip.label <- sub("_.*$", "", genus_tree$tip.label)
genus_tree$tip.label[genus_tree$tip.label == "Acoelorrhaphe"] <- "Acoelorraphe"
genus_tree$node.label <- NULL
if (ape::Ntip(genus_tree) != 184L || anyDuplicated(genus_tree$tip.label)) {
  stop("The genus tree must contain 184 unique Arecaceae genera.", call. = FALSE)
}

genus_coverage[, `:=`(
  Genus_current = trim_to_na(Genus_current),
  Studied_species = to_numeric(Studied_species),
  Total_species = to_numeric(Total_species),
  Representation_rate = to_numeric(Representation_rate)
)]
genus_tree_coverage <- merge(
  data.table(Genus_current = genus_tree$tip.label),
  genus_coverage,
  by = "Genus_current",
  all.x = TRUE,
  sort = FALSE
)
genus_tree_coverage[is.na(Studied_species), Studied_species := 0]
genus_tree_coverage[is.na(Representation_rate), Representation_rate := 0]
if (any(!is.finite(genus_tree_coverage$Total_species))) {
  stop("Genus richness is missing for one or more genera in the phylogeny.", call. = FALSE)
}
genus_tree_coverage[, Coverage_status := ifelse(
  Studied_species > 0,
  "Chemically documented",
  "No chemical evidence"
)]

tdwg_levels <- c(
  "AFRICA", "ASIA-TEMPERATE", "ASIA-TROPICAL", "AUSTRALASIA",
  "EUROPE", "NORTHERN AMERICA", "PACIFIC", "SOUTHERN AMERICA"
)
tdwg_labels <- c(
  "Africa", "Temperate Asia", "Tropical Asia", "Australasia",
  "Europe", "Northern America", "Pacific", "Southern America"
)
continent_labels <- c(
  "Africa", "Asia", "Europe", "North America", "South America", "Oceania"
)

genus_tdwg_long <- species[
  Genus_current %in% genus_tree$tip.label & !is.na(POWO_native_continents),
  .(
    TDWG_region = toupper(trimws(unlist(strsplit(
      POWO_native_continents,
      "\\s*\\|\\s*"
    ))))
  ),
  by = Genus_current
]
genus_tdwg_long <- unique(genus_tdwg_long[TDWG_region %in% tdwg_levels])

genus_region_long <- species[
  Genus_current %in% genus_tree$tip.label & !is.na(POWO_native_regions_L2),
  .(
    Region_L2 = trimws(unlist(strsplit(
      POWO_native_regions_L2,
      "\\s*\\|\\s*"
    )))
  ),
  by = Genus_current
]
genus_region_long <- unique(genus_region_long[nzchar(Region_L2)])

north_american_southern_regions <- c("Caribbean", "Central America")
south_american_regions <- c(
  "Brazil", "Northern South America", "Western South America",
  "Southern South America"
)

itol_genus_data <- data.table::copy(genus_tree_coverage)
itol_genus_data <- itol_genus_data[
  match(genus_tree$tip.label, Genus_current)
]
itol_genus_data[, `:=`(
  Study_coverage_percent = 100 * Studied_species / Total_species,
  Has_studies = as.integer(Studied_species > 0),
  Continent_Africa = as.integer(
    Genus_current %in% genus_tdwg_long[TDWG_region == "AFRICA", Genus_current]
  ),
  Continent_Asia = as.integer(
    Genus_current %in% genus_tdwg_long[
      TDWG_region %in% c("ASIA-TEMPERATE", "ASIA-TROPICAL"),
      Genus_current
    ]
  ),
  Continent_Europe = as.integer(
    Genus_current %in% genus_tdwg_long[TDWG_region == "EUROPE", Genus_current]
  ),
  Continent_North_America = as.integer(
    Genus_current %in% genus_tdwg_long[
      TDWG_region == "NORTHERN AMERICA",
      Genus_current
    ] |
      Genus_current %in% genus_region_long[
        Region_L2 %in% north_american_southern_regions,
        Genus_current
      ]
  ),
  Continent_South_America = as.integer(
    Genus_current %in% genus_region_long[
      Region_L2 %in% south_american_regions,
      Genus_current
    ]
  ),
  Continent_Oceania = as.integer(
    Genus_current %in% genus_tdwg_long[
      TDWG_region %in% c("AUSTRALASIA", "PACIFIC"),
      Genus_current
    ]
  )
)]

for (tdwg_index in seq_along(tdwg_levels)) {
  tdwg_column <- paste0("TDWG_", gsub("[^A-Z]+", "_", tdwg_levels[tdwg_index]))
  itol_genus_data[, (tdwg_column) := as.integer(
    Genus_current %in% genus_tdwg_long[
      TDWG_region == tdwg_levels[tdwg_index],
      Genus_current
    ]
  )]
}

astral_genus_tree <- genus_tree
cladogram_genus_tree <- genus_tree
cladogram_genus_tree$edge.length <- NULL
ape::write.tree(
  cladogram_genus_tree,
  file.path(ITOL_DIRECTORY, "01_arecaceae_184_genera_cladogram.nwk")
)
ape::write.tree(
  astral_genus_tree,
  file.path(ITOL_DIRECTORY, "S3_arecaceae_184_genera_astral_branch_lengths.nwk")
)

richness_lines <- c(
  "DATASET_SIMPLEBAR",
  "SEPARATOR TAB",
  "DATASET_LABEL\tGenus richness (current WCVP species; log scale)",
  paste0("COLOR\t", semantic_palette[["secondary"]]),
  paste(
    "DATASET_SCALE", 1, 10, 100,
    max(itol_genus_data$Total_species),
    sep = "\t"
  ),
  "WIDTH\t80",
  "MARGIN\t5",
  "LOG_SCALE\t1",
  "HEIGHT_FACTOR\t0.65",
  "SHOW_VALUE\t0",
  "SHOW_LABELS\t0",
  "BORDER_WIDTH\t0",
  "DATA",
  sprintf(
    "%s\t%d",
    itol_genus_data$Genus_current,
    as.integer(itol_genus_data$Total_species)
  )
)
write_itol(richness_lines, "02_itol_genus_richness_logbar.txt")

coverage_lines <- c(
  "DATASET_GRADIENT",
  "SEPARATOR TAB",
  "DATASET_LABEL\tStudy coverage (% of current species)",
  paste0("COLOR\t", semantic_palette[["documented"]]),
  "STRIP_WIDTH\t35",
  "MARGIN\t5",
  "BORDER_WIDTH\t0",
  "AUTO_LEGEND\t0",
  paste0("COLOR_MIN\t", semantic_palette[["zero"]]),
  "USE_MID_COLOR\t1",
  paste0("COLOR_MID\t", semantic_palette[["documented_light"]]),
  paste0("COLOR_MAX\t", semantic_palette[["documented"]]),
  "SHOW_LABELS\t0",
  "LEGEND_TITLE\tStudy coverage (% of current species)",
  "LEGEND_SCALE\t0.8",
  "LEGEND_HORIZONTAL\t0",
  "LEGEND_SHAPES\t1\t1\t1\t1\t1",
  paste(
    "LEGEND_COLORS",
    semantic_palette[["zero"]],
    semantic_palette[["documented_light"]],
    "#72B5A8",
    semantic_palette[["documented"]],
    semantic_palette[["universe"]],
    sep = "\t"
  ),
  "LEGEND_LABELS\t0\t25\t50\t75\t100",
  "DATA",
  sprintf(
    "%s\t%.6f",
    itol_genus_data$Genus_current,
    itol_genus_data$Study_coverage_percent
  )
)
write_itol(coverage_lines, "03_itol_study_coverage_percent.txt")

study_presence_lines <- c(
  "DATASET_BINARY",
  "SEPARATOR TAB",
  "DATASET_LABEL\tStudy presence",
  paste0("COLOR\t", semantic_palette[["documented"]]),
  "FIELD_SHAPES\t2",
  "FIELD_LABELS\tAt least one chemically studied species",
  paste0("FIELD_COLORS\t", semantic_palette[["documented"]]),
  "MARGIN\t5",
  "HEIGHT_FACTOR\t0.8",
  "SHOW_LABELS\t1",
  "DATA",
  sprintf(
    "%s\t%d",
    itol_genus_data$Genus_current,
    itol_genus_data$Has_studies
  )
)
write_itol(study_presence_lines, "S1_itol_study_presence.txt")

continent_columns <- c(
  "Continent_Africa", "Continent_Asia", "Continent_Europe",
  "Continent_North_America", "Continent_South_America", "Continent_Oceania"
)
continent_colors <- c(
  "#D55E00", "#E69F00", "#0072B2", "#CC79A7", "#009E73", "#56B4E9"
)
continent_matrix <- as.matrix(itol_genus_data[, ..continent_columns])
continent_matrix[continent_matrix == 0] <- -1
continent_lines <- c(
  "DATASET_BINARY",
  "SEPARATOR TAB",
  "DATASET_LABEL\tNative continents (six groups)",
  "COLOR\t#555555",
  paste("FIELD_SHAPES", paste(rep(1, 6), collapse = "\t"), sep = "\t"),
  paste("FIELD_LABELS", paste(continent_labels, collapse = "\t"), sep = "\t"),
  paste("FIELD_COLORS", paste(continent_colors, collapse = "\t"), sep = "\t"),
  "MARGIN\t5",
  "SYMBOL_SPACING\t2",
  "HEIGHT_FACTOR\t0.85",
  "SHOW_LABELS\t0",
  "LABEL_ROTATION\t45",
  "LEGEND_TITLE\tNative continents",
  "LEGEND_SCALE\t0.8",
  "LEGEND_HORIZONTAL\t0",
  paste("LEGEND_SHAPES", paste(rep(1, 6), collapse = "\t"), sep = "\t"),
  paste("LEGEND_COLORS", paste(continent_colors, collapse = "\t"), sep = "\t"),
  paste("LEGEND_LABELS", paste(continent_labels, collapse = "\t"), sep = "\t"),
  "DATA",
  apply(
    cbind(itol_genus_data$Genus_current, continent_matrix),
    1,
    paste,
    collapse = "\t"
  )
)
write_itol(continent_lines, "04_itol_native_continents_6.txt")

tdwg_columns <- paste0("TDWG_", gsub("[^A-Z]+", "_", tdwg_levels))
tdwg_colors <- c(
  "#D73027", "#FDAE61", "#FC8D59", "#FFD92F",
  "#4575B4", "#984EA3", "#E6AB02", "#1B9E77"
)
tdwg_matrix <- as.matrix(itol_genus_data[, ..tdwg_columns])
tdwg_matrix[tdwg_matrix == 0] <- -1
tdwg_lines <- c(
  "DATASET_BINARY",
  "SEPARATOR TAB",
  "DATASET_LABEL\tNative TDWG/POWO regions (eight groups)",
  "COLOR\t#555555",
  paste("FIELD_SHAPES", paste(rep(1, 8), collapse = "\t"), sep = "\t"),
  paste("FIELD_LABELS", paste(tdwg_labels, collapse = "\t"), sep = "\t"),
  paste("FIELD_COLORS", paste(tdwg_colors, collapse = "\t"), sep = "\t"),
  "MARGIN\t5",
  "SYMBOL_SPACING\t2",
  "HEIGHT_FACTOR\t0.85",
  "SHOW_LABELS\t1",
  "LABEL_ROTATION\t45",
  "DATA",
  apply(
    cbind(itol_genus_data$Genus_current, tdwg_matrix),
    1,
    paste,
    collapse = "\t"
  )
)
write_itol(tdwg_lines, "S2_itol_native_TDWG_regions_8.txt")

itol_genus_export <- data.table::copy(itol_genus_data)
itol_genus_export[, Native_continents := apply(
  .SD,
  1,
  function(values) paste(continent_labels[as.logical(as.integer(values))], collapse = "; ")
), .SDcols = continent_columns]
write_semicolon(
  itol_genus_export,
  file.path(ITOL_DIRECTORY, "05_genus_level_data.csv")
)

load_wgsrpd3 <- function() {
  map_environment <- new.env(parent = emptyenv())
  suppressWarnings(utils::data("wgsrpd3", package = "rWCVP", envir = map_environment))
  if (!exists("wgsrpd3", envir = map_environment, inherits = FALSE)) {
    stop("The rWCVP package does not provide the wgsrpd3 spatial layer.", call. = FALSE)
  }
  map <- get("wgsrpd3", envir = map_environment)
  if (!inherits(map, "sf")) {
    stop("The wgsrpd3 object is not an sf spatial layer.", call. = FALSE)
  }
  code_candidates <- c(
    "LEVEL3_COD", "LEVEL3_CODE", "level3_cod", "level3_code",
    "area_code_l3", "LEVEL3"
  )
  code_column <- intersect(code_candidates, names(map))[1L]
  if (is.na(code_column)) {
    stop("No WGSRPD level-three code column was found.", call. = FALSE)
  }
  map$Area_code_L3 <- toupper(trimws(as.character(map[[code_column]])))
  map <- sf::st_make_valid(map)
  map <- sf::st_transform(map, 4326)
  map <- map[!sf::st_is_empty(map), c("Area_code_L3")]
  map
}

native_area_long <- species[
  !is.na(POWO_native_codes_L3),
  .(
    Area_code_L3 = toupper(trimws(unlist(strsplit(
      POWO_native_codes_L3,
      "\\s*\\|\\s*"
    ))))
  ),
  by = .(WCVP_taxon_id, Research_studied_binary)
]
native_area_long <- unique(native_area_long[nzchar(Area_code_L3)])
area_coverage <- native_area_long[, .(
  Total_species = uniqueN(WCVP_taxon_id),
  Studied_species = uniqueN(WCVP_taxon_id[Research_studied_binary == 1L])
), by = Area_code_L3]
area_coverage[, Representation_rate := Studied_species / Total_species]
world_map <- merge(
  load_wgsrpd3(),
  area_coverage,
  by = "Area_code_L3",
  all.x = TRUE,
  sort = FALSE
)
world_map$Palm_region <- !is.na(world_map$Total_species)
world_map$Coverage_group <- ifelse(
  !world_map$Palm_region,
  "Outside native palm range",
  ifelse(world_map$Studied_species > 0, "Documented palm region", "Undocumented palm region")
)
write_semicolon(
  sf::st_drop_geometry(world_map),
  file.path(TABLE_DIRECTORY, "Table_geographic_research_coverage.csv")
)

visual_identity <- rbindlist(list(
  data.table(
    Palette = "Semantic",
    Meaning = names(semantic_palette),
    Hex = unname(semantic_palette)
  ),
  data.table(
    Palette = "Chemical domain",
    Meaning = names(domain_palette),
    Hex = unname(domain_palette)
  )
))
write_semicolon(visual_identity, file.path(TABLE_DIRECTORY, "Table_visual_identity.csv"))

term_labels <- c(
  z_log_native_range_L3 = "Native range breadth",
  z_Palm_stature_PC1 = "Palm stature",
  z_Fruit_size_PC1 = "Fruit size",
  z_log_GBIF_10km_cells = "GBIF geographic coverage",
  z_GHS_mean_exposure = "Human exposure",
  Trait_leaf_blade_z = "Maximum leaf-blade length",
  Trait_climbing = "Climbing habit",
  Trait_acaulescent = "Acaulescent habit",
  Trait_erect = "Erect habit",
  Trait_solitary = "Solitary stem",
  Trait_stem_armed = "Stem armature",
  Trait_leaves_armed = "Leaf armature",
  Trait_canopy = "Canopy position",
  Commodity_species_resolved = "Species-resolved commodity",
  Human_use_any = "Documented human use",
  Economic_quantitative_data_available = "Quantitative economic data"
)

forest_panel <- function(data, title, subtitle = NULL, model_shapes = NULL) {
  plot_data <- data.table::copy(data)
  plot_data[, Display_term := factor(
    Display_term,
    levels = rev(unique(Display_term))
  )]
  plot_data[, Evidence := factor(
    ifelse(Supported == 1L, "Supported", "Not supported"),
    levels = c("Supported", "Not supported")
  )]
  has_multiple_models <- !is.null(model_shapes) && "Model_label" %in% names(plot_data)
  plot_data[, Plot_group := if (has_multiple_models) Model_label else "Single model"]
  point_position <- if (has_multiple_models) {
    ggplot2::position_dodge(width = 0.48)
  } else {
    ggplot2::position_identity()
  }
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Odds_ratio,
      y = Display_term,
      colour = Evidence,
      group = Plot_group
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = semantic_palette[["secondary"]],
      linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = CI_low, xmax = CI_high),
      width = 0,
      linewidth = 0.75,
      position = point_position,
      orientation = "y"
    ) +
    ggplot2::scale_x_log10(
      breaks = c(0.25, 0.5, 1, 2, 4, 8, 16),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    ggplot2::scale_colour_manual(values = c(
      Supported = semantic_palette[["documented"]],
      `Not supported` = semantic_palette[["secondary"]]
    )) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Odds ratio for being chemically studied",
      y = NULL,
      colour = "Evidence"
    ) +
    theme_manuscript() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  if (has_multiple_models) {
    plot <- plot +
      ggplot2::geom_point(
        ggplot2::aes(shape = Model_label),
        size = 2.8,
        stroke = 0.8,
        fill = semantic_palette[["background"]],
        position = point_position
      ) +
      ggplot2::scale_shape_manual(values = model_shapes) +
      ggplot2::labs(shape = "Framework")
  } else {
    plot <- plot + ggplot2::geom_point(size = 2.8, position = point_position)
  }
  plot
}

coverage_contrast <- data.table(
  Dimension = factor(
    c("Species", "Genera", "Phylogenetic diversity"),
    levels = c("Species", "Genera", "Phylogenetic diversity")
  ),
  Fraction = c(173 / 2582, 75 / 189, observed_pd / full_pd)
)
coverage_contrast[, Label := scales::percent(Fraction, accuracy = 0.1)]
figure_1a <- ggplot2::ggplot(
  coverage_contrast,
  ggplot2::aes(x = Fraction, y = Dimension)
) +
  ggplot2::geom_col(width = 0.58, fill = semantic_palette[["documented"]]) +
  ggplot2::geom_text(
    ggplot2::aes(label = Label),
    hjust = -0.15,
    family = FIGURE_FONT,
    fontface = "bold",
    colour = semantic_palette[["universe"]],
    size = 3.5
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::label_percent(),
    limits = c(0, 0.5),
    expand = ggplot2::expansion(mult = c(0, 0.03))
  ) +
  ggplot2::labs(
    title = "Coverage spans many genera but few species",
    subtitle = paste0(
      "Faith's PD SES = ", sprintf("%.2f", pd_ses),
      "; randomization P = ", sprintf("%.3f", pd_two_sided_p)
    ),
    x = "Fraction of the Arecaceae universe",
    y = NULL
  ) +
  theme_manuscript() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "none",
    plot.title = ggplot2::element_text(size = 10.5),
    plot.subtitle = ggplot2::element_text(size = 8.5),
    plot.margin = ggplot2::margin(4, 12, 4, 8)
  )

map_documented <- world_map[world_map$Coverage_group == "Documented palm region", ]
map_undocumented <- world_map[world_map$Coverage_group == "Undocumented palm region", ]
figure_1b <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = world_map,
    fill = semantic_palette[["zero"]],
    colour = semantic_palette[["background"]],
    linewidth = 0.08
  ) +
  ggplot2::geom_sf(
    data = map_documented,
    ggplot2::aes(fill = Representation_rate),
    colour = semantic_palette[["background"]],
    linewidth = 0.08
  ) +
  ggplot2::geom_sf(
    data = map_undocumented,
    fill = semantic_palette[["gap"]],
    colour = semantic_palette[["background"]],
    linewidth = 0.08
  ) +
  ggplot2::scale_fill_gradient(
    low = semantic_palette[["documented_light"]],
    high = semantic_palette[["documented"]],
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, max(map_documented$Representation_rate, na.rm = TRUE)),
    oob = scales::squish,
    name = "Studied species\nwithin native region"
  ) +
  ggplot2::coord_sf(crs = "+proj=robin", datum = NA, expand = FALSE) +
  ggplot2::labs(
    title = "Geographic evidence remains uneven across native palm regions",
    subtitle = "Orange areas contain native palms but no chemically studied species"
  ) +
  ggplot2::theme_void(base_family = FIGURE_FONT) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", colour = semantic_palette[["universe"]], size = 11),
    plot.subtitle = ggplot2::element_text(colour = semantic_palette[["secondary"]], size = 9),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  )

figure_1 <- patchwork::wrap_plots(
  figure_1a,
  figure_1b,
  ncol = 1,
  heights = c(0.42, 1.58)
) +
  patchwork::plot_annotation(
    title = "Palm phytochemical knowledge spans the family but remains shallow and geographically uneven",
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = FIGURE_FONT,
        face = "bold",
        size = 16,
        colour = semantic_palette[["universe"]]
      ),
      plot.tag = ggplot2::element_text(
        family = FIGURE_FONT,
        face = "bold",
        size = 13,
        colour = semantic_palette[["universe"]]
      )
    )
  )

figure_2 <- blank_itol_panel(
  "Phylogenetic, geographic and taxonomic distribution of phytochemical research coverage across Arecaceae",
  "Reserved for the genus-level iTOL figure assembled in Affinity"
)

combined_selection <- data.table::copy(main_selection_models[
  Model %in% c("M3_logistic_combined", "M3_phylo_combined") & Term != "(Intercept)"
])
combined_selection[, `:=`(
  Display_term = unname(term_labels[Term]),
  Model_label = ifelse(Model == "M3_phylo_combined", "Phylogenetic", "Taxonomic"),
  Supported = as.integer(P_value < 0.05)
)]
combined_selection <- combined_selection[!is.na(Display_term)]
figure_3a <- forest_panel(
  combined_selection,
  "Selection is structured by biological visibility and data availability",
  "Parallel taxonomic and phylogenetic models separate robust from framework-dependent effects",
  model_shapes = c(Taxonomic = 16, Phylogenetic = 17)
)

trait_selection <- data.table::copy(trait_models)
trait_selection[, `:=`(
  Display_term = Trait,
  Supported = to_binary(Holm_significant)
)]
figure_3b <- forest_panel(
  trait_selection,
  "Individual traits reveal selective attention",
  "Holm correction controls the family of eight trait tests"
)

socioeconomic_selection <- data.table::copy(socioeconomic_models[
  Term %in% c(
    "Commodity_species_resolved",
    "Human_use_any",
    "Economic_quantitative_data_available"
  )
])
socioeconomic_selection[, `:=`(
  Display_term = unname(term_labels[Term]),
  Supported = as.integer(P_value < 0.05)
)]
figure_3c <- forest_panel(
  socioeconomic_selection,
  "Human use is associated with research attention",
  "Bias-reduced models are shown because positive socioeconomic classes are sparse"
)

figure_3 <- patchwork::wrap_plots(
  figure_3a,
  figure_3b,
  figure_3c,
  nrow = 1,
  widths = c(1.25, 1, 0.95)
) +
  patchwork::plot_annotation(
    title = "Palm phytochemistry is a selected sample rather than a neutral survey of diversity",
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = FIGURE_FONT,
        face = "bold",
        size = 16,
        colour = semantic_palette[["universe"]]
      ),
      plot.tag = ggplot2::element_text(
        family = FIGURE_FONT,
        face = "bold",
        size = 13,
        colour = semantic_palette[["universe"]]
      )
    )
  )

domain_plot_data <- data.table::copy(domain_architecture)
domain_plot_data[, `:=`(
  Chemical_entity_count = to_numeric(Chemical_entity_count),
  Exact_structure_count = to_numeric(Exact_structure_count),
  Domain_analysis = factor(
    Domain_analysis,
    levels = Domain_analysis[order(Chemical_entity_count)]
  )
)]
domain_plot_long <- data.table::melt(
  domain_plot_data,
  id.vars = "Domain_analysis",
  measure.vars = c("Chemical_entity_count", "Exact_structure_count"),
  variable.name = "Metric",
  value.name = "Count"
)
domain_plot_long[, Metric := factor(
  Metric,
  levels = c("Chemical_entity_count", "Exact_structure_count"),
  labels = c("Curated entities", "Exact structures")
)]
figure_4a <- ggplot2::ggplot(
  domain_plot_long,
  ggplot2::aes(x = Count, y = Domain_analysis, fill = Domain_analysis, alpha = Metric)
) +
  ggplot2::geom_col(position = "identity", width = 0.72) +
  ggplot2::scale_fill_manual(values = domain_palette, guide = "none") +
  ggplot2::scale_alpha_manual(values = c(`Curated entities` = 0.36, `Exact structures` = 1)) +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    title = "Chemical domains differ in both breadth and structural resolution",
    x = "Canonical chemical entities",
    y = NULL,
    alpha = "Evidence resolution"
  ) +
  theme_manuscript() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

class_plot_data <- data.table::copy(
  class_architecture[Metabolite_scope == 1L][order(-Chemical_entity_count)][1:20]
)
class_plot_data[, Class_analysis := factor(
  Class_analysis,
  levels = rev(Class_analysis)
)]
figure_4b <- ggplot2::ggplot(
  class_plot_data,
  ggplot2::aes(
    x = Genus_count,
    y = Class_analysis,
    size = Chemical_entity_count,
    fill = Domain_analysis
  )
) +
  ggplot2::geom_segment(
    data = class_plot_data,
    ggplot2::aes(
      x = 0,
      xend = Genus_count,
      y = Class_analysis,
      yend = Class_analysis
    ),
    inherit.aes = FALSE,
    colour = semantic_palette[["grid"]],
    linewidth = 0.4,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    shape = 21,
    colour = semantic_palette[["background"]],
    stroke = 0.5,
    alpha = 0.95
  ) +
  ggplot2::scale_fill_manual(values = domain_palette) +
  ggplot2::scale_size_continuous(
    range = c(2.5, 11),
    breaks = c(25, 50, 100, 200, 400),
    name = "Chemical entities"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 50, by = 10),
    limits = c(0, 50),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::scale_y_discrete(labels = label_wrap(28)) +
  ggplot2::labs(
    title = "Abundant chemical classes are not necessarily taxonomically broad",
    subtitle = "The twenty richest resolved classes are shown",
    x = "Genera with evidence",
    y = NULL,
    fill = "Chemical domain"
  ) +
  theme_manuscript() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "right"
  )

heatmap_genera <- intersect(
  genus_coverage[Studied_species >= 2L, Genus_current],
  genus_tree$tip.label
)
if (length(heatmap_genera) < 20L) {
  stop("Too few genera have at least two chemically studied species.", call. = FALSE)
}
heatmap_classes <- class_architecture[
  Metabolite_scope == 1L & !is.na(Class_analysis)
][order(-Chemical_entity_count), head(Class_analysis, 12L)]
genus_class_grid <- data.table::CJ(
  Genus_current = heatmap_genera,
  Class_analysis = heatmap_classes,
  unique = TRUE
)
genus_class_grid <- merge(
  genus_class_grid,
  genus_class[, .(
    Genus_current,
    Class_analysis,
    Species_prevalence_within_genus = to_numeric(Species_prevalence_within_genus)
  )],
  by = c("Genus_current", "Class_analysis"),
  all.x = TRUE,
  sort = FALSE
)
genus_class_grid[is.na(Species_prevalence_within_genus), Species_prevalence_within_genus := 0]
genus_class_matrix_data <- data.table::dcast(
  genus_class_grid,
  Genus_current ~ Class_analysis,
  value.var = "Species_prevalence_within_genus",
  fill = 0
)
genus_class_matrix <- as.matrix(genus_class_matrix_data[, -"Genus_current"])
rownames(genus_class_matrix) <- genus_class_matrix_data$Genus_current

chemical_profile <- genus_class_matrix
chemical_profile_row_sums <- rowSums(chemical_profile)
positive_profile_rows <- chemical_profile_row_sums > 0
chemical_profile[positive_profile_rows, ] <- sqrt(
  chemical_profile[positive_profile_rows, , drop = FALSE] /
    chemical_profile_row_sums[positive_profile_rows]
)
chemical_clustering <- stats::hclust(
  stats::dist(chemical_profile, method = "euclidean"),
  method = "ward.D2"
)
chemical_profile_tree <- ape::as.phylo(chemical_clustering)
ape::write.tree(
  chemical_profile_tree,
  file.path(ITOL_DIRECTORY, "06_chemical_profile_dendrogram.nwk")
)

heatmap_class_key <- class_architecture[
  match(heatmap_classes, Class_analysis),
  .(Class_analysis, Domain_analysis)
]
heatmap_class_key[, Color := unname(domain_palette[Domain_analysis])]
if (any(is.na(heatmap_class_key$Color))) {
  stop("A chemical domain color is missing from the iTOL heatmap.", call. = FALSE)
}

itol_heatmap_matrix <- apply(
  genus_class_matrix,
  2,
  function(values) sprintf("%.6f", values)
)
itol_heatmap_lines <- c(
  "DATASET_HEATMAP",
  "SEPARATOR TAB",
  "DATASET_LABEL\tChemical class prevalence within genera",
  paste0("COLOR\t", semantic_palette[["documented"]]),
  paste("FIELD_LABELS", paste(heatmap_classes, collapse = "\t"), sep = "\t"),
  paste("FIELD_COLORS", paste(heatmap_class_key$Color, collapse = "\t"), sep = "\t"),
  paste0("COLOR_MIN\t", semantic_palette[["zero"]]),
  "USE_MID_COLOR\t1",
  paste0("COLOR_MID\t", semantic_palette[["documented_light"]]),
  paste0("COLOR_MAX\t", semantic_palette[["universe"]]),
  "USER_MIN_VALUE\t0",
  "USER_MID_VALUE\t0.5",
  "USER_MAX_VALUE\t1",
  "MARGIN\t5",
  "BORDER_WIDTH\t0",
  "SHOW_LABELS\t0",
  "DATA",
  apply(
    cbind(rownames(genus_class_matrix), itol_heatmap_matrix),
    1,
    paste,
    collapse = "\t"
  )
)
write_itol(itol_heatmap_lines, "07_itol_genus_class_heatmap.txt")

chemical_profile_export <- data.table::as.data.table(
  genus_class_matrix,
  keep.rownames = "Genus_current"
)
write_semicolon(
  chemical_profile_export,
  file.path(ITOL_DIRECTORY, "08_chemical_profile_matrix.csv")
)
write_semicolon(
  heatmap_class_key,
  file.path(ITOL_DIRECTORY, "09_chemical_class_key.csv")
)

itol_validation <- c(
  "ITOL PACKAGE VALIDATION - ARECACEAE",
  paste("Genus phylogeny tips:", ape::Ntip(genus_tree)),
  paste("Genus richness records:", nrow(itol_genus_data)),
  paste("Genera with chemical evidence:", sum(itol_genus_data$Has_studies)),
  paste("Genera without chemical evidence:", sum(itol_genus_data$Has_studies == 0)),
  paste("Chemical-profile dendrogram tips:", ape::Ntip(chemical_profile_tree)),
  paste("Chemical classes in heatmap:", ncol(genus_class_matrix)),
  paste("Genera in heatmap:", nrow(genus_class_matrix)),
  "Genus phylogeny match: 184/184",
  "Chemical heatmap values: species prevalence within genus",
  "Chemical clustering: Hellinger transformation, Euclidean distance, Ward.D2"
)
write_itol(itol_validation, "VALIDATION.txt")

figure_4c <- blank_itol_panel(
  "Clustered genus-level chemical profiles",
  "Reserved for the chemical dendrogram and heatmap assembled in iTOL and Affinity"
)

figure_4 <- patchwork::wrap_plots(
  figure_4a,
  figure_4b,
  figure_4c,
  design = "AB\nCC",
  heights = c(0.88, 1.6)
) +
  patchwork::plot_annotation(
    title = "Palm phytochemistry is organized into dominant domains and clustered class profiles",
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = FIGURE_FONT,
        face = "bold",
        size = 16,
        colour = semantic_palette[["universe"]]
      ),
      plot.tag = ggplot2::element_text(
        family = FIGURE_FONT,
        face = "bold",
        size = 13,
        colour = semantic_palette[["universe"]]
      )
    )
  )

organ_plot_data <- data.table::melt(
  data.table::copy(organ_architecture),
  id.vars = "Organ_primary",
  measure.vars = c("Article_count", "Chemical_entity_count"),
  variable.name = "Metric",
  value.name = "Count"
)
organ_order <- organ_architecture[order(Article_count), Organ_primary]
organ_plot_data[, `:=`(
  Organ_primary = factor(Organ_primary, levels = organ_order),
  Metric = factor(
    Metric,
    levels = c("Article_count", "Chemical_entity_count"),
    labels = c("Articles", "Chemical entities")
  )
)]
figure_5a <- ggplot2::ggplot(
  organ_plot_data,
  ggplot2::aes(x = Count, y = Organ_primary, fill = Metric)
) +
  ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
  ggplot2::facet_wrap(~Metric, scales = "free_x", nrow = 1) +
  ggplot2::scale_fill_manual(values = c(
    Articles = semantic_palette[["universe"]],
    `Chemical entities` = semantic_palette[["documented"]]
  )) +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    title = "Sampling effort and documented richness concentrate in a few organs",
    x = "Count",
    y = NULL
  ) +
  theme_manuscript() +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

leading_organs <- organ_architecture[order(-Article_count), head(Organ_primary, 10L)]
leading_organ_classes <- class_architecture[
  Metabolite_scope == 1L & !is.na(Class_analysis)
][order(-Article_count), head(Class_analysis, 15L)]
organ_class_grid <- data.table::CJ(
  Organ_primary = leading_organs,
  Class_analysis = leading_organ_classes,
  unique = TRUE
)
organ_class_grid <- merge(
  organ_class_grid,
  organ_class[, .(
    Organ_primary,
    Class_analysis,
    Species_prevalence_within_organ = to_numeric(Species_prevalence_within_organ)
  )],
  by = c("Organ_primary", "Class_analysis"),
  all.x = TRUE,
  sort = FALSE
)
organ_class_grid[is.na(Species_prevalence_within_organ), Species_prevalence_within_organ := 0]
organ_class_grid[, `:=`(
  Organ_primary = factor(Organ_primary, levels = rev(leading_organs)),
  Class_analysis = factor(Class_analysis, levels = leading_organ_classes)
)]
figure_5b <- ggplot2::ggplot(
  organ_class_grid,
  ggplot2::aes(x = Class_analysis, y = Organ_primary, fill = Species_prevalence_within_organ)
) +
  ggplot2::geom_tile(
    colour = semantic_palette[["background"]],
    linewidth = 0.35
  ) +
  ggplot2::scale_fill_gradientn(
    colours = c(
      semantic_palette[["zero"]],
      semantic_palette[["documented_light"]],
      semantic_palette[["documented"]],
      semantic_palette[["universe"]]
    ),
    limits = c(0, 1),
    labels = scales::label_percent(),
    name = "Species prevalence\nwithin organ"
  ) +
  ggplot2::scale_x_discrete(labels = label_wrap(18)) +
  ggplot2::labs(
    title = "Organ choice filters which chemical classes become visible",
    subtitle = "Ten most sampled organs and fifteen most documented classes",
    x = NULL,
    y = NULL
  ) +
  theme_manuscript() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 48, hjust = 1, vjust = 1, size = 7),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right"
  )

method_plot_data <- data.table::copy(method_class_models)
method_plot_data[, `:=`(
  Log2_odds_ratio = pmax(-5, pmin(5, log2(Odds_ratio))),
  Analytical_family = factor(
    Analytical_family,
    levels = c("GC-based", "LC-based", "Isolation/NMR", "Mass spectrometry")
  ),
  Class_analysis = factor(
    Class_analysis,
    levels = rev(leading_classes$Class_analysis)
  )
)]
figure_5c <- ggplot2::ggplot(
  method_plot_data,
  ggplot2::aes(x = Analytical_family, y = Class_analysis, fill = Log2_odds_ratio)
) +
  ggplot2::geom_tile(
    colour = semantic_palette[["background"]],
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    data = method_plot_data[Supported == 1L],
    shape = 8,
    size = 2.2,
    colour = semantic_palette[["text"]]
  ) +
  ggplot2::scale_fill_gradient2(
    low = semantic_palette[["gap"]],
    mid = semantic_palette[["zero"]],
    high = semantic_palette[["documented"]],
    midpoint = 0,
    limits = c(-5, 5),
    oob = scales::squish,
    name = "Method association\nlog2 odds ratio"
  ) +
  ggplot2::scale_y_discrete(labels = label_wrap(25)) +
  ggplot2::labs(
    title = "Analytical methods systematically expose different chemical classes",
    subtitle = "Stars indicate Holm-supported associations in adjusted article-level models",
    x = NULL,
    y = NULL
  ) +
  theme_manuscript() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    panel.grid = ggplot2::element_blank(),
    legend.position = "right"
  )

article_log_mean <- mean(log1p(species_effort$Article_count))
article_log_sd <- stats::sd(log1p(species_effort$Article_count))
prediction_grid <- data.table(
  Article_count = exp(seq(
    log(min(species_effort$Article_count)),
    log(max(species_effort$Article_count)),
    length.out = 240L
  ))
)
prediction_grid[, `:=`(
  Article_z = (log1p(Article_count) - article_log_mean) / article_log_sd,
  Organ_z = 0,
  Method_z = 0
)]
effort_prediction <- stats::predict(
  effort_models$entities,
  newdata = prediction_grid,
  type = "link",
  se.fit = TRUE
)
prediction_grid[, `:=`(
  Expected = exp(effort_prediction$fit),
  CI_low = exp(effort_prediction$fit - 1.96 * effort_prediction$se.fit),
  CI_high = exp(effort_prediction$fit + 1.96 * effort_prediction$se.fit)
)]
species_effort[, Evidence_tier := cut(
  Article_count,
  breaks = c(-Inf, 1, 3, 7, Inf),
  labels = c("Minimal: 1 article", "Limited: 2-3 articles", "Moderate: 4-7 articles", "Broad: 8+ articles")
)]
figure_5d <- ggplot2::ggplot(
  species_effort,
  ggplot2::aes(x = Article_count, y = Entity_richness_observed)
) +
  ggplot2::geom_ribbon(
    data = prediction_grid,
    ggplot2::aes(x = Article_count, ymin = CI_low, ymax = CI_high),
    inherit.aes = FALSE,
    fill = semantic_palette[["documented_light"]],
    alpha = 0.8
  ) +
  ggplot2::geom_line(
    data = prediction_grid,
    ggplot2::aes(x = Article_count, y = Expected),
    inherit.aes = FALSE,
    colour = semantic_palette[["documented"]],
    linewidth = 0.95
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = Evidence_tier),
    shape = 21,
    size = 2.35,
    alpha = 0.82,
    colour = semantic_palette[["background"]],
    stroke = 0.35
  ) +
  ggplot2::scale_x_log10(
    breaks = c(1, 2, 3, 5, 10, 20, 40),
    labels = scales::label_number()
  ) +
  ggplot2::scale_fill_manual(values = c(
    `Minimal: 1 article` = semantic_palette[["missing"]],
    `Limited: 2-3 articles` = semantic_palette[["documented_light"]],
    `Moderate: 4-7 articles` = semantic_palette[["documented"]],
    `Broad: 8+ articles` = semantic_palette[["universe"]]
  )) +
  ggplot2::labs(
    title = "Documented richness rises strongly with bibliographic effort",
    subtitle = "Negative-binomial expectation at average organ and method breadth",
    x = "Articles contributing chemical evidence per species",
    y = "Observed canonical chemical entities",
    fill = "Evidence tier"
  ) +
  theme_manuscript()

figure_5 <- patchwork::wrap_plots(
  figure_5a,
  figure_5b,
  figure_5c,
  figure_5d,
  design = "AB\nCD",
  heights = c(0.95, 1.05)
) +
  patchwork::plot_annotation(
    title = "Organ choice, analytical method and publication effort filter the chemistry that is documented",
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = FIGURE_FONT,
        face = "bold",
        size = 16,
        colour = semantic_palette[["universe"]]
      ),
      plot.tag = ggplot2::element_text(
        family = FIGURE_FONT,
        face = "bold",
        size = 13,
        colour = semantic_palette[["universe"]]
      )
    )
  )

screening_flow <- data.table(
  Stage = factor(
    c(
      "Screened documents",
      "Included articles",
      "Original evidence records",
      "Analytical incidences",
      "Canonical entities",
      "Biological-scope entities"
    ),
    levels = rev(c(
      "Screened documents",
      "Included articles",
      "Original evidence records",
      "Analytical incidences",
      "Canonical entities",
      "Biological-scope entities"
    ))
  ),
  Count = c(389, 250, 6287, 5207, 1677, 1645),
  X = 1,
  Fill = c("Universe", rep("Documented", 5))
)
screening_flow[, Label := paste0(Stage, "\n", scales::comma(Count))]
figure_s1 <- ggplot2::ggplot(screening_flow, ggplot2::aes(x = X, y = Stage)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 1, xend = 1, y = Stage, yend = data.table::shift(Stage, type = "lead")),
    colour = semantic_palette[["grid"]],
    linewidth = 2,
    na.rm = TRUE
  ) +
  ggplot2::geom_label(
    ggplot2::aes(label = Label, fill = Fill),
    family = FIGURE_FONT,
    fontface = "bold",
    colour = semantic_palette[["background"]],
    size = 4,
    label.size = 0,
    label.padding = grid::unit(0.55, "lines")
  ) +
  ggplot2::annotate(
    "label",
    x = 1.62,
    y = 5.55,
    label = "139 records excluded",
    family = FIGURE_FONT,
    size = 3.5,
    fill = semantic_palette[["gap"]],
    colour = semantic_palette[["background"]],
    label.size = 0
  ) +
  ggplot2::scale_fill_manual(values = c(
    Universe = semantic_palette[["universe"]],
    Documented = semantic_palette[["documented"]]
  )) +
  ggplot2::coord_cartesian(xlim = c(0.45, 2.05), clip = "off") +
  ggplot2::labs(
    title = "Evidence screening, consolidation and chemical standardization",
    subtitle = "Counts refer to distinct analytical levels and should not be interpreted as a single attrition series"
  ) +
  ggplot2::theme_void(base_family = FIGURE_FONT) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14, colour = semantic_palette[["universe"]]),
    plot.subtitle = ggplot2::element_text(size = 10, colour = semantic_palette[["secondary"]]),
    legend.position = "none",
    plot.margin = ggplot2::margin(12, 60, 12, 12)
  )

deviation_plot_data <- rbindlist(list(
  species_effort[order(Entity_conditional_score)][1:15],
  species_effort[order(-Entity_conditional_score)][1:15]
), use.names = TRUE, fill = TRUE)
deviation_plot_data <- unique(deviation_plot_data, by = "WCVP_taxon_id")
deviation_plot_data[, `:=`(
  Species_current = factor(
    Species_current,
    levels = Species_current[order(Entity_conditional_score)]
  ),
  Direction = factor(
    ifelse(Entity_conditional_score >= 0, "Richer than expected", "Poorer than expected"),
    levels = c("Richer than expected", "Poorer than expected")
  )
)]
figure_s2 <- ggplot2::ggplot(
  deviation_plot_data,
  ggplot2::aes(x = Entity_conditional_score, y = Species_current, fill = Direction)
) +
  ggplot2::geom_vline(
    xintercept = 0,
    colour = semantic_palette[["secondary"]],
    linewidth = 0.45
  ) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::scale_fill_manual(values = c(
    `Richer than expected` = semantic_palette[["documented"]],
    `Poorer than expected` = semantic_palette[["gap"]]
  )) +
  ggplot2::labs(
    title = "Species at the extremes of effort-conditioned chemical richness",
    subtitle = "The fifteen strongest positive and negative deviations are shown",
    x = "Effort-conditioned richness score, log2 observed/expected",
    y = NULL,
    fill = NULL
  ) +
  theme_manuscript() +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(face = "italic")
  )

sensitivity_panel <- function(x_column, x_label, rho) {
  ggplot2::ggplot(
    species_effort,
    ggplot2::aes(x = .data[[x_column]], y = Entity_conditional_score)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      colour = semantic_palette[["grid"]]
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      colour = semantic_palette[["grid"]]
    ) +
    ggplot2::geom_point(
      shape = 21,
      size = 2.2,
      fill = semantic_palette[["documented"]],
      colour = semantic_palette[["background"]],
      stroke = 0.35,
      alpha = 0.78
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      colour = semantic_palette[["universe"]],
      fill = semantic_palette[["documented_light"]],
      linewidth = 0.75
    ) +
    ggplot2::annotate(
      "label",
      x = -Inf,
      y = Inf,
      hjust = -0.08,
      vjust = 1.2,
      label = paste0("Spearman rho = ", sprintf("%.2f", rho)),
      family = FIGURE_FONT,
      size = 3.4,
      colour = semantic_palette[["universe"]],
      fill = semantic_palette[["background"]],
      label.size = 0
    ) +
    ggplot2::labs(
      x = x_label,
      y = "Entity-based conditional score"
    ) +
    theme_manuscript()
}

figure_s3a <- sensitivity_panel(
  "Exact_conditional_score",
  "Exact-structure conditional score",
  sensitivity_correlations[Comparison == "Entity versus exact-structure conditional score", Spearman_rho]
)
figure_s3b <- sensitivity_panel(
  "Strict_conditional_score",
  "Strict-entity conditional score",
  sensitivity_correlations[Comparison == "Entity versus strict-entity conditional score", Spearman_rho]
)
figure_s3c <- sensitivity_panel(
  "Class_conditional_score",
  "Chemical-class conditional score",
  sensitivity_correlations[Comparison == "Entity versus chemical-class conditional score", Spearman_rho]
)
figure_s3 <- patchwork::wrap_plots(
  figure_s3a,
  figure_s3b,
  figure_s3c,
  nrow = 1
) +
  patchwork::plot_annotation(
    title = "Effort-conditioned rankings are robust to alternative definitions of chemical richness",
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = FIGURE_FONT,
        face = "bold",
        size = 15,
        colour = semantic_palette[["universe"]]
      ),
      plot.tag = ggplot2::element_text(
        family = FIGURE_FONT,
        face = "bold",
        size = 13,
        colour = semantic_palette[["universe"]]
      )
    )
  )

rarefaction_plot_data <- data.table::copy(
  species_effort[Article_count >= 2L & is.finite(Two_article_expected)]
)
rarefaction_plot_data <- rarefaction_plot_data[order(Two_article_expected)]
rarefaction_plot_data[, Species_rank := seq_len(.N)]
figure_s4 <- ggplot2::ggplot(
  rarefaction_plot_data,
  ggplot2::aes(x = Species_rank, y = Two_article_expected)
) +
  ggplot2::geom_linerange(
    ggplot2::aes(ymin = Two_article_low, ymax = Two_article_high),
    colour = semantic_palette[["documented_light"]],
    linewidth = 0.7,
    alpha = 0.85
  ) +
  ggplot2::geom_point(
    shape = 21,
    size = 2.5,
    fill = semantic_palette[["documented"]],
    colour = semantic_palette[["background"]],
    stroke = 0.35
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = Entity_richness_observed),
    shape = 1,
    size = 2.1,
    colour = semantic_palette[["universe"]],
    stroke = 0.65
  ) +
  ggplot2::labs(
    title = "Two-article richness estimates retain substantial within-species uncertainty",
    subtitle = "Filled points and intervals show the mean and 95% pairwise range; open points show full observed richness",
    x = "Species ordered by mean two-article richness",
    y = "Canonical chemical entities"
  ) +
  theme_manuscript() +
  ggplot2::theme(legend.position = "none")

export_manifest <- rbindlist(list(
  save_vector_figure(
    figure_1,
    MAIN_FIGURE_DIRECTORY,
    "Figure_01_taxonomic_and_geographic_coverage",
    16,
    9.5
  ),
  save_vector_figure(
    figure_2,
    MAIN_FIGURE_DIRECTORY,
    "Figure_02_phylogenetic_research_coverage",
    13.2,
    10.1
  ),
  save_vector_figure(
    figure_3,
    MAIN_FIGURE_DIRECTORY,
    "Figure_03_research_selection",
    16,
    7.8
  ),
  save_vector_figure(
    figure_4,
    MAIN_FIGURE_DIRECTORY,
    "Figure_04_chemical_architecture",
    16,
    14
  ),
  save_vector_figure(
    figure_5,
    MAIN_FIGURE_DIRECTORY,
    "Figure_05_observation_filters",
    16,
    12
  ),
  save_vector_figure(
    figure_s1,
    SUPPLEMENTARY_FIGURE_DIRECTORY,
    "Figure_S01_evidence_screening_flow",
    9,
    10
  ),
  save_vector_figure(
    figure_s2,
    SUPPLEMENTARY_FIGURE_DIRECTORY,
    "Figure_S02_effort_conditioned_deviations",
    10,
    11
  ),
  save_vector_figure(
    figure_s3,
    SUPPLEMENTARY_FIGURE_DIRECTORY,
    "Figure_S03_effort_score_sensitivity",
    15,
    5.8
  ),
  save_vector_figure(
    figure_s4,
    SUPPLEMENTARY_FIGURE_DIRECTORY,
    "Figure_S04_two_article_rarefaction",
    11,
    7
  )
), use.names = TRUE, fill = TRUE)

write_semicolon(export_manifest, file.path(TABLE_DIRECTORY, "figure_export_manifest.csv"))

normalized_project_directory <- paste0(
  normalizePath(PROJECT_DIRECTORY, winslash = "/", mustWork = FALSE),
  "/"
)
normalized_input_paths <- normalizePath(input_paths, winslash = "/", mustWork = FALSE)
required_input_manifest <- data.table(
  Input = names(input_paths),
  Repository_path = sub(normalized_project_directory, "", normalized_input_paths, fixed = TRUE),
  Exists = file.exists(input_paths),
  Bytes = file.info(input_paths)$size
)
write_semicolon(
  required_input_manifest,
  file.path(TABLE_DIRECTORY, "required_input_manifest.csv")
)

expected_itol_files <- c(
  "01_arecaceae_184_genera_cladogram.nwk",
  "02_itol_genus_richness_logbar.txt",
  "03_itol_study_coverage_percent.txt",
  "04_itol_native_continents_6.txt",
  "05_genus_level_data.csv",
  "06_chemical_profile_dendrogram.nwk",
  "07_itol_genus_class_heatmap.txt",
  "08_chemical_profile_matrix.csv",
  "09_chemical_class_key.csv",
  "S1_itol_study_presence.txt",
  "S2_itol_native_TDWG_regions_8.txt",
  "S3_arecaceae_184_genera_astral_branch_lengths.nwk",
  "VALIDATION.txt"
)
itol_export_manifest <- data.table(
  File = expected_itol_files,
  Path = file.path(ITOL_DIRECTORY, expected_itol_files)
)
itol_export_manifest[, `:=`(
  Exists = file.exists(Path),
  Bytes = file.info(Path)$size
)]
write_semicolon(
  itol_export_manifest,
  file.path(TABLE_DIRECTORY, "Table_iTOL_export_manifest.csv")
)

pipeline_audit <- data.table(
  Test = c(
    "Normalized analytical incidences",
    "Included articles",
    "Species universe",
    "Studied species",
    "Species-tree tips",
    "Studied species in species tree",
    "Genus-tree tips",
    "Method-by-class coefficients",
    "Chemical-profile dendrogram tips",
    "Main figures",
    "Supplementary figures",
    "iTOL files",
    "All iTOL files are non-empty",
    "All PDF exports are non-empty",
    "All SVG exports are non-empty",
    "All PNG exports are non-empty"
  ),
  Observed = c(
    nrow(incidence),
    uniqueN(articles$Article_ID),
    nrow(species),
    uniqueN(main_incidence$WCVP_taxon_id),
    ape::Ntip(species_tree),
    length(studied_tree_tips),
    ape::Ntip(genus_tree),
    nrow(method_class_models),
    ape::Ntip(chemical_profile_tree),
    sum(grepl("^Figure_0", export_manifest$Figure)),
    sum(grepl("^Figure_S", export_manifest$Figure)),
    sum(itol_export_manifest$Exists),
    as.integer(
      all(itol_export_manifest$Exists) &&
        all(itol_export_manifest$Bytes[itol_export_manifest$Exists] > 0)
    ),
    as.integer(all(export_manifest$PDF_bytes > 0)),
    as.integer(all(export_manifest$SVG_bytes > 0)),
    as.integer(all(export_manifest$PNG_bytes > 0))
  ),
  Expected = c(
    5207, 250, 2582, 173, 2233, 172, 184, 48,
    length(heatmap_genera), 5, 4, 13, 1, 1, 1, 1
  )
)
pipeline_audit[, Status := ifelse(Observed == Expected, "PASS", "FAIL")]
write_semicolon(pipeline_audit, file.path(TABLE_DIRECTORY, "Table_pipeline_QA.csv"))

if (any(pipeline_audit$Status != "PASS")) {
  stop("The final pipeline audit detected one or more failures.", call. = FALSE)
}

capture.output(
  utils::sessionInfo(),
  file = file.path(LOG_DIRECTORY, "figure_generation_sessionInfo.txt")
)

message("[05/05] Complete. All analyses, tables and editable figures passed the final audit.")
