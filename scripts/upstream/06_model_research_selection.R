# =============================================================================
# ARTIGO ARECACEAE — ETAPA 05B
# Vies taxonomico, geografico, funcional e socioeconomico na selecao de
# especies para estudos quimicos
# Versao: 1.0.3 | 2026-08-13
#
# PRINCIPIOS
# - A resposta e binaria: especie estudada ou nao estudada quimicamente.
# - O universo analitico contem as 2.582 especies WCVP atuais.
# - Especies nao estudadas nao recebem riqueza quimica zero; metricas quimicas
#   sao mascaradas como NA na copia processada desta etapa.
# - A amplitude primaria e o numero de regioes nativas POWO/WGSRPD nivel 3,
#   nao abundancia. A amplitude WCVP total e apenas sensibilidade.
# - Traits ausentes nao sao imputados automaticamente.
# - Uso humano e importancia economica sao indicadores de vies de pesquisa,
#   nao causas demonstradas.
# - Commodity e economia possuem separacao; por isso usam regressao logistica
#   com reducao de vies e tabelas 2 x 2 explicitas.
# =============================================================================

get_script_path <- function() {
  configured_path <- Sys.getenv("ARECACEAE_SCRIPT_PATH", unset = "")
  if (nzchar(configured_path) && file.exists(configured_path)) {
    return(normalizePath(configured_path, winslash = "/"))
  }
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(path)) return(normalizePath(path, winslash = "/"))
  }
  fallback <- file.path(
    getwd(), "05B_analisar_vies_selecao_pesquisa_v1_0_0.R"
  )
  if (file.exists(fallback)) return(normalizePath(fallback, winslash = "/"))
  stop(
    "Nao foi possivel identificar o arquivo do script. Execute com Rscript, ",
    "pelo RStudio ou defina ARECACEAE_SCRIPT_PATH.",
    call. = FALSE
  )
}

SCRIPT_PATH <- get_script_path()
SCRIPT_DIR <- dirname(SCRIPT_PATH)
PROJECT_DIR <- Sys.getenv(
  "ARECACEAE_PROJECT_DIR",
  unset = "G:/Meu Drive/ALUNOS/Mestrado/Ana/review/Arquivos_organizados"
)

required_packages <- c(
  "data.table", "ape", "ggplot2", "phylolm", "sandwich",
  "brglm2", "patchwork", "scales", "svglite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Instale os pacotes ausentes antes de continuar:\ninstall.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "), "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ape)
  library(ggplot2)
  library(patchwork)
})

if (!capabilities("cairo")) {
  stop("O R instalado nao possui suporte Cairo.", call. = FALSE)
}

figure_font <- Sys.getenv(
  "ARECACEAE_FIGURE_FONT",
  unset = if (.Platform$OS.type == "windows") "Arial" else "sans"
)

first_existing <- function(paths, label) {
  found <- paths[file.exists(paths)]
  if (length(found) == 0L) {
    stop(
      label, " nao encontrado. Caminhos testados:\n- ",
      paste(paths, collapse = "\n- "),
      call. = FALSE
    )
  }
  normalizePath(found[1L], winslash = "/")
}

analysis_candidates <- c(
  file.path(PROJECT_DIR, "09_article_analysis", "article_analysis_v1_0_2"),
  file.path(PROJECT_DIR, "09_article_analysis"),
  dirname(SCRIPT_DIR),
  SCRIPT_DIR
)

SPECIES_FILE <- first_existing(
  c(
    file.path(
      analysis_candidates, "data_processed",
      "02_species_analysis_with_GBIF_GHS.csv"
    )
  ),
  "Base integrada de especies com GBIF-GHS"
)

TREE_FILE <- first_existing(
  c(
    file.path(
      analysis_candidates, "data_processed",
      "01_species_tree_WCVP_pruned.tre"
    )
  ),
  "Arvore WCVP podada"
)

OUT_DIR <- file.path(
  PROJECT_DIR, "09_article_analysis",
  "stage05B_research_selection_bias"
)
DATA_DIR <- file.path(OUT_DIR, "data_processed")
TABLES_DIR <- file.path(OUT_DIR, "tables")
FIGURES_DIR <- file.path(OUT_DIR, "figures")
MODELS_DIR <- file.path(OUT_DIR, "models")
LOGS_DIR <- file.path(OUT_DIR, "logs")
invisible(lapply(
  c(OUT_DIR, DATA_DIR, TABLES_DIR, FIGURES_DIR, MODELS_DIR, LOGS_DIR),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

SCRIPT_COPY_PATH <- file.path(OUT_DIR, basename(SCRIPT_PATH))
script_copy_ok <- if (identical(
  normalizePath(SCRIPT_PATH, winslash = "/", mustWork = FALSE),
  normalizePath(SCRIPT_COPY_PATH, winslash = "/", mustWork = FALSE)
)) {
  TRUE
} else {
  isTRUE(file.copy(
    SCRIPT_PATH, SCRIPT_COPY_PATH, overwrite = TRUE, copy.mode = TRUE
  ))
}
if (!script_copy_ok) {
  stop("Nao foi possivel incluir o script na pasta de resultados.", call. = FALSE)
}

read_semicolon <- function(path) {
  data.table::fread(
    path, sep = ";", encoding = "UTF-8",
    na.strings = c("", "NA"), check.names = FALSE
  )
}

write_semicolon <- function(x, path) {
  data.table::fwrite(
    as.data.table(x), path, sep = ";", bom = TRUE, quote = TRUE, na = ""
  )
  invisible(path)
}

trim_to_na <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x == "NA"] <- NA_character_
  x
}

to_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

z_score <- function(x) {
  x <- to_numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) < 2L) return(out)
  sx <- stats::sd(x[ok])
  if (!is.finite(sx) || sx == 0) return(out)
  out[ok] <- (x[ok] - mean(x[ok])) / sx
  out
}

wilson_interval <- function(successes, total, confidence = 0.95) {
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  p <- successes / total
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half <- z * sqrt(
    p * (1 - p) / total + z^2 / (4 * total^2)
  ) / denominator
  data.table(
    CI_low = pmax(0, center - half),
    CI_high = pmin(1, center + half)
  )
}

started_at <- Sys.time()
set.seed(20260813)

message("[01/10] Lendo e auditando a base integrada...")

species <- read_semicolon(SPECIES_FILE)
species_tree <- ape::read.tree(TREE_FILE)

required_columns <- c(
  "WCVP_taxon_id", "Analysis_species_key", "Species_current",
  "Genus_current", "Research_studied_binary", "POWO_native_L3_count",
  "POWO_native_continents", "WCVP_TDWG_count", "Trait_data_available",
  "Tree_included", "PalmSubfamily", "PalmTribe", "Climbing",
  "Acaulescent", "Erect", "StemSolitary", "StemArmed", "LeavesArmed",
  "UnderstoreyCanopy", "MaxStemHeight_m", "MaxStemDia_cm",
  "Max_Blade_Length_m", "AverageFruitLength_cm",
  "AverageFruitWidth_cm", "GBIF_GHS_data_available",
  "GBIF_unique_10km_cells", "GHS_POP_2025_mean_log1p",
  "Commodity_species_resolved", "Commodity_candidate_broad",
  "Human_use_any", "Economic_quantitative_data_available"
)
missing_columns <- setdiff(required_columns, names(species))
if (length(missing_columns) > 0L) {
  stop(
    "A base integrada nao contem: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

numeric_columns <- c(
  "Research_studied_binary", "POWO_native_L3_count", "WCVP_TDWG_count",
  "Trait_data_available", "Tree_included", "Climbing", "Acaulescent",
  "Erect", "StemSolitary", "StemArmed", "LeavesArmed",
  "MaxStemHeight_m", "MaxStemDia_cm", "Max_Blade_Length_m",
  "AverageFruitLength_cm", "AverageFruitWidth_cm",
  "GBIF_GHS_data_available", "GBIF_unique_10km_cells",
  "GHS_POP_2025_mean_log1p", "Commodity_species_resolved",
  "Commodity_candidate_broad", "Human_use_any",
  "Economic_quantitative_data_available"
)
numeric_columns <- intersect(numeric_columns, names(species))
species[, (numeric_columns) := lapply(.SD, to_numeric), .SDcols = numeric_columns]
species[, `:=`(
  WCVP_taxon_id = sub("\\.0$", "", trim_to_na(WCVP_taxon_id)),
  Analysis_species_key = trim_to_na(Analysis_species_key),
  Species_current = trim_to_na(Species_current),
  Genus_current = trim_to_na(Genus_current),
  PalmSubfamily = trim_to_na(PalmSubfamily),
  PalmTribe = trim_to_na(PalmTribe),
  POWO_native_continents = trim_to_na(POWO_native_continents),
  UnderstoreyCanopy = tolower(trim_to_na(UnderstoreyCanopy))
)]

# PalmSubfamily foi originalmente anexada com PalmTraits e, por isso, sua
# ausencia nao deve ser interpretada como ausencia taxonomica. Primeiro usamos
# a correspondencia genero-subfamilia inequivoca observada na propria base.
# Os nove generos sem referencia interna recebem uma classificacao explicita
# segundo a literatura taxonomica de Arecaceae. Referencias: Genera Palmarum 2
# (Dransfield et al., 2008); Baker & Dransfield (2016), DOI 10.1111/boj.12401;
# Bernal & Galeano (2013), DOI 10.11646/phytotaxa.144.2.1; e Sam et al.
# (2023), DOI 10.11646/phytotaxa.613.3.1. Isso nao imputa traits.
species[, PalmSubfamily_source := fifelse(
  !is.na(PalmSubfamily), "PalmTraits", NA_character_
)]
observed_genus_subfamily <- species[
  !is.na(Genus_current) & !is.na(PalmSubfamily),
  .(
    N_subfamilies = uniqueN(PalmSubfamily),
    Subfamily = sort(unique(PalmSubfamily))[1L]
  ),
  by = Genus_current
][N_subfamilies == 1L]
genus_subfamily_map <- setNames(
  observed_genus_subfamily$Subfamily,
  observed_genus_subfamily$Genus_current
)
subfamily_from_genus <- unname(genus_subfamily_map[species$Genus_current])
fill_from_genus <- is.na(species$PalmSubfamily) &
  !is.na(subfamily_from_genus)
species[fill_from_genus, `:=`(
  PalmSubfamily = subfamily_from_genus[fill_from_genus],
  PalmSubfamily_source = "Genus mapping from PalmTraits taxonomy"
)]

taxonomic_subfamily_fallback <- c(
  Acanthorhiza = "Coryphoideae",
  Butyagrus = "Arecoideae",
  Ceratolobus = "Calamoideae",
  Klopstockia = "Arecoideae",
  Martinezia = "Arecoideae",
  Nephrosperma = "Arecoideae",
  Sabinaria = "Coryphoideae",
  Truongsonia = "Arecoideae",
  Wallaceodoxa = "Arecoideae"
)
taxonomic_subfamily_fallback_table <- data.table(
  Genus_current = names(taxonomic_subfamily_fallback),
  PalmSubfamily = unname(taxonomic_subfamily_fallback),
  Reference = c(
    rep(
      paste(
        "Dransfield et al. (2008), Genera Palmarum 2;",
        "Baker & Dransfield (2016), DOI 10.1111/boj.12401"
      ),
      6L
    ),
    "Bernal & Galeano (2013), DOI 10.11646/phytotaxa.144.2.1",
    "Sam et al. (2023), DOI 10.11646/phytotaxa.613.3.1",
    "Heatubun et al. (2014), DOI 10.1007/s12225-014-9525-x"
  )
)
subfamily_from_literature <- unname(
  taxonomic_subfamily_fallback[species$Genus_current]
)
fill_from_literature <- is.na(species$PalmSubfamily) &
  !is.na(subfamily_from_literature)
species[fill_from_literature, `:=`(
  PalmSubfamily = subfamily_from_literature[fill_from_literature],
  PalmSubfamily_source = "Arecaceae taxonomic literature fallback"
)]

if (nrow(species) != 2582L || uniqueN(species$WCVP_taxon_id) != 2582L) {
  stop("Esperadas 2.582 especies WCVP unicas.", call. = FALSE)
}
if (sum(species$Research_studied_binary == 1L, na.rm = TRUE) != 173L) {
  stop("Esperadas 173 especies estudadas.", call. = FALSE)
}
if (!inherits(species_tree, "phylo") || ape::Ntip(species_tree) != 2233L) {
  stop("Esperada arvore WCVP com 2.233 terminais.", call. = FALSE)
}

# O arquivo legado preencheu zeros estruturais para especies nao estudadas.
# Eles nao entram nos modelos de selecao e sao convertidos para NA somente na
# copia processada desta etapa.
legacy_chemistry_columns <- intersect(c(
  "Research_article_count", "Research_evidence_record_count",
  "Chem_richness_entities_observed", "Chem_richness_exact_structures",
  "Chem_richness_classes", "Chem_richness_domains",
  "Chem_organ_primary_count", "Article_count", "Organ_primary_count",
  "Occurrence_or_presence_records", "Chemical_entity_count",
  "Exact_structure_count", "Chemical_domain_count", "Chemical_class_count",
  "Compounds_per_article"
), names(species))
species[
  Research_studied_binary == 0L,
  (legacy_chemistry_columns) := lapply(.SD, function(x) NA),
  .SDcols = legacy_chemistry_columns
]

message("[02/10] Construindo eixos funcionais e preditores padronizados...")

build_size_axis <- function(data, columns, axis_name, labels) {
  raw <- as.matrix(data[, lapply(.SD, function(x) log1p(to_numeric(x))),
                        .SDcols = columns])
  complete <- stats::complete.cases(raw)
  scaled <- scale(raw[complete, , drop = FALSE])
  pca <- stats::prcomp(scaled, center = FALSE, scale. = FALSE)
  loading <- pca$rotation[, 1L]
  score <- pca$x[, 1L]
  if (sum(loading) < 0) {
    loading <- -loading
    score <- -score
  }
  full_score <- rep(NA_real_, nrow(data))
  full_score[complete] <- score
  full_score <- z_score(full_score)
  list(
    score = full_score,
    pca = pca,
    loadings = data.table(
      Axis = axis_name,
      Component = labels,
      Original_variable = columns,
      PC1_loading = unname(loading),
      PC1_variance_explained = unname(summary(pca)$importance[2L, 1L]),
      Species_complete = sum(complete)
    )
  )
}

stature_axis <- build_size_axis(
  species,
  c("MaxStemHeight_m", "MaxStemDia_cm"),
  "Palm stature",
  c("Maximum stem height", "Maximum stem diameter")
)
fruit_axis <- build_size_axis(
  species,
  c("AverageFruitLength_cm", "AverageFruitWidth_cm"),
  "Fruit size",
  c("Average fruit length", "Average fruit width")
)

species[, `:=`(
  z_Palm_stature_PC1 = stature_axis$score,
  z_Fruit_size_PC1 = fruit_axis$score,
  z_log_native_range_L3 = z_score(log1p(POWO_native_L3_count)),
  z_log_WCVP_TDWG_sensitivity = z_score(log1p(WCVP_TDWG_count)),
  z_log_GBIF_10km_cells = z_score(fifelse(
    GBIF_GHS_data_available == 1L,
    log1p(GBIF_unique_10km_cells),
    NA_real_
  )),
  z_GHS_mean_exposure = z_score(fifelse(
    GBIF_GHS_data_available == 1L,
    GHS_POP_2025_mean_log1p,
    NA_real_
  )),
  Trait_leaf_blade_z = z_score(log1p(Max_Blade_Length_m)),
  Trait_climbing = fifelse(
    is.na(Climbing), NA_real_, as.numeric(Climbing %in% c(1, 2))
  ),
  Trait_acaulescent = fifelse(
    is.na(Acaulescent), NA_real_, as.numeric(Acaulescent %in% c(1, 2))
  ),
  Trait_erect = fifelse(
    is.na(Erect), NA_real_, as.numeric(Erect %in% c(1, 2))
  ),
  Trait_solitary = fifelse(
    is.na(StemSolitary), NA_real_, as.numeric(StemSolitary %in% c(1, 2))
  ),
  Trait_stem_armed = fifelse(
    is.na(StemArmed), NA_real_, as.numeric(StemArmed %in% c(1, 2))
  ),
  Trait_leaves_armed = fifelse(
    is.na(LeavesArmed), NA_real_, as.numeric(LeavesArmed %in% c(1, 2))
  ),
  Trait_canopy = fifelse(
    UnderstoreyCanopy == "canopy", 1,
    fifelse(UnderstoreyCanopy == "understorey", 0, NA_real_)
  )
)]

continent_map <- c(
  "AFRICA" = "Native_Africa",
  "ASIA-TEMPERATE" = "Native_Asia_temperate",
  "ASIA-TROPICAL" = "Native_Asia_tropical",
  "AUSTRALASIA" = "Native_Australasia",
  "EUROPE" = "Native_Europe",
  "NORTHERN AMERICA" = "Native_Northern_America",
  "PACIFIC" = "Native_Pacific",
  "SOUTHERN AMERICA" = "Native_Southern_America"
)
for (continent in names(continent_map)) {
  target <- unname(continent_map[[continent]])
  species[, (target) := as.integer(
    !is.na(POWO_native_continents) &
      grepl(continent, POWO_native_continents, fixed = TRUE)
  )]
}

axis_loadings <- rbindlist(
  list(stature_axis$loadings, fruit_axis$loadings),
  use.names = TRUE
)
axis_loadings[, Interpretation := fifelse(
  Axis == "Palm stature",
  "Higher values indicate taller palms with thicker stems.",
  "Higher values indicate longer and wider fruits."
)]

message("[03/10] Resumindo cobertura taxonomica e geografica...")

summarise_coverage <- function(data, group_column) {
  out <- data[!is.na(get(group_column)), .(
    Studied_species = sum(Research_studied_binary == 1L),
    Total_species = .N
  ), by = group_column]
  out[, Representation_rate := Studied_species / Total_species]
  intervals <- wilson_interval(out$Studied_species, out$Total_species)
  cbind(out, intervals)
}

subfamily_coverage <- summarise_coverage(species, "PalmSubfamily")
setnames(subfamily_coverage, "PalmSubfamily", "Group")
subfamily_coverage[, Coverage_type := "Palm subfamily"]

continent_long <- species[
  !is.na(POWO_native_continents),
  .(Group = trimws(unlist(strsplit(
    POWO_native_continents, "\\|"
  )))),
  by = .(WCVP_taxon_id, Research_studied_binary)
]
continent_long <- unique(
  continent_long,
  by = c("WCVP_taxon_id", "Group")
)
continent_coverage <- continent_long[, .(
  Studied_species = sum(Research_studied_binary == 1L),
  Total_species = .N
), by = Group]
continent_coverage[, Representation_rate := Studied_species / Total_species]
continent_coverage <- cbind(
  continent_coverage,
  wilson_interval(
    continent_coverage$Studied_species,
    continent_coverage$Total_species
  )
)
continent_coverage[, Coverage_type := "Native continent"]

genus_coverage <- species[, .(
  Studied_species = sum(Research_studied_binary == 1L),
  Total_species = .N,
  Representation_rate = mean(Research_studied_binary == 1L)
), by = Genus_current][order(-Representation_rate, -Studied_species)]

continent_screen <- continent_coverage[, .(
  Group,
  Eligible_for_model = as.integer(
    Total_species >= 50L & Studied_species >= 5L
  ),
  Total_species,
  Studied_species
)]
continent_screen[, Predictor := unname(continent_map[Group])]
eligible_continent_predictors <- continent_screen[
  Eligible_for_model == 1L & !is.na(Predictor),
  Predictor
]

message("[04/10] Ajustando modelos logisticos e filogeneticos principais...")

auc_rank <- function(y, probability) {
  ok <- is.finite(y) & is.finite(probability)
  y <- y[ok]
  probability <- probability[ok]
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  (sum(rank(probability)[y == 1L]) - n1 * (n1 + 1) / 2) /
    (n1 * n0)
}

tidy_coefficients <- function(
  estimates, standard_errors, p_values, model_name, framework,
  model_class, formula_text
) {
  out <- data.table(
    Term = names(estimates),
    Estimate = as.numeric(estimates),
    Std_error = as.numeric(standard_errors),
    Statistic = as.numeric(estimates / standard_errors),
    P_value = as.numeric(p_values)
  )
  out[, `:=`(
    Model = model_name,
    Framework = framework,
    Model_class = model_class,
    Formula = formula_text,
    Odds_ratio = exp(Estimate),
    CI_low = exp(Estimate - 1.96 * Std_error),
    CI_high = exp(Estimate + 1.96 * Std_error)
  )]
  out
}

fit_binomial_glm_with_retry <- function(formula, data) {
  controls <- list(
    stats::glm.control(epsilon = 1e-8, maxit = 200),
    stats::glm.control(epsilon = 1e-8, maxit = 1000)
  )
  fit <- NULL
  for (attempt in seq_along(controls)) {
    fit <- stats::glm(
      formula,
      data = data,
      family = stats::binomial(),
      na.action = stats::na.omit,
      model = TRUE,
      y = TRUE,
      x = TRUE,
      control = controls[[attempt]]
    )
    if (isTRUE(fit$converged)) break
  }
  attr(fit, "fit_attempts") <- attempt
  fit
}

fit_logistic_cluster <- function(formula, data, model_name, framework) {
  model_data <- as.data.frame(copy(data))
  rownames(model_data) <- seq_len(nrow(model_data))
  fit <- fit_binomial_glm_with_retry(formula, model_data)
  used_rows <- as.integer(rownames(stats::model.frame(fit)))
  cluster <- model_data$Genus_current[used_rows]
  robust_vcov <- sandwich::vcovCL(
    fit, cluster = cluster, type = "HC1"
  )
  estimates <- stats::coef(fit)
  standard_errors <- sqrt(diag(robust_vcov))
  p_values <- 2 * stats::pnorm(
    abs(estimates / standard_errors), lower.tail = FALSE
  )
  formula_text <- paste(deparse(formula), collapse = " ")
  coefficients <- tidy_coefficients(
    estimates, standard_errors, p_values,
    model_name, framework,
    "Logistic; genus-clustered HC1 SE",
    formula_text
  )
  fitted_probability <- stats::fitted(fit)
  fit_summary <- data.table(
    Model = model_name,
    Framework = framework,
    Model_class = "Logistic; genus-clustered HC1 SE",
    N_species = stats::nobs(fit),
    N_events = sum(fit$y == 1L),
    AIC = stats::AIC(fit),
    AUC_in_sample = auc_rank(fit$y, fitted_probability),
    Brier_score = mean((fit$y - fitted_probability)^2),
    Phylogenetic_alpha = NA_real_,
    Converged = isTRUE(fit$converged),
    Iterations = as.integer(fit$iter),
    Fit_attempts = as.integer(attr(fit, "fit_attempts")),
    Optimizer_code = NA_integer_,
    Alpha_warning = NA_integer_,
    AIC_comparison_note = paste(
      "Compare AIC only among models fitted to identical species sets."
    )
  )
  list(
    model = fit,
    robust_vcov = robust_vcov,
    coefficients = coefficients,
    fit = fit_summary
  )
}

fit_phylogenetic_logistic <- function(
  formula, data, phy, model_name, framework
) {
  variables <- all.vars(formula)
  model_data <- data[
    Tree_included == 1L &
      !is.na(Analysis_species_key) &
      complete.cases(data[, ..variables])
  ]
  model_data <- unique(model_data, by = "Analysis_species_key")
  keep <- intersect(phy$tip.label, model_data$Analysis_species_key)
  model_tree <- ape::keep.tip(phy, keep)
  model_data <- model_data[match(
    model_tree$tip.label, Analysis_species_key
  )]
  model_frame <- as.data.frame(model_data)
  rownames(model_frame) <- model_frame$Analysis_species_key
  fit <- phylolm::phyloglm(
    formula,
    data = model_frame,
    phy = model_tree,
    method = "logistic_MPLE",
    btol = 10
  )
  coefficient_matrix <- summary(fit)$coefficients
  estimates <- coefficient_matrix[, 1L]
  standard_errors <- coefficient_matrix[, 2L]
  p_values <- coefficient_matrix[, 4L]
  coefficients <- tidy_coefficients(
    estimates, standard_errors, p_values,
    model_name, framework,
    "Phylogenetic logistic MPLE",
    paste(deparse(formula), collapse = " ")
  )
  fit_summary <- data.table(
    Model = model_name,
    Framework = framework,
    Model_class = "Phylogenetic logistic MPLE",
    N_species = nrow(model_data),
    N_events = sum(model_data$Research_studied_binary == 1L),
    AIC = tryCatch(stats::AIC(fit), error = function(e) NA_real_),
    AUC_in_sample = NA_real_,
    Brier_score = NA_real_,
    Phylogenetic_alpha = if (!is.null(fit$alpha)) fit$alpha else NA_real_,
    Converged = isTRUE(fit$convergence == 0L),
    Iterations = NA_integer_,
    Fit_attempts = 1L,
    Optimizer_code = if (!is.null(fit$convergence)) {
      as.integer(fit$convergence)
    } else {
      NA_integer_
    },
    Alpha_warning = if (!is.null(fit$alphaWarn)) {
      as.integer(fit$alphaWarn)
    } else {
      NA_integer_
    },
    AIC_comparison_note = paste(
      "Compare AIC only among models fitted to identical species sets."
    )
  )
  list(
    model = fit,
    tree = model_tree,
    coefficients = coefficients,
    fit = fit_summary
  )
}

formula_geography <- stats::as.formula(paste(
  "Research_studied_binary ~ z_log_native_range_L3 +",
  paste(eligible_continent_predictors, collapse = " + ")
))
formula_morphology <- Research_studied_binary ~
  z_log_native_range_L3 + z_Palm_stature_PC1 + z_Fruit_size_PC1
formula_accessibility <- Research_studied_binary ~
  z_log_native_range_L3 + z_log_GBIF_10km_cells + z_GHS_mean_exposure
formula_combined <- Research_studied_binary ~
  z_log_native_range_L3 + z_Palm_stature_PC1 + z_Fruit_size_PC1 +
  z_log_GBIF_10km_cells + z_GHS_mean_exposure

main_specs <- list(
  M0_logistic_geography = list(
    formula_geography, "Geographic availability", FALSE
  ),
  M1_logistic_morphology = list(
    formula_morphology, "Morphology", FALSE
  ),
  M1_phylo_morphology = list(
    formula_morphology, "Morphology", TRUE
  ),
  M2_logistic_accessibility = list(
    formula_accessibility, "Accessibility and exposure", FALSE
  ),
  M2_phylo_accessibility = list(
    formula_accessibility, "Accessibility and exposure", TRUE
  ),
  M3_logistic_combined = list(
    formula_combined, "Combined sensitivity", FALSE
  ),
  M3_phylo_combined = list(
    formula_combined, "Combined sensitivity", TRUE
  )
)

main_models <- list()
main_coefficients <- list()
main_fit <- list()
for (model_name in names(main_specs)) {
  specification <- main_specs[[model_name]]
  fitted <- if (isTRUE(specification[[3L]])) {
    fit_phylogenetic_logistic(
      specification[[1L]], species, species_tree,
      model_name, specification[[2L]]
    )
  } else {
    fit_logistic_cluster(
      specification[[1L]], species,
      model_name, specification[[2L]]
    )
  }
  main_models[[model_name]] <- fitted$model
  main_coefficients[[model_name]] <- fitted$coefficients
  main_fit[[model_name]] <- fitted$fit
}
main_coefficients <- rbindlist(
  main_coefficients, use.names = TRUE, fill = TRUE
)
main_fit <- rbindlist(main_fit, use.names = TRUE, fill = TRUE)

# Sensibilidade da amplitude: os dois indicadores sao comparados com as mesmas
# covariaveis continentais e exatamente o mesmo conjunto de especies. A
# distribuicao WCVP total pode incluir amplitude introduzida e nao substitui o
# indicador nativo POWO.
range_sensitivity_data <- species[
  complete.cases(z_log_native_range_L3, z_log_WCVP_TDWG_sensitivity)
]
formula_range_native_matched <- stats::as.formula(paste(
  "Research_studied_binary ~ z_log_native_range_L3 +",
  paste(eligible_continent_predictors, collapse = " + ")
))
formula_range_wcvp_matched <- stats::as.formula(paste(
  "Research_studied_binary ~ z_log_WCVP_TDWG_sensitivity +",
  paste(eligible_continent_predictors, collapse = " + ")
))
range_sensitivity_specs <- list(
  RANGE_native_L3_matched = formula_range_native_matched,
  RANGE_WCVP_total_matched = formula_range_wcvp_matched
)
range_sensitivity_models <- list()
range_sensitivity_coefficients <- list()
range_sensitivity_fit <- list()
for (model_name in names(range_sensitivity_specs)) {
  fitted <- fit_logistic_cluster(
    range_sensitivity_specs[[model_name]],
    range_sensitivity_data,
    model_name,
    "Range sensitivity; matched sample"
  )
  range_sensitivity_models[[model_name]] <- fitted$model
  range_sensitivity_coefficients[[model_name]] <- fitted$coefficients
  range_sensitivity_fit[[model_name]] <- fitted$fit
}
range_sensitivity_coefficients <- rbindlist(
  range_sensitivity_coefficients, use.names = TRUE, fill = TRUE
)
range_sensitivity_fit <- rbindlist(
  range_sensitivity_fit, use.names = TRUE, fill = TRUE
)
range_sensitivity_coefficients <- merge(
  range_sensitivity_coefficients,
  range_sensitivity_fit[, .(Model, N_species, N_events)],
  by = "Model", all.x = TRUE, sort = FALSE
)

message("[05/10] Ajustando modelos individuais de traits...")

trait_specifications <- data.table(
  Trait = c(
    "Maximum leaf-blade length",
    "Climbing habit",
    "Acaulescent habit",
    "Erect habit",
    "Solitary stem",
    "Stem armature",
    "Leaf armature",
    "Canopy position"
  ),
  Predictor = c(
    "Trait_leaf_blade_z",
    "Trait_climbing",
    "Trait_acaulescent",
    "Trait_erect",
    "Trait_solitary",
    "Trait_stem_armed",
    "Trait_leaves_armed",
    "Trait_canopy"
  ),
  Scale = c(
    "Per 1 SD increase in log1p-transformed length",
    rep("Present or variable versus absent", 6L),
    "Canopy versus understorey"
  )
)

trait_models <- list()
trait_results <- list()
trait_fit <- list()
for (i in seq_len(nrow(trait_specifications))) {
  trait_name <- trait_specifications$Trait[i]
  predictor <- trait_specifications$Predictor[i]
  formula <- stats::as.formula(paste(
    "Research_studied_binary ~ z_log_native_range_L3 +",
    predictor
  ))
  model_name <- paste0("T", sprintf("%02d", i), "_", predictor)
  fitted <- fit_logistic_cluster(
    formula, species, model_name, "Trait-specific"
  )
  result <- fitted$coefficients[Term == predictor]
  result[, `:=`(
    Trait = trait_name,
    Trait_scale = trait_specifications$Scale[i],
    N_species = fitted$fit$N_species,
    N_events = fitted$fit$N_events
  )]
  trait_models[[model_name]] <- fitted$model
  trait_results[[model_name]] <- result
  trait_fit[[model_name]] <- fitted$fit
}
trait_results <- rbindlist(trait_results, use.names = TRUE, fill = TRUE)
trait_fit <- rbindlist(trait_fit, use.names = TRUE, fill = TRUE)
trait_results[, P_value_Holm := stats::p.adjust(P_value, method = "holm")]
trait_results[, Holm_significant := as.integer(P_value_Holm < 0.05)]

message("[06/10] Tratando uso humano e economia com reducao de vies...")

fit_bias_reduced <- function(formula, data, model_name, label) {
  # Para regressao logistica, AS_mean e MPL_Jeffreys produzem o mesmo ajuste
  # de Firth. A ultima tentativa usa a parametrizacao penalizada equivalente,
  # que pode ser numericamente mais estavel sob separacao completa.
  controls <- list(
    brglm2::brglmControl(
      epsilon = 1e-8, maxit = 500, type = "AS_mean", slowit = 1
    ),
    brglm2::brglmControl(
      epsilon = 1e-8, maxit = 2000, type = "AS_mean", slowit = 0.5,
      max_step_factor = 24
    ),
    brglm2::brglmControl(
      epsilon = 1e-8, maxit = 5000, type = "AS_mean", slowit = 0.1,
      max_step_factor = 48
    ),
    brglm2::brglmControl(
      epsilon = 1e-8, maxit = 5000, type = "MPL_Jeffreys", slowit = 0.5,
      max_step_factor = 48, a = 0.5
    )
  )
  estimation_types <- c(
    "AS_mean", "AS_mean", "AS_mean", "MPL_Jeffreys_equivalent"
  )
  # Coeficientes iniciais nulos evitam que a rotina tente usar como ponto de
  # partida uma estimativa de maxima verossimilhanca divergente sob separacao.
  model_data <- as.data.frame(data)
  model_frame <- stats::model.frame(
    formula,
    data = model_data,
    na.action = stats::na.omit
  )
  start_values <- rep(
    0,
    ncol(stats::model.matrix(formula, data = model_frame))
  )
  fit <- NULL
  for (attempt in seq_along(controls)) {
    fit <- stats::glm(
      formula,
      data = model_data,
      family = stats::binomial("logit"),
      method = brglm2::brglmFit,
      na.action = stats::na.omit,
      model = TRUE,
      y = TRUE,
      start = start_values,
      control = controls[[attempt]]
    )
    if (isTRUE(fit$converged)) break
  }
  attr(fit, "fit_attempts") <- attempt
  attr(fit, "estimation_type") <- estimation_types[attempt]
  coefficient_matrix <- summary(fit)$coefficients
  estimates <- coefficient_matrix[, 1L]
  standard_errors <- coefficient_matrix[, 2L]
  p_values <- coefficient_matrix[, 4L]
  coefficients <- tidy_coefficients(
    estimates, standard_errors, p_values,
    model_name, "Human use and economy",
    "Mean bias-reduced logistic",
    paste(deparse(formula), collapse = " ")
  )
  coefficients[, Socioeconomic_definition := label]
  fit_summary <- data.table(
    Model = model_name,
    Framework = "Human use and economy",
    Model_class = "Mean bias-reduced logistic",
    N_species = stats::nobs(fit),
    N_events = sum(fit$y == 1L),
    AIC = stats::AIC(fit),
    AUC_in_sample = auc_rank(fit$y, stats::fitted(fit)),
    Brier_score = mean((fit$y - stats::fitted(fit))^2),
    Phylogenetic_alpha = NA_real_,
    Converged = isTRUE(fit$converged),
    Iterations = as.integer(fit$iter),
    Fit_attempts = as.integer(attr(fit, "fit_attempts")),
    Bias_reduction_type = as.character(attr(fit, "estimation_type")),
    Optimizer_code = NA_integer_,
    Alpha_warning = NA_integer_,
    AIC_comparison_note = paste(
      "Compare AIC only among models fitted to identical species sets."
    )
  )
  list(model = fit, coefficients = coefficients, fit = fit_summary)
}

socio_specs <- list(
  S1_commodity_resolved = list(
    Research_studied_binary ~
      z_log_native_range_L3 + Commodity_species_resolved,
    "Species-resolved commodity flag",
    "Commodity_species_resolved"
  ),
  S2_human_use_any = list(
    Research_studied_binary ~
      z_log_native_range_L3 + Human_use_any,
    "Resolved or broad human-use candidate",
    "Human_use_any"
  ),
  S3_economic_quantitative = list(
    Research_studied_binary ~
      z_log_native_range_L3 + Economic_quantitative_data_available,
    "Species with quantitative economic data",
    "Economic_quantitative_data_available"
  )
)

socio_models <- list()
socio_coefficients <- list()
socio_fit <- list()
socio_fisher <- list()
for (model_name in names(socio_specs)) {
  specification <- socio_specs[[model_name]]
  fitted <- fit_bias_reduced(
    specification[[1L]], species, model_name, specification[[2L]]
  )
  socio_models[[model_name]] <- fitted$model
  socio_coefficients[[model_name]] <- fitted$coefficients
  socio_fit[[model_name]] <- fitted$fit

  flag_name <- specification[[3L]]
  contingency <- table(
    factor(species[[flag_name]], levels = c(0, 1)),
    factor(species$Research_studied_binary, levels = c(0, 1))
  )
  fisher <- stats::fisher.test(contingency)
  a <- contingency[2L, 2L]
  b <- contingency[2L, 1L]
  c_value <- contingency[1L, 2L]
  d <- contingency[1L, 1L]
  log_or <- log(
    ((a + 0.5) * (d + 0.5)) /
      ((b + 0.5) * (c_value + 0.5))
  )
  log_or_se <- sqrt(sum(1 / (c(a, b, c_value, d) + 0.5)))
  socio_fisher[[model_name]] <- data.table(
    Model = model_name,
    Definition = specification[[2L]],
    Flag_column = flag_name,
    Flagged_studied = a,
    Flagged_unstudied = b,
    Unflagged_studied = c_value,
    Unflagged_unstudied = d,
    Fisher_odds_ratio = unname(fisher$estimate),
    Fisher_P_value = fisher$p.value,
    Haldane_odds_ratio = exp(log_or),
    Haldane_CI_low = exp(log_or - 1.96 * log_or_se),
    Haldane_CI_high = exp(log_or + 1.96 * log_or_se)
  )
}
socio_coefficients <- rbindlist(
  socio_coefficients, use.names = TRUE, fill = TRUE
)
socio_fit <- rbindlist(socio_fit, use.names = TRUE, fill = TRUE)
socio_fisher <- rbindlist(socio_fisher, use.names = TRUE, fill = TRUE)

message("[07/10] Construindo diagnosticos e controles de qualidade...")

all_model_fit <- rbindlist(list(
  main_fit,
  range_sensitivity_fit,
  trait_fit,
  socio_fit
), use.names = TRUE, fill = TRUE)

convergence_audit <- all_model_fit[, .(
  Model,
  Framework,
  Model_class,
  N_species,
  N_events,
  Converged,
  Iterations,
  Fit_attempts,
  Bias_reduction_type,
  Optimizer_code,
  Alpha_warning
)]
failed_convergence_models <- convergence_audit[
  !Converged | is.na(Converged), Model
]
phylogenetic_alpha_warnings <- convergence_audit[
  Model_class == "Phylogenetic logistic MPLE" &
    !is.na(Alpha_warning) & Alpha_warning != 0L,
  Model
]

all_main_coefficients <- rbindlist(list(
  main_coefficients,
  range_sensitivity_coefficients
), use.names = TRUE, fill = TRUE)

correlation_variables <- c(
  "z_log_native_range_L3", "z_Palm_stature_PC1", "z_Fruit_size_PC1",
  "z_log_GBIF_10km_cells", "z_GHS_mean_exposure"
)
predictor_correlations <- as.data.table(
  stats::cor(
    species[, ..correlation_variables],
    use = "pairwise.complete.obs",
    method = "spearman"
  ),
  keep.rownames = "Predictor_1"
)
predictor_correlations <- melt(
  predictor_correlations,
  id.vars = "Predictor_1",
  variable.name = "Predictor_2",
  value.name = "Spearman_rho"
)

qa_row <- function(indicator, observed, expected, note = "") {
  status <- if (is.na(expected)) {
    "INFO"
  } else if (isTRUE(all.equal(
    as.numeric(observed), as.numeric(expected)
  ))) {
    "PASS"
  } else {
    "FAIL"
  }
  data.table(
    Indicator = indicator,
    Observed = observed,
    Expected = expected,
    Status = status,
    Note = note
  )
}

unstudied_chemistry_nonmissing <- species[
  Research_studied_binary == 0L,
  sum(!is.na(unlist(.SD))),
  .SDcols = legacy_chemistry_columns
]
duplicated_current_binomials <- species[
  !is.na(Species_current), .N, by = Species_current
][N > 1L, .N]
range_matched_denominators <- range_sensitivity_fit[
  , uniqueN(paste(N_species, N_events, sep = ":")) == 1L
]

qa <- rbindlist(list(
  qa_row("WCVP species universe", nrow(species), 2582),
  qa_row(
    "Unique WCVP species identifiers",
    uniqueN(species$WCVP_taxon_id), 2582
  ),
  qa_row(
    "Chemically studied species",
    species[Research_studied_binary == 1L, .N], 173
  ),
  qa_row(
    "Unstudied species",
    species[Research_studied_binary == 0L, .N], 2409
  ),
  qa_row(
    "Unstudied species with nonmissing chemical metrics after masking",
    unstudied_chemistry_nonmissing, 0,
    "Unstudied species must not receive chemical richness zero"
  ),
  qa_row("Phylogenetic tree tips", ape::Ntip(species_tree), 2233),
  qa_row(
    "Studied species represented in the tree",
    species[
      Research_studied_binary == 1L & Tree_included == 1L,
      .N
    ],
    172
  ),
  qa_row(
    "Species with PalmTraits",
    species[Trait_data_available == 1L, .N], 2236
  ),
  qa_row(
    "Species with taxonomic subfamily assigned",
    species[!is.na(PalmSubfamily), .N], 2582,
    paste(
      "Independent of functional-trait availability; PalmTraits, genus",
      "mapping and explicit taxonomic fallbacks are distinguished."
    )
  ),
  qa_row(
    "Studied species with taxonomic subfamily assigned",
    species[
      Research_studied_binary == 1L & !is.na(PalmSubfamily), .N
    ],
    173
  ),
  qa_row(
    "Stature-axis complete species",
    sum(!is.na(species$z_Palm_stature_PC1)), 1787
  ),
  qa_row(
    "Fruit-axis complete species",
    sum(!is.na(species$z_Fruit_size_PC1)), 1849
  ),
  qa_row(
    "Species with native GBIF-GHS data",
    species[GBIF_GHS_data_available == 1L, .N], 2198
  ),
  qa_row(
    "Main models fitted", length(main_models), 7
  ),
  qa_row(
    "Matched range-sensitivity models fitted",
    length(range_sensitivity_models), 2
  ),
  qa_row(
    "Matched range-sensitivity denominators identical",
    range_matched_denominators, 1,
    "Both range definitions use the same species and continent covariates"
  ),
  qa_row(
    "Trait-specific models fitted", length(trait_models), 8
  ),
  qa_row(
    "Bias-reduced socioeconomic models fitted",
    length(socio_models), 3
  ),
  qa_row(
    "All fitted models converged, including phylogenetic models",
    length(failed_convergence_models) == 0L, 1,
    if (length(failed_convergence_models) == 0L) {
      "Convergence verified model by model"
    } else {
      paste("Non-converged:", paste(failed_convergence_models, collapse = ", "))
    }
  ),
  qa_row(
    "Phylogenetic alpha-boundary warnings",
    length(phylogenetic_alpha_warnings), NA_real_,
    if (length(phylogenetic_alpha_warnings) == 0L) {
      "No alphaWarn flags returned by phylolm"
    } else {
      paste(
        "Inspect alphaWarn for:",
        paste(phylogenetic_alpha_warnings, collapse = ", ")
      )
    }
  ),
  qa_row(
    "Duplicated current binomials",
    duplicated_current_binomials, NA_real_,
    paste(
      "WCVP taxon identifiers remain the primary keys; the duplicate",
      "combines accepted and provisionally accepted records."
    )
  ),
  qa_row(
    "Analysis script copied to output package",
    file.exists(SCRIPT_COPY_PATH) &&
      identical(
        unname(tools::md5sum(SCRIPT_PATH)),
        unname(tools::md5sum(SCRIPT_COPY_PATH))
      ),
    1
  ),
  qa_row(
    "Resolved commodity flags among unstudied species",
    species[
      Commodity_species_resolved == 1L &
        Research_studied_binary == 0L,
      .N
    ],
    0,
    "Separation is expected and handled by bias reduction"
  ),
  qa_row(
    "Maximum absolute predictor correlation",
    max(abs(predictor_correlations[
      Predictor_1 != Predictor_2,
      Spearman_rho
    ]), na.rm = TRUE),
    NA_real_,
    "Pairwise Spearman diagnostic"
  )
), use.names = TRUE, fill = TRUE)

if (any(qa$Status == "FAIL")) {
  write_semicolon(qa, file.path(TABLES_DIR, "Table_QA_05B.csv"))
  write_semicolon(
    convergence_audit,
    file.path(TABLES_DIR, "Table_QA_05B_model_convergence.csv")
  )
  stop(
    "A auditoria da etapa 05B encontrou falhas. Consulte Table_QA_05B.csv.",
    call. = FALSE
  )
}

message("[08/10] Gravando dados, tabelas e modelos...")

provenance <- data.table(
  Input = c("Integrated species table", "Pruned WCVP tree", "Analysis script"),
  Path = c(SPECIES_FILE, TREE_FILE, SCRIPT_PATH),
  MD5 = unname(tools::md5sum(c(SPECIES_FILE, TREE_FILE, SCRIPT_PATH)))
)

write_semicolon(
  species,
  file.path(DATA_DIR, "05B_species_selection_analysis.csv")
)
write_semicolon(
  axis_loadings,
  file.path(TABLES_DIR, "Table_05B_01_functional_size_axes.csv")
)
write_semicolon(
  subfamily_coverage,
  file.path(TABLES_DIR, "Table_05B_02_subfamily_coverage.csv")
)
write_semicolon(
  continent_coverage,
  file.path(TABLES_DIR, "Table_05B_03_native_continent_coverage.csv")
)
write_semicolon(
  genus_coverage,
  file.path(TABLES_DIR, "Table_05B_04_genus_coverage.csv")
)
write_semicolon(
  all_main_coefficients,
  file.path(TABLES_DIR, "Table_05B_05_main_selection_models.csv")
)
write_semicolon(
  main_fit,
  file.path(TABLES_DIR, "Table_05B_06_main_model_fit.csv")
)
write_semicolon(
  trait_results,
  file.path(TABLES_DIR, "Table_05B_07_trait_specific_models.csv")
)
write_semicolon(
  socio_coefficients,
  file.path(TABLES_DIR, "Table_05B_08_socioeconomic_models.csv")
)
write_semicolon(
  socio_fisher,
  file.path(TABLES_DIR, "Table_05B_09_socioeconomic_2x2_tests.csv")
)
write_semicolon(
  range_sensitivity_coefficients,
  file.path(TABLES_DIR, "Table_05B_10_range_sensitivity.csv")
)
write_semicolon(
  continent_screen,
  file.path(TABLES_DIR, "Table_05B_11_continent_model_screen.csv")
)
write_semicolon(
  predictor_correlations,
  file.path(TABLES_DIR, "Table_05B_12_predictor_correlations.csv")
)
write_semicolon(
  provenance,
  file.path(TABLES_DIR, "Table_05B_13_input_provenance.csv")
)
write_semicolon(
  taxonomic_subfamily_fallback_table,
  file.path(TABLES_DIR, "Table_05B_14_taxonomic_subfamily_fallback.csv")
)
write_semicolon(qa, file.path(TABLES_DIR, "Table_QA_05B.csv"))
write_semicolon(
  convergence_audit,
  file.path(TABLES_DIR, "Table_QA_05B_model_convergence.csv")
)

saveRDS(
  list(
    main_models = main_models,
    trait_models = trait_models,
    socioeconomic_models = socio_models,
    range_sensitivity_models = range_sensitivity_models,
    stature_PCA = stature_axis$pca,
    fruit_PCA = fruit_axis$pca,
    script_version = "1.0.3"
  ),
  file.path(MODELS_DIR, "05B_selection_bias_models.rds"),
  compress = "xz"
)

message("[09/10] Produzindo figuras vetoriais...")

theme_article <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size, base_family = figure_font) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3),
      plot.subtitle = element_text(colour = "#4E5555"),
      plot.caption = element_text(
        colour = "#4E5555", size = base_size - 1, hjust = 0
      ),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(10, 18, 10, 10)
    )
}

save_figure <- function(plot, filename, width, height) {
  pdf_path <- file.path(FIGURES_DIR, paste0(filename, ".pdf"))
  svg_path <- file.path(FIGURES_DIR, paste0(filename, ".svg"))
  png_path <- file.path(FIGURES_DIR, paste0(filename, ".png"))

  grDevices::cairo_pdf(
    pdf_path, width = width, height = height,
    family = figure_font, bg = "white", onefile = FALSE
  )
  print(plot)
  grDevices::dev.off()

  svglite::svglite(
    svg_path, width = width, height = height, bg = "white",
    system_fonts = list(sans = figure_font)
  )
  print(plot)
  grDevices::dev.off()

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      png_path, plot = plot, device = ragg::agg_png,
      width = width, height = height, units = "in", dpi = 600,
      bg = "white"
    )
  } else {
    ggsave(
      png_path, plot = plot,
      width = width, height = height, units = "in", dpi = 600,
      bg = "white"
    )
  }
  data.table(
    Figure = filename,
    PDF_bytes = file.info(pdf_path)$size,
    SVG_bytes = file.info(svg_path)$size,
    PNG_bytes = file.info(png_path)$size
  )
}

coverage_plot <- rbindlist(list(
  subfamily_coverage[, .(
    Coverage_type, Group, Studied_species, Total_species,
    Representation_rate, CI_low, CI_high
  )],
  continent_coverage[, .(
    Coverage_type, Group, Studied_species, Total_species,
    Representation_rate, CI_low, CI_high
  )]
), use.names = TRUE)
coverage_plot[, Label := sprintf(
  "%s/%s", Studied_species, Total_species
)]
coverage_plot[, Group_plot := reorder(
  Group, Representation_rate
)]

figure_05B_01 <- ggplot(
  coverage_plot,
  aes(100 * Representation_rate, Group_plot)
) +
  geom_errorbarh(
    aes(xmin = 100 * CI_low, xmax = 100 * CI_high),
    height = 0.16, colour = "#557A95"
  ) +
  geom_point(
    aes(size = Total_species), colour = "#2F6F8F", alpha = 0.9
  ) +
  geom_text(
    aes(label = Label), nudge_y = 0.22,
    family = figure_font, size = 3, colour = "#4E5555"
  ) +
  facet_wrap(~Coverage_type, scales = "free_y", ncol = 1) +
  scale_x_continuous(
    labels = scales::label_percent(scale = 1, accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_size_area(max_size = 8) +
  labs(
    title = "Chemical-study representation is taxonomically and geographically uneven",
    subtitle = paste(
      "Points are representation rates with Wilson 95% intervals;",
      "labels show studied/available species."
    ),
    caption = paste(
      sprintf(
        "Subfamily coverage: %s species (%s studied); native-continent coverage: %s species (%s studied).",
        species[!is.na(PalmSubfamily), .N],
        species[!is.na(PalmSubfamily) & Research_studied_binary == 1L, .N],
        species[!is.na(POWO_native_continents), .N],
        species[
          !is.na(POWO_native_continents) & Research_studied_binary == 1L,
          .N
        ]
      ),
      sprintf(
        "%s species (%s studied) lack native-continent data.",
        species[is.na(POWO_native_continents), .N],
        species[
          is.na(POWO_native_continents) & Research_studied_binary == 1L,
          .N
        ]
      ),
      "Species may occur on more than one native continent; range is not abundance.",
      sep = "\n"
    ),
    x = "Species represented in the chemical corpus",
    y = NULL,
    size = "Available species"
  ) +
  theme_article() +
  theme(legend.position = "right")

term_labels <- c(
  z_log_native_range_L3 = "Native range breadth",
  z_Palm_stature_PC1 = "Palm stature axis",
  z_Fruit_size_PC1 = "Fruit size axis",
  z_log_GBIF_10km_cells = "Native GBIF 10-km cells",
  z_GHS_mean_exposure = "Human-population exposure"
)
forest <- main_coefficients[
  Framework %chin% c("Morphology", "Accessibility and exposure") &
    Term != "(Intercept)" &
    is.finite(Odds_ratio) & is.finite(CI_low) & is.finite(CI_high)
]
forest[, Predictor := unname(term_labels[Term])]
forest <- forest[!is.na(Predictor)]
forest[, Model_display := fifelse(
  grepl("Phylogenetic", Model_class),
  "Phylogenetic logistic", "Logistic; genus-clustered SE"
)]
forest[, Predictor := factor(
  Predictor,
  levels = rev(c(
    "Native range breadth", "Palm stature axis", "Fruit size axis",
    "Native GBIF 10-km cells", "Human-population exposure"
  ))
)]
model_denominator <- function(model_name) {
  row <- main_fit[Model == model_name]
  sprintf("n = %s; studied = %s", row$N_species, row$N_events)
}
framework_labels <- c(
  Morphology = paste0(
    "Morphology\nLogistic: ",
    model_denominator("M1_logistic_morphology"),
    "; phylogenetic: ",
    model_denominator("M1_phylo_morphology")
  ),
  `Accessibility and exposure` = paste0(
    "Accessibility and exposure\nLogistic: ",
    model_denominator("M2_logistic_accessibility"),
    "; phylogenetic: ",
    model_denominator("M2_phylo_accessibility")
  )
)

figure_05B_02 <- ggplot(
  forest,
  aes(Odds_ratio, Predictor, colour = Model_display)
) +
  geom_vline(
    xintercept = 1, linetype = "dashed", colour = "#737979"
  ) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.15, position = position_dodge(width = 0.45)
  ) +
  geom_point(
    size = 2.8, position = position_dodge(width = 0.45)
  ) +
  facet_wrap(
    ~Framework, scales = "free_y", ncol = 1,
    labeller = as_labeller(framework_labels)
  ) +
  scale_x_log10() +
  scale_colour_manual(values = c(
    "Logistic; genus-clustered SE" = "#2F6F8F",
    "Phylogenetic logistic" = "#B75A4A"
  )) +
  labs(
    title = "Correlates of representation in the chemical corpus",
    subtitle = "Odds ratios per one standard deviation with 95% confidence intervals.",
    caption = paste0(
      "Morphology and accessibility are fitted in separate models because their complete-case denominators differ.\n",
      "Associations are not causal."
    ),
    x = "Odds ratio (log scale)",
    y = NULL,
    colour = "Model"
  ) +
  theme_article() +
  theme(legend.position = "bottom")

trait_plot <- copy(trait_results)
trait_plot[, Trait := factor(
  Trait, levels = Trait[order(Odds_ratio)]
)]
trait_plot[, Holm_result := fifelse(
  Holm_significant == 1L,
  "Holm-adjusted p < 0.05",
  "Holm-adjusted p >= 0.05"
)]

figure_05B_03 <- ggplot(
  trait_plot,
  aes(Odds_ratio, Trait, colour = Holm_result)
) +
  geom_vline(
    xintercept = 1, linetype = "dashed", colour = "#737979"
  ) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.16, linewidth = 0.7
  ) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_colour_manual(values = c(
    "Holm-adjusted p < 0.05" = "#176B66",
    "Holm-adjusted p >= 0.05" = "#A7ADAD"
  )) +
  labs(
    title = "Functional traits associated with chemical-study selection",
    subtitle = paste(
      "Separate logistic models adjust for native range breadth and use",
      "genus-clustered robust standard errors."
    ),
    caption = paste(
      "Continuous leaf length is scaled per SD; binary effects compare",
      "present/variable with absent, and canopy with understorey."
    ),
    x = "Odds ratio (log scale)",
    y = NULL,
    colour = NULL
  ) +
  theme_article() +
  theme(legend.position = "bottom")

socio_plot <- socio_coefficients[
  Term %chin% c(
    "Commodity_species_resolved", "Human_use_any",
    "Economic_quantitative_data_available"
  )
]
socio_labels <- c(
  Commodity_species_resolved = "Resolved commodity species",
  Human_use_any = "Any mapped human use",
  Economic_quantitative_data_available = "Quantitative economic data"
)
socio_plot[, Indicator := unname(socio_labels[Term])]
socio_plot <- merge(
  socio_plot,
  socio_fisher[, .(
    Model, Flagged_studied, Flagged_unstudied
  )],
  by = "Model",
  all.x = TRUE,
  sort = FALSE
)
socio_plot[, Indicator := factor(
  Indicator, levels = rev(socio_labels)
)]

figure_05B_04 <- ggplot(
  socio_plot,
  aes(Odds_ratio, Indicator)
) +
  geom_vline(
    xintercept = 1, linetype = "dashed", colour = "#737979"
  ) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.16, colour = "#B75A4A"
  ) +
  geom_point(size = 3.3, colour = "#B75A4A") +
  geom_text(
    aes(label = paste0(
      Flagged_studied, " studied / ", Flagged_unstudied, " unstudied"
    )),
    nudge_y = 0.18, hjust = 0,
    family = figure_font, size = 3, colour = "#4E5555"
  ) +
  scale_x_log10(expand = expansion(mult = c(0.04, 0.5))) +
  labs(
    title = "Human-use and economic mapping are enriched among studied species",
    subtitle = paste0(
      "Mean bias-reduced logistic models adjust for native range breadth;\n",
      "bias reduction is required because the indicators show separation."
    ),
    caption = paste(
      "Flags represent mapped research-selection correlates, not causal",
      "effects or species-level economic values."
    ),
    x = "Bias-reduced odds ratio (log scale)",
    y = NULL
  ) +
  theme_article()

figure_exports <- rbindlist(list(
  save_figure(
    figure_05B_01,
    "Figure_05B_01_taxonomic_geographic_coverage", 10.3, 8.2
  ),
  save_figure(
    figure_05B_02,
    "Figure_05B_02_main_selection_effects", 9.0, 6.5
  ),
  save_figure(
    figure_05B_03,
    "Figure_05B_03_trait_specific_effects", 9.0, 6.2
  ),
  save_figure(
    figure_05B_04,
    "Figure_05B_04_human_use_economic_enrichment", 9.0, 5.2
  )
), use.names = TRUE)

write_semicolon(
  figure_exports,
  file.path(TABLES_DIR, "Table_QA_05B_figure_exports.csv")
)

message("[10/10] Gravando log e sessionInfo...")

capture.output(
  sessionInfo(),
  file = file.path(LOGS_DIR, "05B_sessionInfo.txt")
)

finished_at <- Sys.time()
log_lines <- c(
  "ETAPA 05B — VIES DE SELECAO PARA ESTUDOS QUIMICOS",
  "Versao: 1.0.3",
  paste0("Inicio: ", format(started_at, "%Y-%m-%d %H:%M:%S")),
  paste0("Fim: ", format(finished_at, "%Y-%m-%d %H:%M:%S")),
  paste0("Duracao_min: ", round(as.numeric(difftime(
    finished_at, started_at, units = "mins"
  )), 2)),
  paste0("Especies_WCVP: ", nrow(species)),
  paste0("Especies_estudadas: ", sum(species$Research_studied_binary == 1L)),
  paste0("Main_models: ", length(main_models)),
  paste0("Range_sensitivity_models: ", length(range_sensitivity_models)),
  paste0("Trait_models: ", length(trait_models)),
  paste0("Socioeconomic_models: ", length(socio_models)),
  paste0("QA_falhas: ", qa[Status == "FAIL", .N]),
  "CONCLUSAO: correlatos de selecao para pesquisa; associacao nao causal"
)
writeLines(
  log_lines,
  file.path(LOGS_DIR, "05B_execution_log.txt"),
  useBytes = TRUE
)

message("ETAPA 05B CONCLUIDA")
message("Resultados salvos em: ", OUT_DIR)
