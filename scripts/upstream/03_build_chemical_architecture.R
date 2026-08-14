# =============================================================================
# ARTIGO ARECACEAE — 04. ARQUITETURA QUIMICA DO CORPUS
# Versao: 1.0.2
#
# Unidades mantidas separadas:
#   Article_ID          = esforco bibliografico
#   Evidence_record_ID  = registro original de presenca
#   Chemical_entity_ID  = unidade de riqueza quimica observada
#   WCVP_taxon_id       = especie corrente
#   Organ_primary       = orgao harmonizado
#
# A etapa cria uma camada de incidencia deduplicada por
# artigo x taxon x orgao x entidade. Os registros originais nao sao alterados.
# =============================================================================

get_script_directory <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  script_file <- grep("^--file=", arguments, value = TRUE)
  if (length(script_file) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", script_file), winslash = "/")))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(path)) return(dirname(normalizePath(path, winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

SCRIPT_DIR <- get_script_directory()
config_candidates <- c(
  file.path(SCRIPT_DIR, "00_configurar_projeto_Arecaceae_v1_0_2.R"),
  file.path(dirname(SCRIPT_DIR), "00_configurar_projeto_Arecaceae_v1_0_2.R")
)
config_file <- config_candidates[file.exists(config_candidates)][1L]
if (is.na(config_file)) {
  message(
    "Configuracao externa nao encontrada; usando o caminho padrao interno."
  )
  options(stringsAsFactors = FALSE, scipen = 999, warn = 1)
  PROJECT_DIR <- Sys.getenv(
    "ARECACEAE_PROJECT_DIR",
    unset = "G:/Meu Drive/ALUNOS/Mestrado/Ana/review/Arquivos_organizados"
  )
  ANALYSIS_DIR <- file.path(PROJECT_DIR, "09_article_analysis")
  DATA_PROCESSED_DIR <- file.path(ANALYSIS_DIR, "data_processed")
  BASE_MASTER_FILE <- file.path(
    PROJECT_DIR, "08_manuscript", "Base_mestra_Arecaceae_250_artigos_v1.xlsx"
  )
  SPECIES_BASE_OUTPUT <- file.path(
    DATA_PROCESSED_DIR, "01_species_analysis_base.csv"
  )
  assert_files_exist <- function(paths, labels = basename(paths)) {
    missing <- !file.exists(paths)
    if (any(missing)) {
      stop(
        "Arquivo(s) ausente(s):\n",
        paste0("- ", labels[missing], ": ", paths[missing], collapse = "\n"),
        call. = FALSE
      )
    }
    invisible(TRUE)
  }
  require_packages <- function(packages) {
    missing <- packages[
      !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing) > 0L) {
      stop(
        "Instale os pacotes ausentes: install.packages(c(",
        paste(sprintf("\"%s\"", missing), collapse = ", "), "))",
        call. = FALSE
      )
    }
    invisible(TRUE)
  }
  normalize_species_name <- function(x) {
    x <- trimws(as.character(x))
    x <- gsub("[[:space:]]+", " ", x)
    x[x == ""] <- NA_character_
    x
  }
  to_numeric_safe <- function(x) suppressWarnings(as.numeric(as.character(x)))
  collapse_sorted <- function(x) {
    x <- sort(unique(trimws(as.character(x))))
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) NA_character_ else paste(x, collapse = " | ")
  }
  write_csv_semicolon <- function(x, path) {
    data.table::fwrite(
      data.table::as.data.table(x), path,
      sep = ";", na = "", bom = TRUE, quote = TRUE
    )
    invisible(path)
  }
  write_session_info <- function(path) {
    capture.output(sessionInfo(), file = path)
    invisible(path)
  }
} else {
  source(config_file)
}

set.seed(20260813)
require_packages(c("data.table", "readxl", "ggplot2", "scales", "patchwork"))
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

STAGE04_DIR <- file.path(ANALYSIS_DIR, "stage04_chemical_architecture")
STAGE04_DATA_DIR <- file.path(STAGE04_DIR, "data_processed")
STAGE04_TABLES_DIR <- file.path(STAGE04_DIR, "tables")
STAGE04_FIGURES_DIR <- file.path(STAGE04_DIR, "figures")
STAGE04_LOGS_DIR <- file.path(STAGE04_DIR, "logs")
invisible(lapply(
  c(STAGE04_DIR, STAGE04_DATA_DIR, STAGE04_TABLES_DIR,
    STAGE04_FIGURES_DIR, STAGE04_LOGS_DIR),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

assert_files_exist(
  c(BASE_MASTER_FILE, SPECIES_BASE_OUTPUT),
  c("base mestra quimica", "tabela analitica de especies da etapa 01")
)

started_at <- Sys.time()
message("[01/07] Lendo e validando as unidades originais...")

articles <- as.data.table(readxl::read_excel(
  BASE_MASTER_FILE, sheet = "Articles_250", .name_repair = "minimal"
))
evidence <- as.data.table(readxl::read_excel(
  BASE_MASTER_FILE, sheet = "Evidence_presence", .name_repair = "minimal"
))
entities <- as.data.table(readxl::read_excel(
  BASE_MASTER_FILE, sheet = "Chemical_entities", .name_repair = "minimal"
))
species_universe <- data.table::fread(
  SPECIES_BASE_OUTPUT, sep = ";", encoding = "UTF-8", na.strings = c("", "NA")
)

required_article_columns <- c("Article_ID", "Year", "Title")
required_evidence_columns <- c(
  "Evidence_record_ID", "Source_record_ID", "Article_ID", "DOI", "Taxon_ID",
  "Taxon_standardized", "Genus", "Species_standardized", "Taxon_type",
  "Species_analysis_ready", "Organ_primary", "Compound_ID",
  "Chemical_entity_ID", "Evidence_grade", "Methods"
)
required_entity_columns <- c(
  "Chemical_entity_ID", "Entity_canonical_name", "Chemical_domain",
  "Chemical_class_harmonized", "Classification_confidence", "Sulfated",
  "Exact_structure_ID", "Structure_status", "Validated_PubChem_CID",
  "Validated_InChIKey", "Molecular_formula", "Molecular_weight"
)
required_species_columns <- c(
  "Analysis_species_key", "WCVP_taxon_id", "Species_current", "Genus_current",
  "Research_studied_binary", "Chemistry_source_taxa_collapsed",
  "Chemistry_name_original", "Research_article_count",
  "Chem_richness_entities_observed", "Chem_richness_exact_structures",
  "Chem_richness_classes", "Chem_organ_primary_count"
)

check_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(label, " nao contem: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}
check_columns(articles, required_article_columns, "Articles_250")
check_columns(evidence, required_evidence_columns, "Evidence_presence")
check_columns(entities, required_entity_columns, "Chemical_entities")
check_columns(species_universe, required_species_columns, "01_species_analysis_base")

if (nrow(articles) != 250L || uniqueN(articles$Article_ID) != 250L) {
  stop("Esperados 250 Article_ID unicos.", call. = FALSE)
}
if (nrow(evidence) != 6287L || uniqueN(evidence$Evidence_record_ID) != 6287L) {
  stop("Esperados 6.287 Evidence_record_ID unicos.", call. = FALSE)
}
if (nrow(entities) != 1682L || uniqueN(entities$Chemical_entity_ID) != 1682L) {
  stop("Esperadas 1.682 entidades quimicas unicas.", call. = FALSE)
}
if (nrow(species_universe) != 2582L ||
    uniqueN(species_universe$WCVP_taxon_id) != 2582L) {
  stop("Esperadas 2.582 especies WCVP unicas.", call. = FALSE)
}

trim_to_na <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x == "NA"] <- NA_character_
  x
}

is_yes <- function(x) {
  tolower(trimws(as.character(x))) %chin% c("yes", "sim", "1", "true")
}

best_evidence_grade <- function(x) {
  x <- toupper(trim_to_na(x))
  grade_order <- c("A", "B", "C", "D")
  available <- grade_order[grade_order %chin% x]
  if (length(available) == 0L) NA_character_ else available[1L]
}

safe_unique_n <- function(x) uniqueN(x[!is.na(x) & trimws(as.character(x)) != ""])

articles[, Article_ID := trim_to_na(Article_ID)]
evidence[, `:=`(
  Evidence_record_ID = trim_to_na(Evidence_record_ID),
  Article_ID = trim_to_na(Article_ID),
  Taxon_ID = trim_to_na(Taxon_ID),
  Taxon_standardized = trim_to_na(Taxon_standardized),
  Genus = trim_to_na(Genus),
  Species_standardized = normalize_species_name(Species_standardized),
  Organ_primary = trim_to_na(Organ_primary),
  Compound_ID = trim_to_na(Compound_ID),
  Chemical_entity_ID = trim_to_na(Chemical_entity_ID)
)]
entities[, `:=`(
  Chemical_entity_ID = trim_to_na(Chemical_entity_ID),
  Exact_structure_ID = trim_to_na(Exact_structure_ID),
  Chemical_domain = trim_to_na(Chemical_domain),
  Chemical_class_harmonized = trim_to_na(Chemical_class_harmonized),
  Classification_confidence = trim_to_na(Classification_confidence)
)]
species_universe[, `:=`(
  WCVP_taxon_id = sub("\\.0$", "", as.character(WCVP_taxon_id)),
  Species_current = normalize_species_name(Species_current),
  Genus_current = trim_to_na(Genus_current),
  Chemistry_name_original = trim_to_na(Chemistry_name_original)
)]

if (any(!evidence$Article_ID %chin% articles$Article_ID)) {
  stop("Ha registros de evidencia com Article_ID ausente de Articles_250.", call. = FALSE)
}
if (any(!evidence$Chemical_entity_ID %chin% entities$Chemical_entity_ID)) {
  stop("Ha registros de evidencia sem entidade no dicionario quimico.", call. = FALSE)
}

message("[02/07] Aplicando a taxonomia corrente e a hierarquia quimica...")

# Crosswalk reverso: as 177 unidades-fonte do corpus convergem para 173 especies
# correntes. Esta tabela reaproveita as resolucoes taxonomicas congeladas na
# etapa 01, sem consultar novamente servicos externos.
studied_current <- species_universe[Research_studied_binary == 1L]
if (nrow(studied_current) != 173L) {
  stop("Esperadas 173 especies correntes estudadas na etapa 01.", call. = FALSE)
}

source_to_current <- studied_current[, {
  source_names <- trimws(unlist(strsplit(Chemistry_name_original, "\\|")))
  source_names <- source_names[nzchar(source_names)]
  .(Species_standardized = source_names)
}, by = .(
  Analysis_species_key, WCVP_taxon_id, Species_current, Genus_current,
  Chemistry_source_taxa_collapsed
)]

if (nrow(source_to_current) != 177L ||
    uniqueN(source_to_current$Species_standardized) != 177L) {
  stop("O crosswalk deveria conter 177 unidades taxonomicas-fonte unicas.", call. = FALSE)
}

evidence <- merge(
  evidence, source_to_current,
  by = "Species_standardized", all.x = TRUE, sort = FALSE
)
evidence[, Species_ready_source := is_yes(Species_analysis_ready)]

unmapped_species <- evidence[
  Species_ready_source & is.na(WCVP_taxon_id),
  unique(Species_standardized)
]
if (length(unmapped_species) > 0L) {
  stop(
    "Registros prontos para especie sem destino corrente: ",
    paste(unmapped_species, collapse = " | "), call. = FALSE
  )
}

evidence[, `:=`(
  Analysis_taxon_key = fifelse(
    Species_ready_source,
    paste0("WCVP_", WCVP_taxon_id),
    paste0("SOURCE_", Taxon_ID)
  ),
  Taxon_current_analysis = fifelse(
    Species_ready_source, Species_current,
    fifelse(!is.na(Taxon_standardized), Taxon_standardized, Species_standardized)
  ),
  Genus_current_analysis = fifelse(
    Species_ready_source, Genus_current, Genus
  ),
  Species_current_analysis = fifelse(
    Species_ready_source, Species_current, NA_character_
  ),
  Taxon_analysis_level = fifelse(
    Species_ready_source, "Current species", Taxon_type
  )
)]

# Pequena correcao de nomenclatura na camada analitica. O nome original e
# preservado no dicionario; apenas duas categorias equivalentes sao fundidas.
entities[, `:=`(
  Domain_analysis = Chemical_domain,
  Class_analysis = fifelse(
    Chemical_class_harmonized %chin% c("Volatile", "Volatile organic compound"),
    "Volatile organic compounds",
    Chemical_class_harmonized
  ),
  Exact_structure_resolved = as.integer(!is.na(Exact_structure_ID)),
  Class_resolved = as.integer(
    Classification_confidence %chin% c("High", "Medium") &
      Chemical_class_harmonized != "Other or unresolved chemical entity"
  ),
  Metabolite_scope = as.integer(Chemical_domain != "Inorganic/non-metabolite"),
  Sulfated_flag = as.integer(is_yes(Sulfated))
)]

entity_join <- entities[, .(
  Chemical_entity_ID,
  Entity_name_analysis = Entity_canonical_name,
  Domain_analysis,
  Class_original_analysis = Chemical_class_harmonized,
  Class_analysis,
  Class_confidence_analysis = Classification_confidence,
  Class_resolved,
  Metabolite_scope,
  Sulfated_flag,
  Exact_structure_ID_analysis = Exact_structure_ID,
  Exact_structure_resolved,
  Structure_status_analysis = Structure_status,
  Validated_PubChem_CID,
  Validated_InChIKey,
  Molecular_formula,
  Molecular_weight
)]
evidence[entity_join, on = "Chemical_entity_ID", `:=`(
  Entity_name_analysis = i.Entity_name_analysis,
  Domain_analysis = i.Domain_analysis,
  Class_original_analysis = i.Class_original_analysis,
  Class_analysis = i.Class_analysis,
  Class_confidence_analysis = i.Class_confidence_analysis,
  Class_resolved = i.Class_resolved,
  Metabolite_scope = i.Metabolite_scope,
  Sulfated_flag = i.Sulfated_flag,
  Exact_structure_ID_analysis = i.Exact_structure_ID_analysis,
  Exact_structure_resolved = i.Exact_structure_resolved,
  Structure_status_analysis = i.Structure_status_analysis,
  Validated_PubChem_CID_analysis = i.Validated_PubChem_CID,
  Validated_InChIKey_analysis = i.Validated_InChIKey,
  Molecular_formula_analysis = i.Molecular_formula,
  Molecular_weight_analysis = i.Molecular_weight
)]

if (anyNA(evidence$Domain_analysis) || anyNA(evidence$Class_analysis)) {
  stop("Falha ao anexar a classificacao quimica a todos os registros.", call. = FALSE)
}

# A Etapa 01 contou a classe armazenada em cada linha de evidencia. A Etapa 04
# usa a classe canonica, unica por Chemical_entity_ID, proveniente do dicionario
# consolidado. As diferencas sao documentadas, mas nao representam perda de
# artigos, entidades ou estruturas.
evidence[, `:=`(
  Class_evidence_original = trim_to_na(Chemical_class_harmonized),
  Class_entity_canonical = Class_original_analysis
)]
class_reclassification <- evidence[
  fcoalesce(Class_evidence_original, "<NA>") !=
    fcoalesce(Class_entity_canonical, "<NA>"),
  .(
    Evidence_record_count = .N,
    Article_count = uniqueN(Article_ID),
    Source_species_unit_count = safe_unique_n(Species_standardized),
    Current_species_count = safe_unique_n(Species_current_analysis),
    Evidence_record_IDs = collapse_sorted(Evidence_record_ID)
  ),
  by = .(
    Chemical_entity_ID,
    Entity_name_analysis,
    Class_evidence_original,
    Class_entity_canonical,
    Class_confidence_analysis
  )
]
setorder(
  class_reclassification,
  -Evidence_record_count, Chemical_entity_ID, Class_evidence_original
)

# As familias analiticas sao multi-rotulo: um registro pode usar, por exemplo,
# LC-MS e NMR. Por isso as porcentagens por metodo nao devem somar 100%.
evidence[, Analytical_text := tolower(fcoalesce(trim_to_na(Methods), ""))]
evidence[, `:=`(
  Method_GC = as.integer(grepl(
    "(^|[^a-z])gc([^a-z]|$)|gas chromat|headspace|fatty[- ]acid profile",
    Analytical_text, perl = TRUE
  )),
  Method_LC = as.integer(grepl(
    "hplc|uplc|uhplc|(^|[^a-z])lc([^a-z]|$)|liquid chromat",
    Analytical_text, perl = TRUE
  )),
  Method_NMR_isolation = as.integer(grepl(
    "nmr|hmbc|hsqc|cosy|noesy|structure elucidation|isolat",
    Analytical_text, perl = TRUE
  )),
  Method_MS = as.integer(grepl(
    "mass spect|(^|[^a-z])ms([^a-z]|$)|qtof|orbitrap|maldi|esi|apci",
    Analytical_text, perl = TRUE
  )),
  Method_other_assay = as.integer(grepl(
    "spectrophot|colorimet|biochemical|assay|enzyme|titration|gravimet|ion chromat",
    Analytical_text, perl = TRUE
  ))
)]
evidence[, Method_unspecified := as.integer(
  Method_GC + Method_LC + Method_NMR_isolation + Method_MS +
    Method_other_assay == 0L
)]

message("[03/07] Construindo a incidencia deduplicada...")

incidence <- evidence[, .(
  Evidence_record_count = .N,
  Evidence_record_IDs = collapse_sorted(Evidence_record_ID),
  Source_record_IDs = collapse_sorted(Source_record_ID),
  Compound_IDs = collapse_sorted(Compound_ID),
  Best_evidence_grade = best_evidence_grade(Evidence_grade),
  Methods_combined = collapse_sorted(Methods),
  Method_GC = max(Method_GC, na.rm = TRUE),
  Method_LC = max(Method_LC, na.rm = TRUE),
  Method_NMR_isolation = max(Method_NMR_isolation, na.rm = TRUE),
  Method_MS = max(Method_MS, na.rm = TRUE),
  Method_other_assay = max(Method_other_assay, na.rm = TRUE),
  Method_unspecified = max(Method_unspecified, na.rm = TRUE)
), by = .(
  Article_ID, DOI, Analysis_taxon_key, WCVP_taxon_id,
  Taxon_current_analysis, Genus_current_analysis,
  Species_current_analysis, Taxon_analysis_level, Species_ready_source,
  Organ_primary, Chemical_entity_ID, Entity_name_analysis,
  Domain_analysis, Class_original_analysis, Class_analysis, Class_confidence_analysis,
  Class_resolved, Metabolite_scope, Sulfated_flag,
  Exact_structure_ID_analysis, Exact_structure_resolved,
  Structure_status_analysis, Validated_PubChem_CID_analysis,
  Validated_InChIKey_analysis, Molecular_formula_analysis,
  Molecular_weight_analysis
)]
setorder(incidence, Article_ID, Analysis_taxon_key, Organ_primary, Chemical_entity_ID)
incidence[, Incidence_record_ID := sprintf("I%05d", .I)]
setcolorder(incidence, c("Incidence_record_ID", setdiff(names(incidence), "Incidence_record_ID")))

write_csv_semicolon(
  incidence,
  file.path(STAGE04_DATA_DIR, "04_chemical_incidence_article_taxon_organ_entity.csv")
)
write_csv_semicolon(
  source_to_current,
  file.path(STAGE04_DATA_DIR, "04_source_taxon_to_current_species_crosswalk.csv")
)

message("[04/07] Calculando arquitetura por dominio, classe, orgao e taxon...")

total_current_species <- 173L
total_articles <- 250L

summarise_architecture <- function(data, group_columns) {
  data[, .(
    Article_count = uniqueN(Article_ID),
    Current_species_count = safe_unique_n(Species_current_analysis),
    Genus_count = safe_unique_n(Genus_current_analysis),
    Organ_count = safe_unique_n(Organ_primary),
    Incidence_record_count = .N,
    Evidence_record_count = sum(Evidence_record_count),
    Chemical_entity_count = uniqueN(Chemical_entity_ID),
    Exact_structure_count = safe_unique_n(Exact_structure_ID_analysis),
    Resolved_class_count = safe_unique_n(Class_analysis[Class_resolved == 1L])
  ), by = group_columns]
}

domain_summary <- summarise_architecture(incidence, "Domain_analysis")
domain_summary[, `:=`(
  Entity_fraction_total = Chemical_entity_count / uniqueN(incidence$Chemical_entity_ID),
  Species_prevalence_current = Current_species_count / total_current_species,
  Article_prevalence = Article_count / total_articles
)]
setorder(domain_summary, -Chemical_entity_count, Domain_analysis)

# A tabela por classe exclui entidades de baixa confianca. Elas continuam na
# riqueza total e nos dominios, mas nao sao forcadas em uma classe especifica.
class_summary <- summarise_architecture(
  incidence[Class_resolved == 1L],
  c("Domain_analysis", "Class_analysis", "Metabolite_scope")
)
class_summary[, Class_resolved := 1L]
class_summary[, `:=`(
  Entity_fraction_total = Chemical_entity_count / uniqueN(incidence$Chemical_entity_ID),
  Species_prevalence_current = Current_species_count / total_current_species,
  Article_prevalence = Article_count / total_articles,
  Entities_per_article_observed = Chemical_entity_count / Article_count,
  Entities_per_species_observed = fifelse(
    Current_species_count > 0L,
    Chemical_entity_count / Current_species_count,
    NA_real_
  )
)]
setorder(class_summary, -Chemical_entity_count, Class_analysis)

organ_summary <- summarise_architecture(incidence, "Organ_primary")
organ_summary[, `:=`(
  Entity_fraction_total = Chemical_entity_count / uniqueN(incidence$Chemical_entity_ID),
  Species_fraction_current = Current_species_count / total_current_species,
  Article_fraction = Article_count / total_articles,
  Entities_per_article_observed = Chemical_entity_count / Article_count
)]
setorder(organ_summary, -Article_count, Organ_primary)

species_incidence <- incidence[Species_ready_source == TRUE]
species_summary <- species_incidence[, .(
  Article_count = uniqueN(Article_ID),
  Organ_count = safe_unique_n(Organ_primary),
  Incidence_record_count = .N,
  Evidence_record_count = sum(Evidence_record_count),
  Chemical_entity_count = uniqueN(Chemical_entity_ID),
  Exact_structure_count = safe_unique_n(Exact_structure_ID_analysis),
  Chemical_domain_count = safe_unique_n(Domain_analysis),
  Chemical_class_count_source_definition = safe_unique_n(Class_original_analysis),
  Chemical_class_count_resolved_analysis = safe_unique_n(
    Class_analysis[Class_resolved == 1L]
  ),
  Sulfated_entity_count = uniqueN(Chemical_entity_ID[Sulfated_flag == 1L])
), by = .(
  Analysis_species_key = Analysis_taxon_key,
  WCVP_taxon_id,
  Species_current = Species_current_analysis,
  Genus_current = Genus_current_analysis
)]
species_summary[, `:=`(
  Entities_per_article_observed = Chemical_entity_count / Article_count,
  Resolved_classes_per_article_observed =
    Chemical_class_count_resolved_analysis / Article_count
)]
setorder(species_summary, -Article_count, -Chemical_entity_count, Species_current)

genus_summary <- species_incidence[, .(
  Article_count = uniqueN(Article_ID),
  Current_species_count = safe_unique_n(Species_current_analysis),
  Organ_count = safe_unique_n(Organ_primary),
  Incidence_record_count = .N,
  Evidence_record_count = sum(Evidence_record_count),
  Chemical_entity_count = uniqueN(Chemical_entity_ID),
  Exact_structure_count = safe_unique_n(Exact_structure_ID_analysis),
  Chemical_domain_count = safe_unique_n(Domain_analysis),
  Chemical_class_count_source_definition = safe_unique_n(Class_original_analysis),
  Chemical_class_count_resolved_analysis = safe_unique_n(
    Class_analysis[Class_resolved == 1L]
  ),
  Sulfated_entity_count = uniqueN(Chemical_entity_ID[Sulfated_flag == 1L])
), by = .(Genus_current = Genus_current_analysis)]
genus_summary[, `:=`(
  Entities_per_article_observed = Chemical_entity_count / Article_count,
  Resolved_classes_per_article_observed =
    Chemical_class_count_resolved_analysis / Article_count
)]
setorder(genus_summary, -Article_count, -Chemical_entity_count, Genus_current)

# Incidencia especie x orgao x classe. Cada especie contribui no maximo uma
# vez por celula, independentemente do numero de artigos, compostos ou sinonimos.
species_organ_denominator <- unique(
  species_incidence[, .(Analysis_taxon_key, Organ_primary)]
)[, .(Species_studied_in_organ = uniqueN(Analysis_taxon_key)), by = Organ_primary]

organ_class_species <- unique(species_incidence[
  Class_resolved == 1L & Metabolite_scope == 1L,
  .(Analysis_taxon_key, Organ_primary, Domain_analysis, Class_analysis)
])
organ_class_summary <- organ_class_species[, .(
  Species_class_count = uniqueN(Analysis_taxon_key)
), by = .(Organ_primary, Domain_analysis, Class_analysis)]
organ_class_summary <- merge(
  organ_class_summary, species_organ_denominator,
  by = "Organ_primary", all.x = TRUE, sort = FALSE
)

organ_class_effort <- species_incidence[
  Class_resolved == 1L & Metabolite_scope == 1L,
  .(
    Article_count = uniqueN(Article_ID),
    Chemical_entity_count = uniqueN(Chemical_entity_ID),
    Incidence_record_count = .N
  ),
  by = .(Organ_primary, Domain_analysis, Class_analysis)
]
organ_class_summary <- merge(
  organ_class_summary, organ_class_effort,
  by = c("Organ_primary", "Domain_analysis", "Class_analysis"),
  all.x = TRUE, sort = FALSE
)
organ_class_summary[, Species_prevalence_within_organ :=
  Species_class_count / Species_studied_in_organ]
setorder(organ_class_summary, Organ_primary, -Species_class_count, Class_analysis)

genus_denominator <- unique(species_incidence[, .(
  Genus_current = Genus_current_analysis,
  Analysis_taxon_key
)])[, .(Studied_species_in_genus = uniqueN(Analysis_taxon_key)), by = Genus_current]

genus_class_species <- unique(species_incidence[
  Class_resolved == 1L & Metabolite_scope == 1L,
  .(
    Genus_current = Genus_current_analysis,
    Analysis_taxon_key,
    Domain_analysis,
    Class_analysis
  )
])
genus_class_summary <- genus_class_species[, .(
  Species_class_count = uniqueN(Analysis_taxon_key)
), by = .(Genus_current, Domain_analysis, Class_analysis)]
genus_class_summary <- merge(
  genus_class_summary, genus_denominator,
  by = "Genus_current", all.x = TRUE, sort = FALSE
)
genus_class_effort <- species_incidence[
  Class_resolved == 1L & Metabolite_scope == 1L,
  .(
    Article_count = uniqueN(Article_ID),
    Chemical_entity_count = uniqueN(Chemical_entity_ID),
    Incidence_record_count = .N
  ),
  by = .(
    Genus_current = Genus_current_analysis,
    Domain_analysis,
    Class_analysis
  )
]
genus_class_summary <- merge(
  genus_class_summary, genus_class_effort,
  by = c("Genus_current", "Domain_analysis", "Class_analysis"),
  all.x = TRUE, sort = FALSE
)
genus_class_summary[, Species_prevalence_within_genus :=
  Species_class_count / Studied_species_in_genus]
setorder(genus_class_summary, Genus_current, -Species_class_count, Class_analysis)

message("[05/07] Auditando identificacao, metodos e compostos sulfatados...")

structure_resolution <- entities[, .(
  Chemical_entity_count = .N
), by = .(Category = Structure_status)]
structure_resolution[, Dimension := "Structure resolution"]
structure_evidence <- incidence[, .(
  Incidence_record_count = .N,
  Evidence_record_count = sum(Evidence_record_count),
  Article_count = uniqueN(Article_ID),
  Current_species_count = safe_unique_n(Species_current_analysis)
), by = .(Category = Structure_status_analysis)]
structure_resolution <- merge(
  structure_resolution, structure_evidence,
  by = "Category", all = TRUE, sort = FALSE
)

class_resolution <- entities[, .(
  Chemical_entity_count = .N
), by = .(Category = Classification_confidence)]
class_resolution[, Dimension := "Class assignment confidence"]
class_evidence <- incidence[, .(
  Incidence_record_count = .N,
  Evidence_record_count = sum(Evidence_record_count),
  Article_count = uniqueN(Article_ID),
  Current_species_count = safe_unique_n(Species_current_analysis)
), by = .(Category = Class_confidence_analysis)]
class_resolution <- merge(
  class_resolution, class_evidence,
  by = "Category", all = TRUE, sort = FALSE
)
identification_resolution <- rbindlist(
  list(structure_resolution, class_resolution), use.names = TRUE, fill = TRUE
)
identification_resolution[, Entity_fraction :=
  Chemical_entity_count / uniqueN(entities$Chemical_entity_ID)]
setcolorder(
  identification_resolution,
  c("Dimension", "Category", setdiff(names(identification_resolution), c("Dimension", "Category")))
)

method_columns <- c(
  "Method_GC", "Method_LC", "Method_NMR_isolation", "Method_MS",
  "Method_other_assay", "Method_unspecified"
)
method_labels <- c(
  Method_GC = "GC-based",
  Method_LC = "LC-based",
  Method_NMR_isolation = "Isolation/NMR",
  Method_MS = "Mass spectrometry",
  Method_other_assay = "Other biochemical/targeted assay",
  Method_unspecified = "Unspecified/other"
)
method_long <- melt(
  incidence[Class_resolved == 1L & Metabolite_scope == 1L],
  id.vars = c(
    "Article_ID", "Species_current_analysis", "Chemical_entity_ID",
    "Domain_analysis", "Class_analysis"
  ),
  measure.vars = method_columns,
  variable.name = "Method_code", value.name = "Method_present"
)[Method_present == 1L]
method_long[, Analytical_family := unname(method_labels[as.character(Method_code)])]

method_class_summary <- method_long[, .(
  Article_count = uniqueN(Article_ID),
  Current_species_count = safe_unique_n(Species_current_analysis),
  Chemical_entity_count = uniqueN(Chemical_entity_ID)
), by = .(Domain_analysis, Class_analysis, Analytical_family)]
class_article_denominator <- incidence[
  Class_resolved == 1L & Metabolite_scope == 1L,
  .(Class_article_count = uniqueN(Article_ID)),
  by = .(Domain_analysis, Class_analysis)
]
method_class_summary <- merge(
  method_class_summary, class_article_denominator,
  by = c("Domain_analysis", "Class_analysis"), all.x = TRUE, sort = FALSE
)
method_class_summary[, Percent_class_articles_using_method :=
  Article_count / Class_article_count]
setorder(method_class_summary, Class_analysis, -Article_count, Analytical_family)

sulfated_entities <- entities[Sulfated_flag == 1L]
sulfated_breadth <- incidence[Sulfated_flag == 1L, .(
  Article_count = uniqueN(Article_ID),
  Current_species_count = safe_unique_n(Species_current_analysis),
  Genus_count = safe_unique_n(Genus_current_analysis),
  Organ_count = safe_unique_n(Organ_primary),
  Incidence_record_count = .N,
  Evidence_record_count = sum(Evidence_record_count)
), by = Chemical_entity_ID]
sulfated_table <- merge(
  sulfated_entities[, .(
    Chemical_entity_ID, Entity_canonical_name, Chemical_domain,
    Class_original = Chemical_class_harmonized,
    Class_analysis, Classification_confidence,
    Exact_structure_ID, Structure_status, Validated_PubChem_CID,
    Validated_InChIKey, Molecular_formula, Molecular_weight
  )],
  sulfated_breadth, by = "Chemical_entity_ID", all.x = TRUE, sort = FALSE
)
setorder(sulfated_table, -Article_count, -Current_species_count, Entity_canonical_name)

metric_definitions <- data.table(
  Metric = c(
    "Evidence record", "Analytical incidence record", "Observed chemical richness",
    "Exact-structure richness", "Resolved-class richness", "Species prevalence",
    "Article prevalence", "Method association"
  ),
  Unit = c(
    "Evidence_record_ID",
    "Article_ID x analysis taxon x Organ_primary x Chemical_entity_ID",
    "Unique Chemical_entity_ID",
    "Unique non-missing Exact_structure_ID",
    "Unique Class_analysis among High/Medium-confidence assignments",
    "Unique current species reporting a domain/class divided by 173",
    "Unique articles reporting a domain/class divided by 250",
    "Unique articles with a method flag and a class; method flags are multi-label"
  ),
  Interpretation = c(
    "Rastreabilidade da matriz; repeated synonyms/entities are retained.",
    "Incidence layer used to prevent repeated source rows from inflating counts.",
    "Primary diversity measure; includes unresolved annotations as distinct observed entities.",
    "Sensitivity subset only; it does not represent all observed chemistry.",
    "Used for class-composition analyses; low-confidence unresolved entities are not forced into classes.",
    "Taxonomic breadth, not chemical richness or abundance.",
    "Bibliographic breadth, not biological prevalence.",
    "Analytical context; percentages across methods can sum to more than 100%."
  )
)

dataset_scope <- data.table(
  Indicator = c(
    "Informative articles",
    "Original evidence records",
    "Deduplicated analytical incidence records",
    "Observed chemical entities",
    "Entities with resolved class",
    "Entities with exact structure",
    "Source species units",
    "Current species after WCVP reconciliation",
    "Current genera represented",
    "Primary organ categories",
    "Resolved analytical chemical classes after volatile-class fusion",
    "Sulfated chemical entities"
  ),
  Value = c(
    uniqueN(articles$Article_ID),
    uniqueN(evidence$Evidence_record_ID),
    nrow(incidence),
    uniqueN(entities$Chemical_entity_ID),
    sum(entities$Class_resolved == 1L),
    sum(entities$Exact_structure_resolved == 1L),
    uniqueN(source_to_current$Species_standardized),
    uniqueN(species_summary$WCVP_taxon_id),
    uniqueN(species_summary$Genus_current),
    uniqueN(incidence$Organ_primary),
    uniqueN(entities[Class_resolved == 1L, Class_analysis]),
    sum(entities$Sulfated_flag == 1L)
  ),
  Unit = c(
    "Article_ID", "Evidence_record_ID", "Incidence_record_ID",
    "Chemical_entity_ID", "Chemical_entity_ID", "Exact_structure_ID",
    "Source species name", "WCVP_taxon_id", "Genus_current", "Organ_primary",
    "Class_analysis", "Chemical_entity_ID"
  )
)
dataset_scope[, Fraction_of_entities := fifelse(
  Indicator %chin% c("Entities with resolved class", "Entities with exact structure"),
  Value / 1682, NA_real_
)]

# Reconciliacao independente com os totais por especie congelados na etapa 01.
species_reference <- studied_current[, .(
  WCVP_taxon_id,
  Reference_articles = to_numeric_safe(Research_article_count),
  Reference_entities = to_numeric_safe(Chem_richness_entities_observed),
  Reference_exact = to_numeric_safe(Chem_richness_exact_structures),
  Reference_classes = to_numeric_safe(Chem_richness_classes),
  Reference_organs = to_numeric_safe(Chem_organ_primary_count)
)]
species_reconciliation <- merge(
  species_summary, species_reference,
  by = "WCVP_taxon_id", all = TRUE, sort = FALSE
)
species_reconciliation[, `:=`(
  Difference_articles = Article_count - Reference_articles,
  Difference_entities = Chemical_entity_count - Reference_entities,
  Difference_exact = Exact_structure_count - Reference_exact,
  Difference_classes = Chemical_class_count_source_definition - Reference_classes,
  Difference_organs = Organ_count - Reference_organs
)]
reconciliation_core_failures <- species_reconciliation[
  abs(Difference_articles) > 0 |
    abs(Difference_entities) > 0 |
    abs(Difference_exact) > 0 |
    abs(Difference_organs) > 0
]
reconciliation_class_changes <- species_reconciliation[
  abs(Difference_classes) > 0
]

qa_row <- function(indicator, observed, expected = NA_real_, interpretation = "") {
  status <- if (is.na(expected)) {
    "INFO"
  } else if (isTRUE(all.equal(as.numeric(observed), as.numeric(expected)))) {
    "PASS"
  } else {
    "FAIL"
  }
  data.table(
    Indicator = indicator, Observed = observed, Expected = expected,
    Status = status, Interpretation = interpretation
  )
}

qa <- rbindlist(list(
  qa_row("Articles", uniqueN(articles$Article_ID), 250, "One row per informative article"),
  qa_row("Evidence records", uniqueN(evidence$Evidence_record_ID), 6287, "Original records preserved"),
  qa_row("Chemical entities", uniqueN(entities$Chemical_entity_ID), 1682, "Observed entities preserved"),
  qa_row("Class-resolved entities", sum(entities$Class_resolved == 1L), 1287, "High or medium confidence"),
  qa_row("Exact structures", sum(entities$Exact_structure_resolved == 1L), 791, "Validated provisional exact structures"),
  qa_row("Source species units", uniqueN(source_to_current$Species_standardized), 177, "Legacy/source names"),
  qa_row("Current studied species", uniqueN(species_summary$WCVP_taxon_id), 173, "After WCVP reconciliation"),
  qa_row("Species-ready evidence mapped", sum(evidence$Species_ready_source & !is.na(evidence$WCVP_taxon_id)), 6088, "All species-ready records"),
  qa_row("Unmapped species-ready evidence", sum(evidence$Species_ready_source & is.na(evidence$WCVP_taxon_id)), 0, "Must remain zero"),
  qa_row("Deduplicated incidence records", nrow(incidence), 5207, "Article x taxon x organ x entity"),
  qa_row("Evidence rows collapsed into incidence", nrow(evidence) - nrow(incidence), 1080, "Not deleted; traceable via Evidence_record_IDs"),
  qa_row("Current genera represented", uniqueN(species_summary$Genus_current), 75, "Species-resolved layer"),
  qa_row("Primary organ categories", uniqueN(incidence$Organ_primary), 10, "Harmonized organs"),
  qa_row("Resolved analytical classes", uniqueN(entities[Class_resolved == 1L, Class_analysis]), 44, "After volatile-class fusion"),
  qa_row("Sulfated entities", sum(entities$Sulfated_flag == 1L), 26, "Observed entity unit"),
  qa_row("Species reconciliation core failures", nrow(reconciliation_core_failures), 0, "Articles, entities, exact structures and organs"),
  qa_row("Species class-count changes after canonical reclassification", nrow(reconciliation_class_changes), 24, "Expected change from evidence-row class to entity-dictionary class"),
  qa_row("Evidence records with canonical class change", sum(class_reclassification$Evidence_record_count), 192, "Documented in Table_QA_04_class_reclassification"),
  qa_row("Entities with canonical class change", uniqueN(class_reclassification$Chemical_entity_ID), 62, "Canonical entity dictionary is authoritative"),
  qa_row("Source species units affected by canonical class change", safe_unique_n(evidence[
    fcoalesce(Class_evidence_original, "<NA>") != fcoalesce(Class_entity_canonical, "<NA>"),
    Species_standardized
  ]), 56, "No species, article or entity is removed"),
  qa_row("Chemical domains", uniqueN(entities$Domain_analysis), 11, "Complete entity dictionary")
), use.names = TRUE)

if (any(qa$Status == "FAIL")) {
  write_csv_semicolon(qa, file.path(STAGE04_TABLES_DIR, "Table_QA_04_chemical_architecture.csv"))
  write_csv_semicolon(
    species_reconciliation,
    file.path(STAGE04_TABLES_DIR, "Table_QA_04_species_reconciliation.csv")
  )
  write_csv_semicolon(
    class_reclassification,
    file.path(STAGE04_TABLES_DIR, "Table_QA_04_class_reclassification.csv")
  )
  stop("A auditoria da etapa 04 encontrou falhas. Consulte Table_QA_04.", call. = FALSE)
}

write_csv_semicolon(dataset_scope, file.path(STAGE04_TABLES_DIR, "Table_04_01_dataset_scope.csv"))
write_csv_semicolon(domain_summary, file.path(STAGE04_TABLES_DIR, "Table_04_02_domain_architecture.csv"))
write_csv_semicolon(class_summary, file.path(STAGE04_TABLES_DIR, "Table_04_03_class_architecture.csv"))
write_csv_semicolon(organ_summary, file.path(STAGE04_TABLES_DIR, "Table_04_04_organ_architecture.csv"))
write_csv_semicolon(organ_class_summary, file.path(STAGE04_TABLES_DIR, "Table_04_05_organ_class_incidence.csv"))
write_csv_semicolon(genus_summary, file.path(STAGE04_TABLES_DIR, "Table_04_06_genus_architecture.csv"))
write_csv_semicolon(genus_class_summary, file.path(STAGE04_TABLES_DIR, "Table_04_07_genus_class_incidence.csv"))
write_csv_semicolon(species_summary, file.path(STAGE04_TABLES_DIR, "Table_04_08_species_architecture.csv"))
write_csv_semicolon(identification_resolution, file.path(STAGE04_TABLES_DIR, "Table_04_09_identification_resolution.csv"))
write_csv_semicolon(sulfated_table, file.path(STAGE04_TABLES_DIR, "Table_04_10_sulfated_entities.csv"))
write_csv_semicolon(method_class_summary, file.path(STAGE04_TABLES_DIR, "Table_04_11_method_class_association.csv"))
write_csv_semicolon(metric_definitions, file.path(STAGE04_TABLES_DIR, "Table_04_12_metric_definitions.csv"))
write_csv_semicolon(qa, file.path(STAGE04_TABLES_DIR, "Table_QA_04_chemical_architecture.csv"))
write_csv_semicolon(
  species_reconciliation,
  file.path(STAGE04_TABLES_DIR, "Table_QA_04_species_reconciliation.csv")
)
write_csv_semicolon(
  class_reclassification,
  file.path(STAGE04_TABLES_DIR, "Table_QA_04_class_reclassification.csv")
)

message("[06/07] Produzindo figuras candidatas para o artigo...")

domain_palette <- c(
  "Phenylpropanoids and polyketides" = "#3B7A57",
  "Terpenoids and steroids" = "#D17C36",
  "Lipids and lipid-like molecules" = "#44799C",
  "Other specialized metabolites" = "#8B5E83",
  "Carbohydrates" = "#B69B35",
  "Alkaloids and nitrogen compounds" = "#A6483D",
  "Amino acids and peptides" = "#6E8B3D",
  "Other primary metabolites" = "#5E7D7E",
  "Vitamins and cofactors" = "#C08A9A",
  "Inorganic/non-metabolite" = "#8A8A8A",
  "Other/unresolved" = "#B7B7B7"
)

theme_article <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 4),
      plot.subtitle = element_text(color = "#555555", size = base_size + 1),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(12, 24, 12, 12)
    )
}

save_figure <- function(plot, filename, width, height) {
  ggsave(
    file.path(STAGE04_FIGURES_DIR, paste0(filename, ".pdf")),
    plot = plot, width = width, height = height, units = "in", bg = "white"
  )
  ggsave(
    file.path(STAGE04_FIGURES_DIR, paste0(filename, ".png")),
    plot = plot, width = width, height = height, units = "in",
    dpi = 320, bg = "white"
  )
}

# Figure 04.01 — riqueza de entidades e amplitude taxonomica por dominio.
domain_plot_data <- copy(domain_summary)
domain_plot_data[, Domain_analysis := factor(
  Domain_analysis,
  levels = domain_plot_data[order(Chemical_entity_count), Domain_analysis]
)]
domain_plot_data[, Annotation_only_count :=
  Chemical_entity_count - Exact_structure_count]
domain_resolution_long <- melt(
  domain_plot_data,
  id.vars = "Domain_analysis",
  measure.vars = c("Exact_structure_count", "Annotation_only_count"),
  variable.name = "Resolution", value.name = "Entity_count"
)
domain_resolution_long[, Resolution := factor(
  Resolution,
  levels = c("Annotation_only_count", "Exact_structure_count"),
  labels = c("Observed annotation only", "Provisional exact structure")
)]

p_domain_entities <- ggplot(
  domain_resolution_long,
  aes(x = Domain_analysis, y = Entity_count, fill = Resolution)
) +
  geom_col(width = 0.72) +
  geom_text(
    data = domain_plot_data,
    aes(x = Domain_analysis, y = Chemical_entity_count, label = Chemical_entity_count),
    inherit.aes = FALSE, hjust = -0.12, size = 3.3
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c(
    "Observed annotation only" = "#C9D2D3",
    "Provisional exact structure" = "#0F6B66"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(
    title = "Chemical diversity and taxonomic breadth tell different stories",
    subtitle = "Entity richness includes unresolved annotations; structure-resolved richness is a subset",
    x = NULL, y = "Observed chemical entities", fill = "Identification"
  ) +
  theme_article()

p_domain_species <- ggplot(
  domain_plot_data,
  aes(x = Domain_analysis, y = Species_prevalence_current)
) +
  geom_col(width = 0.72, fill = "#D17C36") +
  geom_text(
    aes(label = paste0(Current_species_count, "/", total_current_species)),
    hjust = -0.12, size = 3.3
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(x = NULL, y = "Current species reporting the domain") +
  theme_article() +
  theme(
    plot.title = element_blank(), plot.subtitle = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )

figure_04_01 <- p_domain_entities + p_domain_species +
  plot_layout(widths = c(1.45, 1), guides = "collect") &
  theme(legend.position = "bottom")
save_figure(figure_04_01, "Figure_04_01_domain_architecture", 13, 7.8)

# Figure 04.02 — orgaos estudados, esforco e riqueza observada.
organ_plot_data <- copy(organ_summary)
organ_plot_data[, Organ_primary := factor(
  Organ_primary,
  levels = organ_plot_data[order(Article_count), Organ_primary]
)]

p_organ_articles <- ggplot(organ_plot_data, aes(Organ_primary, Article_count)) +
  geom_col(fill = "#44799C", width = 0.72) +
  geom_text(aes(label = Article_count), hjust = -0.15, size = 3.4) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Fruit and seed studies dominate the documented palm chemistry",
    subtitle = "Article effort and observed entity richness are shown separately",
    x = NULL, y = "Informative articles"
  ) +
  theme_article()

p_organ_entities <- ggplot(organ_plot_data, aes(Organ_primary, Chemical_entity_count)) +
  geom_col(fill = "#0F6B66", width = 0.72) +
  geom_text(
    aes(label = paste0(Chemical_entity_count, " entities; ", Current_species_count, " spp.")),
    hjust = -0.08, size = 3.2
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(x = NULL, y = "Observed chemical entities") +
  theme_article() +
  theme(
    plot.title = element_blank(), plot.subtitle = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )

figure_04_02 <- p_organ_articles + p_organ_entities +
  plot_layout(widths = c(1, 1.25))
save_figure(figure_04_02, "Figure_04_02_organ_sampling_and_richness", 13, 7.5)

# Figure 04.03 — prevalencia de classes dentro de cada orgao.
eligible_organs <- species_organ_denominator[
  Species_studied_in_organ >= 5L &
    !Organ_primary %chin% c("Mixed organs", "Other/unclear"),
  Organ_primary
]
class_breadth <- organ_class_species[, .(
  Total_species = uniqueN(Analysis_taxon_key)
), by = .(Domain_analysis, Class_analysis)]
top_classes <- head(
  class_breadth[order(-Total_species, Class_analysis)], 18L
)$Class_analysis

heatmap_organ <- organ_class_summary[
  Organ_primary %chin% eligible_organs & Class_analysis %chin% top_classes
]
organ_levels <- species_organ_denominator[
  Organ_primary %chin% eligible_organs
][order(Species_studied_in_organ), Organ_primary]
class_levels <- class_breadth[
  Class_analysis %chin% top_classes
][order(-Total_species), Class_analysis]
heatmap_organ[, `:=`(
  Organ_primary = factor(Organ_primary, levels = organ_levels),
  Class_analysis = factor(Class_analysis, levels = rev(class_levels))
)]

figure_04_03 <- ggplot(
  heatmap_organ,
  aes(x = Class_analysis, y = Organ_primary)
) +
  geom_point(
    aes(size = Species_class_count, color = Species_prevalence_within_organ),
    alpha = 0.9
  ) +
  scale_color_gradient(
    low = "#DCE7E5", high = "#0F6B66",
    labels = percent_format(accuracy = 1)
  ) +
  scale_size_area(max_size = 12, breaks = pretty_breaks(n = 4)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Organ-specific chemical profiles after incidence standardization",
    subtitle = "Each current species contributes at most once to an organ-by-class cell",
    x = NULL, y = NULL,
    color = "Species prevalence\nwithin organ",
    size = "Species"
  ) +
  theme_article(10.5) +
  theme(
    axis.text.x = element_text(angle = 52, hjust = 1, vjust = 1),
    panel.grid.major = element_line(color = "#E6E6E6"),
    legend.position = "right"
  )
save_figure(figure_04_03, "Figure_04_03_organ_class_incidence", 14, 8.5)

# Figure 04.04 — classes ricas versus classes amplamente distribuidas.
class_scatter <- class_summary[
  Class_resolved == 1L & Metabolite_scope == 1L &
    Current_species_count > 0L & Chemical_entity_count > 0L
]
label_classes <- unique(c(
  head(class_scatter[order(-Chemical_entity_count)], 8L)$Class_analysis,
  head(class_scatter[order(-Current_species_count)], 8L)$Class_analysis
))

figure_04_04 <- ggplot(
  class_scatter,
  aes(
    x = Current_species_count,
    y = Chemical_entity_count,
    size = Article_count,
    color = Domain_analysis
  )
) +
  geom_point(alpha = 0.78) +
  geom_text(
    data = class_scatter[Class_analysis %chin% label_classes],
    aes(label = Class_analysis),
    size = 3.1, hjust = -0.08, vjust = -0.25,
    check_overlap = TRUE, show.legend = FALSE
  ) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 173)) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 200)) +
  scale_size_area(max_size = 13, breaks = pretty_breaks(n = 4)) +
  scale_color_manual(values = domain_palette) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Chemical classes differ in richness and taxonomic breadth",
    subtitle = "Descriptive values; point size represents the number of informative articles",
    x = "Current species reporting the class (log scale)",
    y = "Observed chemical entities (log scale)",
    size = "Articles", color = "Chemical domain"
  ) +
  theme_article() +
  theme(legend.position = "right")
save_figure(figure_04_04, "Figure_04_04_class_breadth_vs_richness", 12.5, 8.5)

# Figure 04.05 — dependencia da classe em relacao ao metodo analitico.
top_method_classes <- head(
  class_summary[
    Class_resolved == 1L & Metabolite_scope == 1L
  ][order(-Article_count, Class_analysis)],
  20L
)$Class_analysis
method_plot <- method_class_summary[Class_analysis %chin% top_method_classes]
method_class_order <- class_summary[
  Class_analysis %chin% top_method_classes
][order(Article_count), Class_analysis]
method_plot[, `:=`(
  Class_analysis = factor(Class_analysis, levels = method_class_order),
  Analytical_family = factor(
    Analytical_family,
    levels = c(
      "GC-based", "LC-based", "Mass spectrometry", "Isolation/NMR",
      "Other biochemical/targeted assay", "Unspecified/other"
    )
  )
)]

figure_04_05 <- ggplot(
  method_plot,
  aes(x = Analytical_family, y = Class_analysis,
      fill = Percent_class_articles_using_method)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Article_count), size = 3) +
  scale_fill_gradient(
    low = "#F2F4F4", high = "#A6483D",
    labels = percent_format(accuracy = 1), limits = c(0, 1)
  ) +
  labs(
    title = "Analytical methods filter the chemistry that becomes visible",
    subtitle = "Numbers are articles; method flags are multi-label and percentages need not sum to 100%",
    x = NULL, y = NULL,
    fill = "Class articles\nusing method"
  ) +
  theme_article(10.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank()
  )
save_figure(figure_04_05, "Figure_04_05_method_class_association", 12.5, 9)

# Figure 04.06 — generos mais estudados e amplitude de classes.
top_genera <- head(
  genus_summary[order(-Article_count, -Current_species_count)], 20L
)$Genus_current
top_genus_classes <- head(
  class_breadth[order(-Total_species, Class_analysis)], 15L
)$Class_analysis
genus_heatmap <- genus_class_summary[
  Genus_current %chin% top_genera & Class_analysis %chin% top_genus_classes
]
genus_levels <- genus_summary[Genus_current %chin% top_genera][
  order(Article_count), Genus_current
]
genus_class_levels <- class_breadth[Class_analysis %chin% top_genus_classes][
  order(-Total_species), Class_analysis
]
genus_heatmap[, `:=`(
  Genus_current = factor(Genus_current, levels = genus_levels),
  Class_analysis = factor(Class_analysis, levels = rev(genus_class_levels))
)]

figure_04_06 <- ggplot(
  genus_heatmap,
  aes(x = Class_analysis, y = Genus_current)
) +
  geom_point(
    aes(size = Species_class_count, color = Species_prevalence_within_genus),
    alpha = 0.9
  ) +
  scale_color_gradient(
    low = "#E6EDF2", high = "#44799C",
    labels = percent_format(accuracy = 1)
  ) +
  scale_size_area(max_size = 11, breaks = pretty_breaks(n = 4)) +
  labs(
    title = "Chemical-class coverage is uneven among the most studied genera",
    subtitle = "Species incidence, not raw occurrence counts",
    x = NULL, y = NULL,
    color = "Studied species\nin genus",
    size = "Species"
  ) +
  theme_article(10.5) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1),
    panel.grid.major = element_line(color = "#E8E8E8")
  )
save_figure(figure_04_06, "Figure_04_06_genus_class_incidence", 13.5, 9)

message("[07/07] Salvando checkpoint e log...")

results_object <- list(
  version = "1.0.2",
  generated_at = Sys.time(),
  source_to_current = source_to_current,
  incidence = incidence,
  dataset_scope = dataset_scope,
  domain_summary = domain_summary,
  class_summary = class_summary,
  organ_summary = organ_summary,
  organ_class_summary = organ_class_summary,
  genus_summary = genus_summary,
  genus_class_summary = genus_class_summary,
  species_summary = species_summary,
  identification_resolution = identification_resolution,
  sulfated_entities = sulfated_table,
  method_class_summary = method_class_summary,
  class_reclassification = class_reclassification,
  qa = qa
)
saveRDS(
  results_object,
  file.path(STAGE04_DATA_DIR, "04_chemical_architecture_v1_0_2.rds"),
  compress = "xz"
)

ended_at <- Sys.time()
execution_log <- c(
  "ARTIGO ARECACEAE — ETAPA 04",
  "Script_version: 1.0.2",
  paste0("Inicio: ", format(started_at, "%Y-%m-%d %H:%M:%S")),
  paste0("Fim: ", format(ended_at, "%Y-%m-%d %H:%M:%S")),
  paste0("Duracao_min: ", round(as.numeric(difftime(ended_at, started_at, units = "mins")), 2)),
  paste0("Articles: ", uniqueN(articles$Article_ID)),
  paste0("Evidence_records: ", uniqueN(evidence$Evidence_record_ID)),
  paste0("Incidence_records: ", nrow(incidence)),
  paste0("Chemical_entities: ", uniqueN(entities$Chemical_entity_ID)),
  paste0("Class_resolved_entities: ", sum(entities$Class_resolved == 1L)),
  paste0("Exact_structures: ", sum(entities$Exact_structure_resolved == 1L)),
  paste0("Current_species: ", uniqueN(species_summary$WCVP_taxon_id)),
  paste0("Sulfated_entities: ", sum(entities$Sulfated_flag == 1L)),
  paste0("QA_failures: ", sum(qa$Status == "FAIL"))
)
writeLines(
  execution_log,
  file.path(STAGE04_LOGS_DIR, "04_execution_log.txt"),
  useBytes = TRUE
)
write_session_info(file.path(STAGE04_LOGS_DIR, "04_sessionInfo.txt"))

message("ETAPA 04 CONCLUIDA")
message("Resultados: ", STAGE04_DIR)
