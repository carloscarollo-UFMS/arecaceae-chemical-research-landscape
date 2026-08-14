# =============================================================================
# ARTIGO ARECACEAE — 02. GBIF + GHS-POP POR ESPECIE
# Versao: 1.0.2
#
# Unidade de entrada GBIF: registro de ocorrencia.
# Unidade espacial intermediaria: especie WCVP x coordenada unica.
# Unidade de exposicao humana: especie WCVP x celula GHS-POP de 1 km.
# Unidade de saida: um WCVP_taxon_id por linha.
#
# Principios metodologicos:
#   - registros explicitamente ausentes sao excluidos;
#   - somente bases de registro admissiveis sao mantidas;
#   - coordenadas invalidas e incerteza conhecida >10 km sao excluidas;
#   - duplicacoes da mesma especie na mesma coordenada/celula nao multiplicam
#     artificialmente a populacao;
#   - os pontos sao restringidos a areas WGSRPD L3 em que a especie e nativa
#     segundo POWO/WCVP;
#   - soma de populacao e mantida apenas como descritor. A metrica principal
#     para os modelos sera a media de log1p(populacao) nas celulas unicas.
#
# Checkpoints permitem retomar o processamento sem reler 889 MB.
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
source(file.path(SCRIPT_DIR, "00_configurar_projeto_Arecaceae_v1_0_2.R"))

# ---- Parametros auditaveis ---------------------------------------------------

FORCE_REBUILD <- FALSE
RUN_COORDINATECLEANER <- TRUE
MAX_COORDINATE_UNCERTAINTY_M <- 10000

ALLOWED_BASIS_OF_RECORD <- c(
  "HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION",
  "MACHINE_OBSERVATION", "MATERIAL_SAMPLE", "OCCURRENCE"
)

EXPLICIT_NON_NATIVE_PATTERN <- paste(
  c("INTRODUCED", "INVASIVE", "NATURALISED", "NATURALIZED", "MANAGED", "CULTIVATED"),
  collapse = "|"
)

EXCLUDED_COORDINATE_ISSUE_PATTERN <- paste(
  c("ZERO_COORDINATE", "COORDINATE_INVALID", "PRESUMED_SWAPPED_COORDINATE"),
  collapse = "|"
)

required <- c("data.table", "sf", "terra", "CoordinateCleaner")
require_packages(required)
assert_files_exist(
  c(
    GBIF_FILE, GHS_POP_FILE, POWO_NATIVE_FILE, SPECIES_BASE_RDS,
    WCVP_SAFE_NAME_DICTIONARY_FILE
  ),
  c(
    "GBIF", "GHS-POP", "POWO", "tabela por especie da etapa 01",
    "dicionario seguro de nomes WCVP"
  )
)

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(terra)
  library(CoordinateCleaner)
})

started_at <- Sys.time()

checkpoint_1 <- file.path(CHECKPOINT_DIR, "02_01_gbif_filtered_matched_v1_0_2.rds")
checkpoint_2 <- file.path(CHECKPOINT_DIR, "02_02_gbif_unique_coordinates_v1_0_2.rds")
checkpoint_3 <- file.path(CHECKPOINT_DIR, "02_03_gbif_native_coordinates_v1_0_2.rds")
checkpoint_4 <- file.path(CHECKPOINT_DIR, "02_04_species_ghs_cells_v1_0_2.rds")
checkpoint_manifest <- file.path(CHECKPOINT_DIR, "02_checkpoint_manifest_v1_0_2.rds")

input_fingerprint <- list(
  files = data.table(
    file = normalizePath(
      c(
        GBIF_FILE, GHS_POP_FILE, POWO_NATIVE_FILE, SPECIES_BASE_RDS,
        WCVP_SAFE_NAME_DICTIONARY_FILE
      ),
      winslash = "/"
    ),
    size = as.numeric(file.info(c(
      GBIF_FILE, GHS_POP_FILE, POWO_NATIVE_FILE, SPECIES_BASE_RDS,
      WCVP_SAFE_NAME_DICTIONARY_FILE
    ))$size),
    mtime = as.character(file.info(c(
      GBIF_FILE, GHS_POP_FILE, POWO_NATIVE_FILE, SPECIES_BASE_RDS,
      WCVP_SAFE_NAME_DICTIONARY_FILE
    ))$mtime)
  ),
  parameters = list(
    PIPELINE_VERSION = "1.0.2",
    TAXONOMY_RULESET = "WCVP_current_plus_safe_synonyms_plus_2_manual_resolutions",
    MAX_COORDINATE_UNCERTAINTY_M = MAX_COORDINATE_UNCERTAINTY_M,
    RUN_COORDINATECLEANER = RUN_COORDINATECLEANER,
    ALLOWED_BASIS_OF_RECORD = ALLOWED_BASIS_OF_RECORD,
    EXPLICIT_NON_NATIVE_PATTERN = EXPLICIT_NON_NATIVE_PATTERN,
    EXCLUDED_COORDINATE_ISSUE_PATTERN = EXCLUDED_COORDINATE_ISSUE_PATTERN
  )
)

if (file.exists(checkpoint_manifest) && !FORCE_REBUILD) {
  old_fingerprint <- readRDS(checkpoint_manifest)
  if (!identical(old_fingerprint, input_fingerprint)) {
    warning("Os arquivos de entrada mudaram. Os checkpoints serao reconstruidos.")
    FORCE_REBUILD <- TRUE
  }
}
if (FORCE_REBUILD) {
  # Remove somente checkpoints gerados por esta etapa. Assim, uma interrupcao
  # durante a reconstrucao nunca permite retomar dados com regras antigas.
  unlink(c(checkpoint_1, checkpoint_2, checkpoint_3, checkpoint_4), force = TRUE)
}
saveRDS(input_fingerprint, checkpoint_manifest)

species_base <- as.data.table(readRDS(SPECIES_BASE_RDS))
species_base[, WCVP_taxon_id := as.character(WCVP_taxon_id)]

# ---- Dicionario conservador GBIF -> WCVP ------------------------------------

normalize_species_key <- function(x) {
  tolower(normalize_species_name(x))
}

# Resolucao manual apenas dos dois binomios sem autoria que afetam diretamente
# o corpus quimico. O restante continua submetido as regras conservadoras do
# dicionario WCVP; nomes com mais de um destino sao excluidos.
manual_alias <- data.table(
  GBIF_name_normalized = c("astrocaryum aculeatum", "phoenix pusilla"),
  WCVP_taxon_id = c("17583", "152708"),
  alias_source = "manual_resolution_documented",
  alias_priority = 0L
)
manual_alias <- merge(
  manual_alias,
  species_base[, .(
    WCVP_taxon_id, Analysis_species_key, Species_current, Genus_current
  )],
  by = "WCVP_taxon_id",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(manual_alias$Analysis_species_key)) {
  stop("Os IDs das resolucoes manuais nao existem na tabela da etapa 01.")
}

current_alias <- species_base[, .(
  GBIF_name_normalized = normalize_species_key(Species_current),
  WCVP_taxon_id,
  Analysis_species_key,
  Species_current,
  Genus_current,
  alias_source = "current_WCVP_name",
  alias_priority = 1L
)]

palmtraits_alias <- species_base[
  !is.na(PalmTraits_name),
  .(
    GBIF_name_normalized = normalize_species_key(PalmTraits_name),
    WCVP_taxon_id,
    Analysis_species_key,
    Species_current,
    Genus_current,
    alias_source = "PalmTraits_legacy_name",
    alias_priority = 2L
  )
]

wcvp_safe_names <- data.table::fread(
  WCVP_SAFE_NAME_DICTIONARY_FILE,
  encoding = "UTF-8",
  na.strings = c("", "NA")
)
required_safe_columns <- c(
  "supplied_name_key", "match_status", "matched_accepted_plant_name_id",
  "matched_accepted_taxon_name", "requires_author_or_manual_review"
)
missing_safe_columns <- setdiff(required_safe_columns, names(wcvp_safe_names))
if (length(missing_safe_columns) > 0L) {
  stop(
    "O dicionario WCVP seguro nao contem: ",
    paste(missing_safe_columns, collapse = ", "),
    call. = FALSE
  )
}
wcvp_safe_names[, WCVP_taxon_id := sub(
  "\\.0$", "", as.character(matched_accepted_plant_name_id)
)]
wcvp_safe_names[, manual_review_flag := to_integer01(
  requires_author_or_manual_review
)]
wcvp_safe_names[, GBIF_name_normalized := normalize_species_key(supplied_name_key)]
wcvp_safe_names <- wcvp_safe_names[
  manual_review_flag == 0L &
    match_status %chin% c(
      "resolved_current_accepted_name_priority",
      "resolved_unique_synonym"
    ) &
    !is.na(GBIF_name_normalized) & !is.na(WCVP_taxon_id)
]

safe_alias <- merge(
  wcvp_safe_names[, .(
    GBIF_name_normalized,
    WCVP_taxon_id,
    alias_source = paste0("WCVP_safe_", match_status),
    alias_priority = 3L
  )],
  species_base[, .(
    WCVP_taxon_id, Analysis_species_key, Species_current, Genus_current
  )],
  by = "WCVP_taxon_id",
  all = FALSE,
  sort = FALSE
)

alias_map <- rbindlist(
  list(current_alias, palmtraits_alias, safe_alias),
  use.names = TRUE,
  fill = TRUE
)
alias_map <- alias_map[!is.na(GBIF_name_normalized)]

# As duas decisoes documentadas tem precedencia total e substituem qualquer
# destino alternativo produzido por um homonimo sem autoria.
alias_map <- alias_map[
  !GBIF_name_normalized %chin% manual_alias$GBIF_name_normalized
]
alias_map <- rbindlist(list(manual_alias, alias_map), use.names = TRUE, fill = TRUE)

# Nomes atuais tem precedencia sobre aliases legados.
setorder(alias_map, GBIF_name_normalized, alias_priority)
alias_audit <- alias_map[, .(
  n_WCVP_targets = uniqueN(WCVP_taxon_id),
  targets = collapse_sorted(WCVP_taxon_id),
  current_names = collapse_sorted(Species_current)
), by = GBIF_name_normalized]

ambiguous_aliases <- alias_audit[n_WCVP_targets > 1L]
alias_map <- alias_map[
  !GBIF_name_normalized %in% ambiguous_aliases$GBIF_name_normalized
]
alias_map <- unique(alias_map, by = "GBIF_name_normalized")
alias_map[, alias_priority := NULL]

alias_dictionary_qa <- alias_map[, .(
  Alias_count = .N,
  WCVP_target_count = uniqueN(WCVP_taxon_id)
), by = alias_source][order(alias_source)]
write_csv_semicolon(
  alias_dictionary_qa,
  file.path(TABLES_DIR, "Table_QA_02_GBIF_alias_dictionary.csv")
)

write_csv_semicolon(
  ambiguous_aliases,
  file.path(TABLES_DIR, "Table_QA_02_ambiguous_GBIF_aliases.csv")
)

# ---- Etapa 1: leitura seletiva, filtros e taxonomia --------------------------

if (file.exists(checkpoint_1) && !FORCE_REBUILD) {
  message("[01/04] Retomando checkpoint GBIF filtrado e reconciliado...")
  stage1 <- readRDS(checkpoint_1)
  gbif_filtered <- stage1$data
  stage_counts <- stage1$counts
  gbif_unmatched_names <- stage1$unmatched_names
} else {
  message("[01/04] Lendo colunas essenciais do GBIF. Esta e a etapa mais pesada...")

  selected_columns <- c(
    "gbifID", "species", "genus", "taxonRank", "countryCode",
    "occurrenceStatus", "basisOfRecord", "establishmentMeans",
    "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
    "speciesKey", "taxonKey", "year", "issue"
  )

  gbif <- data.table::fread(
    GBIF_FILE,
    select = selected_columns,
    na.strings = c("", "NA"),
    showProgress = TRUE,
    encoding = "UTF-8"
  )

  raw_records <- nrow(gbif)
  gbif[, GBIF_name_normalized := normalize_species_key(species)]
  gbif[, occurrenceStatus := toupper(trimws(as.character(occurrenceStatus)))]
  gbif[, basisOfRecord := toupper(trimws(as.character(basisOfRecord)))]
  gbif[, establishmentMeans := toupper(trimws(as.character(establishmentMeans)))]
  gbif[, issue := toupper(as.character(issue))]
  gbif[, decimalLatitude := to_numeric_safe(decimalLatitude)]
  gbif[, decimalLongitude := to_numeric_safe(decimalLongitude)]
  gbif[, coordinateUncertaintyInMeters := to_numeric_safe(coordinateUncertaintyInMeters)]

  gbif <- merge(
    gbif, alias_map,
    by = "GBIF_name_normalized", all.x = TRUE, sort = FALSE
  )

  gbif_unmatched_names <- gbif[
    is.na(WCVP_taxon_id) & !is.na(GBIF_name_normalized),
    .(GBIF_record_count = .N),
    by = GBIF_name_normalized
  ][order(-GBIF_record_count)]

  matched_taxonomy <- sum(!is.na(gbif$WCVP_taxon_id))

  gbif <- gbif[
    !is.na(WCVP_taxon_id) &
      (is.na(occurrenceStatus) | occurrenceStatus == "" | occurrenceStatus == "PRESENT") &
      basisOfRecord %chin% ALLOWED_BASIS_OF_RECORD &
      is.finite(decimalLongitude) & is.finite(decimalLatitude) &
      between(decimalLongitude, -180, 180) &
      between(decimalLatitude, -90, 90) &
      !(decimalLongitude == 0 & decimalLatitude == 0) &
      (is.na(coordinateUncertaintyInMeters) |
         coordinateUncertaintyInMeters <= MAX_COORDINATE_UNCERTAINTY_M) &
      (is.na(establishmentMeans) |
         !grepl(EXPLICIT_NON_NATIVE_PATTERN, establishmentMeans)) &
      (is.na(issue) | !grepl(EXCLUDED_COORDINATE_ISSUE_PATTERN, issue))
  ]

  basic_filtered_records <- nrow(gbif)

  stage_counts <- data.table(
    Stage = c(
      "GBIF raw records",
      "Records matched to current, safe historical or documented manual WCVP name",
      "Records after status, basis, coordinate and uncertainty filters"
    ),
    Records = c(raw_records, matched_taxonomy, basic_filtered_records)
  )

  gbif_filtered <- gbif[, .(
    gbifID = as.character(gbifID),
    WCVP_taxon_id = as.character(WCVP_taxon_id),
    Analysis_species_key,
    Species_current,
    Genus_current,
    GBIF_name_normalized,
    alias_source,
    countryCode,
    decimalLongitude,
    decimalLatitude,
    coordinateUncertaintyInMeters,
    basisOfRecord,
    year
  )]

  saveRDS(
    list(
      data = gbif_filtered,
      counts = stage_counts,
      unmatched_names = gbif_unmatched_names
    ),
    checkpoint_1,
    compress = FALSE
  )
  rm(gbif)
  gc()
}

write_csv_semicolon(
  head(gbif_unmatched_names, 500L),
  file.path(TABLES_DIR, "Table_QA_02_top_unmatched_GBIF_names.csv")
)

# ---- Etapa 2: uma linha por especie x coordenada -----------------------------

if (file.exists(checkpoint_2) && !FORCE_REBUILD) {
  message("[02/04] Retomando checkpoint de coordenadas unicas...")
  unique_coordinates <- readRDS(checkpoint_2)
} else {
  message("[02/04] Removendo duplicacoes por especie e coordenada...")
  unique_coordinates <- gbif_filtered[, .(
    GBIF_records_at_coordinate = .N,
    countryCode = countryCode[which.max(!is.na(countryCode))],
    coordinate_uncertainty_missing = as.integer(all(is.na(coordinateUncertaintyInMeters))),
    coordinate_uncertainty_max_m = suppressWarnings(max(
      coordinateUncertaintyInMeters,
      na.rm = TRUE
    ))
  ), by = .(
    WCVP_taxon_id, Analysis_species_key, Species_current, Genus_current,
    decimalLongitude, decimalLatitude
  )]
  unique_coordinates[!is.finite(coordinate_uncertainty_max_m), coordinate_uncertainty_max_m := NA_real_]

  if (RUN_COORDINATECLEANER) {
    message("Executando CoordinateCleaner nas coordenadas unicas...")
    records_before_cleaner <- nrow(unique_coordinates)
    clean_flags <- CoordinateCleaner::clean_coordinates(
      x = as.data.frame(unique_coordinates),
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      species = "Species_current",
      countries = "countryCode",
      tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros"),
      value = "flagged",
      verbose = FALSE
    )
    clean_flags <- as.logical(clean_flags)
    if (length(clean_flags) != records_before_cleaner) {
      stop(
        "CoordinateCleaner retornou ", length(clean_flags),
        " flags para ", records_before_cleaner, " coordenadas.",
        call. = FALSE
      )
    }
    unique_coordinates <- unique_coordinates[
      which(!is.na(clean_flags) & clean_flags)
    ]
    message(
      "CoordinateCleaner: ", records_before_cleaner, " antes; ",
      nrow(unique_coordinates), " aprovadas; ",
      records_before_cleaner - nrow(unique_coordinates), " removidas."
    )
  }

  saveRDS(unique_coordinates, checkpoint_2, compress = FALSE)
  rm(gbif_filtered)
  gc()
}

# ---- Poligonos WGSRPD nivel 3 ------------------------------------------------

load_wgsrpd3 <- function() {
  wgsrpd3_map <- NULL

  if (requireNamespace("rWCVP", quietly = TRUE)) {
    map_environment <- new.env(parent = emptyenv())
    suppressWarnings(
      utils::data("wgsrpd3", package = "rWCVP", envir = map_environment)
    )
    if (exists("wgsrpd3", envir = map_environment, inherits = FALSE)) {
      wgsrpd3_map <- get("wgsrpd3", envir = map_environment)
    }
  }

  if (is.null(wgsrpd3_map)) {
    assert_files_exist(TDWG3_ZIP_FILE, "shapefile WGSRPD nivel 3")
    unzip_dir <- file.path(CHECKPOINT_DIR, "tdwg_level3_shp")
    dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(TDWG3_ZIP_FILE, exdir = unzip_dir)
    shp_files <- list.files(
      unzip_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE,
      ignore.case = TRUE
    )
    if (length(shp_files) == 0L) stop("Nenhum .shp encontrado em tdwg_level3_shp.zip.")
    wgsrpd3_map <- sf::st_read(shp_files[[1L]], quiet = TRUE)
  }

  if (!inherits(wgsrpd3_map, "sf")) stop("A camada WGSRPD L3 nao e um objeto sf.")
  code_candidates <- c(
    "LEVEL3_COD", "LEVEL3_CODE", "level3_cod", "level3_code",
    "area_code_l3", "LEVEL3"
  )
  code_column <- intersect(code_candidates, names(wgsrpd3_map))[1L]
  if (is.na(code_column)) {
    stop(
      "Coluna de codigo WGSRPD L3 nao encontrada. Colunas disponiveis: ",
      paste(names(wgsrpd3_map), collapse = ", ")
    )
  }
  wgsrpd3_map$area_code_l3 <- toupper(trimws(as.character(wgsrpd3_map[[code_column]])))
  wgsrpd3_map <- sf::st_make_valid(wgsrpd3_map)
  wgsrpd3_map <- sf::st_transform(wgsrpd3_map, 4326)
  wgsrpd3_map <- wgsrpd3_map[!sf::st_is_empty(wgsrpd3_map), c("area_code_l3")]
  wgsrpd3_map
}

# ---- Etapa 3: restricao a distribuicao nativa -------------------------------

if (file.exists(checkpoint_3) && !FORCE_REBUILD) {
  message("[03/04] Retomando checkpoint de coordenadas nativas...")
  native_coordinates <- readRDS(checkpoint_3)
} else {
  message("[03/04] Associando coordenadas a WGSRPD L3 e retendo areas nativas...")
  wgsrpd3_map <- load_wgsrpd3()

  points_sf <- sf::st_as_sf(
    as.data.frame(unique_coordinates),
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326,
    remove = FALSE
  )

  joined <- suppressWarnings(sf::st_join(
    points_sf,
    wgsrpd3_map,
    join = sf::st_intersects,
    left = FALSE
  ))
  joined <- as.data.table(sf::st_drop_geometry(joined))
  joined[, area_code_l3 := toupper(trimws(as.character(area_code_l3)))]

  powo_native <- data.table::fread(
    POWO_NATIVE_FILE,
    sep = ";",
    select = c("plant_name_id", "area_code_l3"),
    encoding = "UTF-8",
    na.strings = c("", "NA")
  )
  powo_native[, WCVP_taxon_id := sub("\\.0$", "", as.character(plant_name_id))]
  powo_native[, area_code_l3 := toupper(trimws(as.character(area_code_l3)))]
  powo_native <- unique(powo_native[, .(WCVP_taxon_id, area_code_l3)])

  native_coordinates <- merge(
    joined,
    powo_native,
    by = c("WCVP_taxon_id", "area_code_l3"),
    all = FALSE,
    sort = FALSE
  )
  native_coordinates <- unique(
    native_coordinates,
    by = c("WCVP_taxon_id", "decimalLongitude", "decimalLatitude")
  )

  if (nrow(native_coordinates) == 0L) {
    stop("Nenhuma coordenada permaneceu apos a restricao a distribuicao nativa.")
  }
  saveRDS(native_coordinates, checkpoint_3, compress = FALSE)
  rm(points_sf, joined, powo_native, wgsrpd3_map)
  gc()
}

# ---- Etapa 4: extracao GHS-POP e deduplicacao por celula --------------------

if (file.exists(checkpoint_4) && !FORCE_REBUILD) {
  message("[04/04] Retomando checkpoint especie x celula GHS...")
  species_cells <- readRDS(checkpoint_4)
} else {
  message("[04/04] Projetando pontos e extraindo GHS-POP 2025...")
  ghs_pop <- terra::rast(GHS_POP_FILE)
  if (terra::nlyr(ghs_pop) != 1L) {
    warning("O GHS-POP possui mais de uma camada; sera usada a primeira.")
    ghs_pop <- ghs_pop[[1L]]
  }

  points_vect <- terra::vect(
    native_coordinates,
    geom = c("decimalLongitude", "decimalLatitude"),
    crs = "EPSG:4326",
    keepgeom = TRUE
  )
  points_mollweide <- terra::project(points_vect, terra::crs(ghs_pop))
  xy <- terra::crds(points_mollweide)
  extracted <- terra::extract(ghs_pop, points_mollweide, ID = FALSE)

  native_coordinates[, GHS_x_m := xy[, 1L]]
  native_coordinates[, GHS_y_m := xy[, 2L]]
  native_coordinates[, GHS_cell_1km := terra::cellFromXY(ghs_pop, xy)]
  native_coordinates[, GHS_POP_2025 := to_numeric_safe(extracted[[1L]])]
  native_coordinates[GHS_POP_2025 < 0, GHS_POP_2025 := NA_real_]
  native_coordinates[, GHS_grid10_x := floor(GHS_x_m / 10000)]
  native_coordinates[, GHS_grid10_y := floor(GHS_y_m / 10000)]

  # Uma especie contribui apenas uma vez por celula de 1 km.
  species_cells <- native_coordinates[
    !is.na(GHS_cell_1km),
    .(
      Analysis_species_key = Analysis_species_key[1L],
      Species_current = Species_current[1L],
      Genus_current = Genus_current[1L],
      area_code_l3 = collapse_sorted(area_code_l3),
      GBIF_records_in_cell = sum(GBIF_records_at_coordinate),
      GBIF_unique_coordinates_in_cell = .N,
      GHS_x_m = GHS_x_m[1L],
      GHS_y_m = GHS_y_m[1L],
      GHS_grid10_x = GHS_grid10_x[1L],
      GHS_grid10_y = GHS_grid10_y[1L],
      GHS_POP_2025 = GHS_POP_2025[which.max(!is.na(GHS_POP_2025))]
    ),
    by = .(WCVP_taxon_id, GHS_cell_1km)
  ]

  saveRDS(species_cells, checkpoint_4, compress = FALSE)
  rm(ghs_pop, points_vect, points_mollweide, extracted, xy)
  gc()
}

# ---- Resumo por especie ------------------------------------------------------

native_point_summary <- native_coordinates[, .(
  GBIF_native_records = sum(GBIF_records_at_coordinate),
  GBIF_native_unique_coordinates = .N,
  GBIF_native_WGSRPD_L3_count = uniqueN(area_code_l3),
  GBIF_native_WGSRPD_L3_codes = collapse_sorted(area_code_l3),
  GBIF_coordinate_uncertainty_missing_share = mean(coordinate_uncertainty_missing == 1L)
), by = .(WCVP_taxon_id)]

gbif_ghs_summary <- species_cells[, .(
  GHS_unique_1km_cells = uniqueN(GHS_cell_1km),
  GBIF_unique_10km_cells = uniqueN(paste(GHS_grid10_x, GHS_grid10_y, sep = "_")),
  GHS_cells_with_population_data = sum(!is.na(GHS_POP_2025)),
  GHS_POP_2025_sum_unique_1km_cells = sum(GHS_POP_2025, na.rm = TRUE),
  GHS_POP_2025_mean = mean(GHS_POP_2025, na.rm = TRUE),
  GHS_POP_2025_median = median(GHS_POP_2025, na.rm = TRUE),
  GHS_POP_2025_p75 = as.numeric(quantile(GHS_POP_2025, 0.75, na.rm = TRUE, names = FALSE)),
  GHS_POP_2025_max = max(GHS_POP_2025, na.rm = TRUE),
  GHS_POP_2025_mean_log1p = mean(log1p(GHS_POP_2025), na.rm = TRUE),
  GHS_POP_2025_median_log1p = median(log1p(GHS_POP_2025), na.rm = TRUE),
  GHS_POP_2025_share_nonzero = mean(GHS_POP_2025 > 0, na.rm = TRUE)
), by = .(WCVP_taxon_id)]

numeric_summary_columns <- setdiff(names(gbif_ghs_summary), "WCVP_taxon_id")
for (column in numeric_summary_columns) {
  set(gbif_ghs_summary, which(!is.finite(gbif_ghs_summary[[column]])), column, NA_real_)
}

gbif_ghs_summary <- merge(
  native_point_summary,
  gbif_ghs_summary,
  by = "WCVP_taxon_id",
  all = TRUE,
  sort = FALSE
)
gbif_ghs_summary[, GBIF_GHS_data_available := as.integer(GHS_unique_1km_cells > 0)]

write_csv_semicolon(gbif_ghs_summary, GBIF_GHS_SUMMARY_OUTPUT)
write_csv_semicolon(
  species_cells,
  file.path(DATA_PROCESSED_DIR, "02_species_GHS_1km_cells.csv")
)

species_integrated <- merge(
  species_base,
  gbif_ghs_summary,
  by = "WCVP_taxon_id",
  all.x = TRUE,
  sort = FALSE
)

zero_columns <- c(
  "GBIF_native_records", "GBIF_native_unique_coordinates",
  "GBIF_native_WGSRPD_L3_count", "GHS_unique_1km_cells",
  "GBIF_unique_10km_cells", "GHS_cells_with_population_data",
  "GBIF_GHS_data_available"
)
zero_columns <- intersect(zero_columns, names(species_integrated))
for (column in zero_columns) {
  set(species_integrated, which(is.na(species_integrated[[column]])), column, 0)
}

write_csv_semicolon(
  species_integrated,
  file.path(DATA_PROCESSED_DIR, "02_species_analysis_with_GBIF_GHS.csv")
)
saveRDS(
  species_integrated,
  file.path(DATA_PROCESSED_DIR, "02_species_analysis_with_GBIF_GHS.rds"),
  compress = "xz"
)

qa <- rbindlist(list(
  stage_counts,
  data.table(
    Stage = c(
      "Unique species x coordinates after CoordinateCleaner",
      "Unique native species x coordinates after WCVP restriction",
      "WCVP species with native GBIF coordinates",
      "Unique species x GHS 1-km cells",
      "WCVP species with GHS-POP data"
    ),
    Records = c(
      nrow(unique_coordinates),
      nrow(native_coordinates),
      uniqueN(native_coordinates$WCVP_taxon_id),
      nrow(species_cells),
      sum(gbif_ghs_summary$GBIF_GHS_data_available == 1L, na.rm = TRUE)
    )
  )
), fill = TRUE)

write_csv_semicolon(qa, file.path(TABLES_DIR, "Table_QA_02_GBIF_GHS_pipeline.csv"))
write_session_info(file.path(LOGS_DIR, "02_sessionInfo.txt"))

elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "mins"))
log_lines <- c(
  "ARTIGO ARECACEAE — ETAPA 02",
  paste("Inicio:", format(started_at, tz = "UTC")),
  paste("Fim:", format(Sys.time(), tz = "UTC")),
  paste("Duracao_min:", round(elapsed, 2)),
  paste("MAX_COORDINATE_UNCERTAINTY_M:", MAX_COORDINATE_UNCERTAINTY_M),
  paste("RUN_COORDINATECLEANER:", RUN_COORDINATECLEANER),
  paste("Species_with_native_GBIF:", uniqueN(native_coordinates$WCVP_taxon_id)),
  paste("Species_with_GHS_POP:", sum(gbif_ghs_summary$GBIF_GHS_data_available == 1L, na.rm = TRUE)),
  paste("Output:", GBIF_GHS_SUMMARY_OUTPUT)
)
writeLines(log_lines, file.path(LOGS_DIR, "02_execution_log.txt"), useBytes = TRUE)

message("Concluido.")
message("Resumo GBIF-GHS: ", GBIF_GHS_SUMMARY_OUTPUT)
message(
  "Tabela integrada: ",
  file.path(DATA_PROCESSED_DIR, "02_species_analysis_with_GBIF_GHS.csv")
)
print(qa)
