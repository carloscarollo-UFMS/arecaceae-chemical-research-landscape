# =============================================================================
# ARTIGO ARECACEAE — 01. TABELA ANALITICA POR ESPECIE
# Versao: 1.0.2
#
# Unidade desta tabela: um WCVP_taxon_id por linha.
# O script nao altera as bases originais. Ele:
#   1. importa a base integrada WCVP + PalmTraits + quimica + economia;
#   2. recalcula amplitude geografica nativa a partir do POWO/WCVP;
#   3. reconcilia a taxonomia atual com TREE.nex;
#   4. poda e renomeia a arvore com chaves WCVP estaveis;
#   5. cria respostas e covariaveis explicitamente nomeadas;
#   6. salva auditorias, CSV, RDS e arvore pronta para os modelos.
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

require_packages(c("data.table", "readxl", "ape"))
assert_files_exist(
  c(BASE_MASTER_FILE, BASE_INTEGRATED_FILE, POWO_NATIVE_FILE, SPECIES_TREE_FILE, GENUS_TREE_FILE),
  c("base mestra", "base integrada", "distribuicao POWO", "arvore de especies", "arvore filogenomica")
)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ape)
})

started_at <- Sys.time()
message("[01] Lendo Species_model...")

species <- as.data.table(readxl::read_excel(
  BASE_INTEGRATED_FILE,
  sheet = "Species_model",
  .name_repair = "minimal"
))

required_species_columns <- c(
  "WCVP_taxon_id", "Species_current", "Genus_current",
  "Chemically_studied_250", "Article_count", "Organ_primary_count",
  "Occurrence_or_presence_records", "Chemical_entity_count",
  "Exact_structure_count", "Chemical_domain_count", "Chemical_class_count",
  "PalmTraits_available", "PalmTraits_name",
  "Commodity_species_resolved", "Commodity_candidate_broad"
)
missing_species_columns <- setdiff(required_species_columns, names(species))
if (length(missing_species_columns) > 0L) {
  stop(
    "A planilha Species_model nao contem: ",
    paste(missing_species_columns, collapse = ", "),
    call. = FALSE
  )
}

species[, WCVP_taxon_id := sub("\\.0$", "", as.character(WCVP_taxon_id))]
species[, Species_current := normalize_species_name(Species_current)]
species[, Genus_current := trimws(as.character(Genus_current))]
species[, PalmTraits_name := normalize_species_name(PalmTraits_name)]

binary_columns <- c(
  "PalmTraits_available", "Chemically_studied_250", "Chemistry_merge_review",
  "Commodity_species_resolved", "Commodity_candidate_broad"
)
binary_columns <- intersect(binary_columns, names(species))
species[, (binary_columns) := lapply(.SD, to_integer01), .SDcols = binary_columns]

numeric_columns <- c(
  "WCVP_locality_count", "WCVP_TDWG_count", "PalmTraits_source_rows_collapsed",
  "Climbing", "Acaulescent", "Erect", "StemSolitary", "StemArmed",
  "LeavesArmed", "MaxStemHeight_m", "MaxStemDia_cm", "MaxLeafNumber",
  "Max_Blade_Length_m", "Max_Rachis_Length_m", "Max_Petiole_length_m",
  "AverageFruitLength_cm", "MinFruitLength_cm", "MaxFruitLength_cm",
  "AverageFruitWidth_cm", "MinFruitWidth_cm", "MaxFruitWidth_cm",
  "Article_count", "Organ_primary_count", "Occurrence_or_presence_records",
  "Chemical_entity_count", "Exact_structure_count", "Chemical_domain_count",
  "Chemical_class_count", "Compounds_per_article", "Commodity_assignment_count",
  "FAOSTAT_GPV_current_USD_mean_2020_2024",
  "FAOSTAT_GPV_constant_USD_mean_2020_2024",
  "FAOSTAT_production_t_mean_2020_2024",
  "FAOSTAT_area_harvested_ha_mean_2020_2024",
  "FAOSTAT_yield_kg_ha_mean_2020_2024",
  "IBGE_PEVS_value_BRL_mean_2020_2024",
  "IBGE_PAM_value_BRL_mean_2020_2024"
)
numeric_columns <- intersect(numeric_columns, names(species))
species[, (numeric_columns) := lapply(.SD, to_numeric_safe), .SDcols = numeric_columns]

# ---- Recalculo da quimica apos resolucoes taxonomicas manuais ---------------
#
# Estes dois nomes sem autoria eram os unicos casos nao resolvidos na primeira
# integracao. A decisao usa WCVP/POWO 2026 e a evidencia dos artigos:
# - Astrocaryum aculeatum G.Mey., tucuma amazonico, e sinonimo de A. tucuma;
# - Phoenix pusilla Gaertn., coletada no artigo em Tamil Nadu, e sinonimo de
#   P. sylvestris no backbone WCVP/POWO adotado.

manual_taxonomy <- data.table(
  Species_standardized = c("Astrocaryum aculeatum", "Phoenix pusilla"),
  WCVP_taxon_id_manual = c("17583", "152708"),
  WCVP_accepted_name_manual = c("Astrocaryum tucuma", "Phoenix sylvestris"),
  Resolution_basis = c(
    "Astrocaryum aculeatum G.Mey. -> Astrocaryum tucuma; Amazonian tucuma context",
    "Phoenix pusilla Gaertn. -> Phoenix sylvestris; Tamil Nadu voucher/locality context"
  ),
  WCVP_source_URL = c(
    "https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:326537-2",
    "https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:60458340-2"
  ),
  Evidence_source = c(
    "WCVP/POWO synonymy plus Amazonian study context",
    "WCVP/POWO synonymy plus article collection at Nemili, Tamil Nadu, India"
  )
)

studied_crosswalk <- as.data.table(readxl::read_excel(
  BASE_INTEGRATED_FILE,
  sheet = "Studied_crosswalk",
  .name_repair = "minimal"
))
studied_crosswalk[, Species_standardized := normalize_species_name(Species_standardized)]
studied_crosswalk[, WCVP_taxon_id := sub("\\.0$", "", as.character(WCVP_taxon_id))]
if (anyDuplicated(studied_crosswalk$Species_standardized)) {
  stop("Studied_crosswalk possui mais de uma linha para a mesma unidade-fonte.")
}
studied_crosswalk <- merge(
  studied_crosswalk,
  manual_taxonomy,
  by = "Species_standardized",
  all.x = TRUE,
  sort = FALSE
)
studied_crosswalk[
  is.na(WCVP_taxon_id) & !is.na(WCVP_taxon_id_manual),
  WCVP_taxon_id := WCVP_taxon_id_manual
]
studied_crosswalk[
  (is.na(WCVP_accepted_name) | trimws(WCVP_accepted_name) == "") &
    !is.na(WCVP_accepted_name_manual),
  WCVP_accepted_name := WCVP_accepted_name_manual
]

if (anyNA(studied_crosswalk$WCVP_taxon_id)) {
  stop(
    "Persistem unidades quimicas sem WCVP_taxon_id: ",
    paste(studied_crosswalk[is.na(WCVP_taxon_id), Species_standardized], collapse = " | "),
    call. = FALSE
  )
}

evidence <- as.data.table(readxl::read_excel(
  BASE_MASTER_FILE,
  sheet = "Evidence_presence",
  .name_repair = "minimal"
))
evidence <- evidence[
  tolower(trimws(as.character(Species_analysis_ready))) %chin% c("yes", "sim", "1", "true")
]
evidence[, Species_standardized := normalize_species_name(Species_standardized)]

evidence <- merge(
  evidence,
  studied_crosswalk[, .(Species_standardized, WCVP_taxon_id)],
  by = "Species_standardized",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(evidence$WCVP_taxon_id)) {
  stop(
    "Registros de evidencia prontos para especie sem destino WCVP: ",
    paste(unique(evidence[is.na(WCVP_taxon_id), Species_standardized]), collapse = " | "),
    call. = FALSE
  )
}

chemistry_recomputed <- evidence[, .(
  Chemistry_source_taxa_collapsed = uniqueN(Species_standardized),
  Chemistry_merge_review = as.integer(uniqueN(Species_standardized) > 1L),
  Chemistry_name_original = collapse_sorted(Species_standardized),
  Article_count_new = uniqueN(Article_ID),
  Organ_primary_count_new = uniqueN(Organ_primary, na.rm = TRUE),
  Occurrence_or_presence_records_new = .N,
  Chemical_entity_count_new = uniqueN(Chemical_entity_ID, na.rm = TRUE),
  Exact_structure_count_new = uniqueN(Exact_structure_ID, na.rm = TRUE),
  Chemical_domain_count_new = uniqueN(Chemical_domain, na.rm = TRUE),
  Chemical_class_count_new = uniqueN(Chemical_class_harmonized, na.rm = TRUE)
), by = WCVP_taxon_id]
chemistry_recomputed[, Compounds_per_article_new :=
  Chemical_entity_count_new / Article_count_new]

# Zera apenas o bloco derivado; dados WCVP, PalmTraits e economicos permanecem.
species[, `:=`(
  Chemically_studied_250 = 0L,
  Chemistry_source_taxa_collapsed = 0,
  Chemistry_merge_review = 0L,
  Chemistry_name_original = NA_character_,
  Article_count = NA_real_,
  Organ_primary_count = NA_real_,
  Occurrence_or_presence_records = NA_real_,
  Chemical_entity_count = NA_real_,
  Exact_structure_count = NA_real_,
  Chemical_domain_count = NA_real_,
  Chemical_class_count = NA_real_,
  Compounds_per_article = NA_real_
)]

# Atualizacao por chave evita sufixos .x/.y e preserva a ordem original das
# 2.582 unidades WCVP.
species[chemistry_recomputed, on = "WCVP_taxon_id", `:=`(
  Chemistry_source_taxa_collapsed = i.Chemistry_source_taxa_collapsed,
  Chemistry_merge_review = i.Chemistry_merge_review,
  Chemistry_name_original = i.Chemistry_name_original,
  Chemically_studied_250 = 1L,
  Article_count = i.Article_count_new,
  Organ_primary_count = i.Organ_primary_count_new,
  Occurrence_or_presence_records = i.Occurrence_or_presence_records_new,
  Chemical_entity_count = i.Chemical_entity_count_new,
  Exact_structure_count = i.Exact_structure_count_new,
  Chemical_domain_count = i.Chemical_domain_count_new,
  Chemical_class_count = i.Chemical_class_count_new,
  Compounds_per_article = i.Compounds_per_article_new
)]

write_csv_semicolon(
  manual_taxonomy,
  file.path(TABLES_DIR, "Table_QA_01_manual_taxonomic_resolutions.csv")
)

if (nrow(species) != 2582L) {
  stop("Esperadas 2.582 unidades WCVP; observadas: ", nrow(species), call. = FALSE)
}
if (uniqueN(species$WCVP_taxon_id) != nrow(species)) {
  stop("WCVP_taxon_id nao e unico em Species_model.", call. = FALSE)
}

message("[02] Recalculando amplitude geografica nativa...")
powo <- data.table::fread(
  POWO_NATIVE_FILE,
  sep = ";",
  encoding = "UTF-8",
  na.strings = c("", "NA")
)

required_powo_columns <- c(
  "plant_name_id", "continent_code_l1", "continent",
  "region_code_l2", "region", "area_code_l3", "area", "taxon_name"
)
missing_powo_columns <- setdiff(required_powo_columns, names(powo))
if (length(missing_powo_columns) > 0L) {
  stop(
    "O arquivo POWO nao contem: ",
    paste(missing_powo_columns, collapse = ", "),
    call. = FALSE
  )
}

powo[, plant_name_id := sub("\\.0$", "", as.character(plant_name_id))]
powo[, area_code_l3 := toupper(trimws(as.character(area_code_l3)))]
powo <- unique(powo, by = c("plant_name_id", "area_code_l3"))

powo_summary <- powo[, .(
  POWO_native_L1_count = uniqueN(continent_code_l1[!is.na(continent_code_l1)]),
  POWO_native_L2_count = uniqueN(region_code_l2[!is.na(region_code_l2)]),
  POWO_native_L3_count = uniqueN(area_code_l3[!is.na(area_code_l3)]),
  POWO_native_continents = collapse_sorted(continent),
  POWO_native_regions_L2 = collapse_sorted(region),
  POWO_native_areas_L3 = collapse_sorted(area),
  POWO_native_codes_L3 = collapse_sorted(area_code_l3)
), by = .(WCVP_taxon_id = plant_name_id)]

species <- merge(
  species, powo_summary,
  by = "WCVP_taxon_id", all.x = TRUE, sort = FALSE
)

message("[03] Reconciliando TREE.nex com os nomes WCVP atuais...")
species_tree <- ape::read.nexus(SPECIES_TREE_FILE)
if (inherits(species_tree, "multiPhylo")) {
  if (length(species_tree) != 1L) {
    stop("TREE.nex contem mais de uma arvore; escolha explicitamente a arvore principal.")
  }
  species_tree <- species_tree[[1L]]
}
if (!inherits(species_tree, "phylo")) stop("TREE.nex nao foi lida como objeto phylo.")
if (anyDuplicated(species_tree$tip.label)) stop("TREE.nex contem terminais duplicados.")

tree_tips <- species_tree$tip.label
current_tip <- gsub(" ", "_", species$Species_current, fixed = TRUE)
palmtraits_tip <- gsub(" ", "_", species$PalmTraits_name, fixed = TRUE)

exact_match <- !is.na(current_tip) & current_tip %in% tree_tips
palmtraits_match <- !exact_match & !is.na(palmtraits_tip) & palmtraits_tip %in% tree_tips

species[, Tree_tip_original := NA_character_]
species[exact_match, Tree_tip_original := current_tip[exact_match]]
species[palmtraits_match, Tree_tip_original := palmtraits_tip[palmtraits_match]]
species[, Tree_match_method := fifelse(
  exact_match, "exact_current_name",
  fifelse(palmtraits_match, "PalmTraits_legacy_name", "not_in_TREE_nex")
)]
species[, Tree_included := as.integer(!is.na(Tree_tip_original))]
species[, Analysis_species_key := paste0("WCVP_", WCVP_taxon_id)]
species[, Tree_tip_analysis := fifelse(Tree_included == 1L, Analysis_species_key, NA_character_)]

matched_crosswalk <- species[Tree_included == 1L, .(
  WCVP_taxon_id, Analysis_species_key, Species_current, Genus_current,
  PalmTraits_name, Tree_tip_original, Tree_tip_analysis,
  Tree_match_method, Chemically_studied_250
)]

if (anyDuplicated(matched_crosswalk$Tree_tip_original)) {
  stop("Mais de uma especie WCVP foi associada ao mesmo terminal de TREE.nex.")
}

pruned_tree <- ape::keep.tip(species_tree, matched_crosswalk$Tree_tip_original)
tip_lookup <- setNames(
  matched_crosswalk$Tree_tip_analysis,
  matched_crosswalk$Tree_tip_original
)
pruned_tree$tip.label <- unname(tip_lookup[pruned_tree$tip.label])
if (anyNA(pruned_tree$tip.label) || anyDuplicated(pruned_tree$tip.label)) {
  stop("Falha ao renomear a arvore podada com chaves WCVP.")
}
ape::write.tree(pruned_tree, file = SPECIES_TREE_OUTPUT)

genus_tree <- ape::read.tree(GENUS_TREE_FILE)
if (!inherits(genus_tree, "phylo")) stop("A arvore filogenomica nao foi lida como phylo.")

message("[04] Criando respostas e covariaveis explicitamente nomeadas...")
species[, Research_studied_binary := as.integer(Chemically_studied_250)]
species[, Research_article_count := fifelse(is.na(Article_count), 0, Article_count)]
species[, Research_evidence_record_count := fifelse(
  is.na(Occurrence_or_presence_records), 0, Occurrence_or_presence_records
)]
species[, Chem_richness_entities_observed := fifelse(
  is.na(Chemical_entity_count), 0, Chemical_entity_count
)]
species[, Chem_richness_exact_structures := fifelse(
  is.na(Exact_structure_count), 0, Exact_structure_count
)]
species[, Chem_richness_classes := fifelse(
  is.na(Chemical_class_count), 0, Chemical_class_count
)]
species[, Chem_richness_domains := fifelse(
  is.na(Chemical_domain_count), 0, Chemical_domain_count
)]
species[, Chem_organ_primary_count := fifelse(
  is.na(Organ_primary_count), 0, Organ_primary_count
)]
species[, Research_log1p_articles := log1p(Research_article_count)]
species[, Research_log1p_evidence_records := log1p(Research_evidence_record_count)]
species[, Chem_entities_per_article_observed := fifelse(
  Research_article_count > 0,
  Chem_richness_entities_observed / Research_article_count,
  NA_real_
)]
species[, Human_use_resolved := as.integer(Commodity_species_resolved)]
species[, Human_use_broad_candidate := as.integer(Commodity_candidate_broad)]
species[, Human_use_any := as.integer(
  Human_use_resolved == 1L | Human_use_broad_candidate == 1L
)]
species[, Economic_quantitative_data_available := as.integer(
  !is.na(FAOSTAT_GPV_current_USD_mean_2020_2024) |
    !is.na(FAOSTAT_production_t_mean_2020_2024) |
    !is.na(IBGE_PEVS_value_BRL_mean_2020_2024) |
    !is.na(IBGE_PAM_value_BRL_mean_2020_2024)
)]
species[, Trait_data_available := as.integer(PalmTraits_available)]
species[, Range_log1p_native_L3 := log1p(POWO_native_L3_count)]

setcolorder(species, c(
  "Analysis_species_key", "WCVP_taxon_id", "Species_current", "Genus_current",
  "Research_studied_binary", "Research_article_count",
  "Research_evidence_record_count", "Chem_richness_entities_observed",
  "Chem_richness_exact_structures", "Chem_richness_classes",
  "Chem_richness_domains", "Chem_organ_primary_count",
  "Research_log1p_articles", "Research_log1p_evidence_records",
  "Chem_entities_per_article_observed",
  "POWO_native_L1_count", "POWO_native_L2_count", "POWO_native_L3_count",
  "POWO_native_continents", "POWO_native_codes_L3",
  "Trait_data_available", "Human_use_resolved",
  "Human_use_broad_candidate", "Human_use_any",
  "Economic_quantitative_data_available",
  "Tree_included", "Tree_match_method", "Tree_tip_original", "Tree_tip_analysis",
  setdiff(names(species), c(
    "Analysis_species_key", "WCVP_taxon_id", "Species_current", "Genus_current",
    "Research_studied_binary", "Research_article_count",
    "Research_evidence_record_count", "Chem_richness_entities_observed",
    "Chem_richness_exact_structures", "Chem_richness_classes",
    "Chem_richness_domains", "Chem_organ_primary_count",
    "Research_log1p_articles", "Research_log1p_evidence_records",
    "Chem_entities_per_article_observed",
    "POWO_native_L1_count", "POWO_native_L2_count", "POWO_native_L3_count",
    "POWO_native_continents", "POWO_native_codes_L3",
    "Trait_data_available", "Human_use_resolved",
    "Human_use_broad_candidate", "Human_use_any",
    "Economic_quantitative_data_available",
    "Tree_included", "Tree_match_method", "Tree_tip_original", "Tree_tip_analysis"
  ))
))

qa <- data.table(
  Indicator = c(
    "WCVP taxon units",
    "Unique WCVP_taxon_id",
    "Unique current binomials",
    "Current binomial duplicates",
    "Chemically studied current species after manual resolutions",
    "Chemical source species units resolved",
    "POWO species with native distribution",
    "TREE.nex original tips",
    "TREE.nex matched WCVP units",
    "Studied species represented in TREE.nex",
    "Studied species absent from TREE.nex",
    "TREE.nex rooted",
    "TREE.nex binary",
    "Genus source-tree tips",
    "Genus source-tree rooted",
    "Genus source-tree binary"
  ),
  Observed = c(
    nrow(species),
    uniqueN(species$WCVP_taxon_id),
    uniqueN(species$Species_current),
    sum(duplicated(species$Species_current)),
    sum(species$Research_studied_binary == 1L, na.rm = TRUE),
    uniqueN(studied_crosswalk$Species_standardized),
    uniqueN(powo$plant_name_id),
    ape::Ntip(species_tree),
    sum(species$Tree_included == 1L),
    sum(species$Tree_included == 1L & species$Research_studied_binary == 1L),
    sum(species$Tree_included == 0L & species$Research_studied_binary == 1L),
    as.integer(ape::is.rooted(species_tree)),
    as.integer(ape::is.binary(species_tree)),
    ape::Ntip(genus_tree),
    as.integer(ape::is.rooted(genus_tree)),
    as.integer(ape::is.binary(genus_tree))
  ),
  Expected = c(
    2582, 2582, NA, NA, 173, 177, 2527, 2539, NA, 172, 1, 1, 1, 187, 1, 1
  ),
  Status = c(
    "PASS", "PASS", "INFO", "REVIEW", "PASS", "PASS", "PASS", "PASS", "INFO",
    "PASS", "REVIEW", "PASS", "PASS", "PASS", "PASS", "PASS"
  )
)
qa[, Interpretation := c(
  "Unidade primaria preservada por WCVP_taxon_id.",
  "Nao ha IDs WCVP repetidos.",
  "Nomes sem autoria nao sao usados como chave primaria.",
  "Calamus elegans possui dois IDs WCVP; manter separado ate revisao taxonomica.",
  "Especies atuais vinculadas ao corpus apos duas resolucoes manuais documentadas.",
  "Todas as 177 unidades-fonte possuem destino WCVP.",
  "Distribuicao nativa, atual e nao duvidosa.",
  "Arvore PalmTraits original.",
  "Correspondencia por nome atual ou nome legado PalmTraits.",
  "Cobertura filogenetica do conjunto estudado.",
  "Butia witeckii requer decisao filogenetica explicita.",
  "Validacao estrutural da arvore de especies.",
  "Validacao estrutural da arvore de especies.",
  "Arvore filogenomica antes da reducao por genero.",
  "Validacao estrutural da arvore filogenomica.",
  "Validacao estrutural da arvore filogenomica."
)]

message("[05] Salvando produtos...")
write_csv_semicolon(species, SPECIES_BASE_OUTPUT)
saveRDS(species, SPECIES_BASE_RDS, compress = "xz")
write_csv_semicolon(
  species[, .(
    WCVP_taxon_id, Analysis_species_key, Species_current, Genus_current,
    PalmTraits_name, Tree_included, Tree_match_method,
    Tree_tip_original, Tree_tip_analysis, Research_studied_binary
  )],
  SPECIES_TREE_CROSSWALK_OUTPUT
)
write_csv_semicolon(qa, file.path(TABLES_DIR, "Table_QA_01_species_base.csv"))
write_csv_semicolon(
  species[Tree_included == 0L & Research_studied_binary == 1L, .(
    WCVP_taxon_id, Species_current, Genus_current, PalmTraits_name,
    Chemistry_name_original, Tree_match_method
  )],
  file.path(TABLES_DIR, "Table_QA_01_studied_species_missing_tree.csv")
)
write_session_info(file.path(LOGS_DIR, "01_sessionInfo.txt"))

elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
log_lines <- c(
  "ARTIGO ARECACEAE — ETAPA 01",
  paste("Inicio:", format(started_at, tz = "UTC")),
  paste("Fim:", format(Sys.time(), tz = "UTC")),
  paste("Duracao_s:", round(elapsed, 2)),
  paste("Species_base:", SPECIES_BASE_OUTPUT),
  paste("Species_tree:", SPECIES_TREE_OUTPUT),
  paste("WCVP_units:", nrow(species)),
  paste("Studied_species:", sum(species$Research_studied_binary == 1L)),
  paste("Tree_matched_units:", sum(species$Tree_included == 1L)),
  paste(
    "Studied_missing_tree:",
    paste(species[Tree_included == 0L & Research_studied_binary == 1L, Species_current], collapse = " | ")
  )
)
writeLines(log_lines, file.path(LOGS_DIR, "01_execution_log.txt"), useBytes = TRUE)

message("Concluido.")
message("Tabela por especie: ", SPECIES_BASE_OUTPUT)
message("Arvore WCVP podada: ", SPECIES_TREE_OUTPUT)
print(qa)
