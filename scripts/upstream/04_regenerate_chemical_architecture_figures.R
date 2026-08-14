# =============================================================================
# ARTIGO ARECACEAE — ETAPA 04C
# Revisao cientifica e regeneracao das seis figuras da arquitetura quimica
# Versao: 1.1.5 | 2026-08-13
#
# Entradas: tabelas CSV curadas da propria pasta 04C (subpasta tables/)
# Saidas:  PDF vetorial via cairo_pdf, SVG editavel e PNG 600 dpi
#
# Este script NAO refaz a curadoria e NAO altera a base mestra ou as tabelas.
# Ele apenas reconstroi as seis figuras a partir dos resultados 04C congelados.
# =============================================================================

get_script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/")))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(path)) return(dirname(normalizePath(path, winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

SCRIPT_DIR <- get_script_directory()
STAGE04C_DIR <- SCRIPT_DIR

if (!dir.exists(file.path(STAGE04C_DIR, "tables"))) {
  candidates <- c(
    dirname(SCRIPT_DIR),
    file.path(SCRIPT_DIR, "stage04C_curadoria_quimica_aplicada_v1_0_1"),
    file.path(dirname(SCRIPT_DIR), "stage04C_curadoria_quimica_aplicada_v1_0_1")
  )
  valid <- candidates[
    vapply(candidates, function(x) dir.exists(file.path(x, "tables")), logical(1))
  ]
  if (length(valid) > 0L) STAGE04C_DIR <- valid[1L]
}

TABLES_DIR <- file.path(STAGE04C_DIR, "tables")
FIGURES_DIR <- file.path(STAGE04C_DIR, "figures")
LOGS_DIR <- file.path(STAGE04C_DIR, "logs")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "data.table", "ggplot2", "scales", "patchwork", "svglite", "ggrepel"
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
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(ggrepel)
})

if (!capabilities("cairo")) {
  stop(
    "Esta instalacao do R nao possui suporte Cairo. Atualize o R antes de gerar os PDFs.",
    call. = FALSE
  )
}

figure_font <- Sys.getenv(
  "ARECACEAE_FIGURE_FONT",
  unset = if (.Platform$OS.type == "windows") "Arial" else "sans"
)

input_files <- c(
  domain = "Table_04C_02_domain_architecture.csv",
  class = "Table_04C_03_class_architecture.csv",
  organ = "Table_04C_04_organ_architecture.csv",
  organ_class = "Table_04C_05_organ_class_incidence.csv",
  genus = "Table_04C_06_genus_architecture.csv",
  genus_class = "Table_04C_07_genus_class_incidence.csv",
  method_class = "Table_04C_10_method_class_association.csv"
)
input_paths <- setNames(
  file.path(TABLES_DIR, unname(input_files)),
  names(input_files)
)

if (any(!file.exists(input_paths))) {
  missing <- input_paths[!file.exists(input_paths)]
  stop(
    "Arquivos obrigatorios ausentes na subpasta tables/:\n- ",
    paste(basename(missing), collapse = "\n- "),
    call. = FALSE
  )
}

read_table <- function(path) {
  data.table::fread(
    path, sep = ";", encoding = "UTF-8", na.strings = c("", "NA"),
    check.names = FALSE
  )
}

domain_table <- read_table(input_paths[["domain"]])
class_table <- read_table(input_paths[["class"]])
organ_table <- read_table(input_paths[["organ"]])
organ_class <- read_table(input_paths[["organ_class"]])
genus_table <- read_table(input_paths[["genus"]])
genus_class <- read_table(input_paths[["genus_class"]])
method_table <- read_table(input_paths[["method_class"]])

required_columns <- list(
  domain = c(
    "Domain_analysis", "Chemical_entity_count", "Exact_structure_count",
    "Current_species_count", "Species_prevalence_current"
  ),
  class = c(
    "Domain_analysis", "Class_analysis", "Current_species_count",
    "Chemical_entity_count", "Article_count"
  ),
  organ = c(
    "Organ_primary", "Article_count", "Current_species_count",
    "Chemical_entity_count"
  ),
  organ_class = c(
    "Organ_primary", "Class_analysis", "Species_class_count",
    "Species_studied_in_organ", "Species_prevalence_within_organ"
  ),
  genus = c(
    "Genus_current", "Current_species_count", "Chemical_entity_count"
  ),
  genus_class = c(
    "Genus_current", "Class_analysis", "Species_class_count",
    "Studied_species_in_genus", "Species_prevalence_within_genus"
  ),
  method_class = c(
    "Class_analysis", "Analytical_family",
    "Percent_class_articles_using_method"
  )
)

tables <- list(
  domain = domain_table,
  class = class_table,
  organ = organ_table,
  organ_class = organ_class,
  genus = genus_table,
  genus_class = genus_class,
  method_class = method_table
)

for (name in names(required_columns)) {
  absent <- setdiff(required_columns[[name]], names(tables[[name]]))
  if (length(absent) > 0L) {
    stop(
      input_files[[name]], " nao contem: ", paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
}

domain_palette <- c(
  "Alkaloids and nitrogen compounds" = "#B2473E",
  "Amino acids and peptides" = "#7C9446",
  "Carbohydrates" = "#C2A23A",
  "Lipids and lipid-like molecules" = "#4A7FA0",
  "Other primary metabolites" = "#66888A",
  "Other specialized metabolites" = "#91668A",
  "Phenylpropanoids and polyketides" = "#397A58",
  "Terpenoids and steroids" = "#D6803A",
  "Vitamins and cofactors" = "#C58B9E"
)

unknown_domains <- setdiff(unique(class_table$Domain_analysis), names(domain_palette))
if (length(unknown_domains) > 0L) {
  extra_colours <- grDevices::hcl.colors(length(unknown_domains), "Dark 3")
  names(extra_colours) <- unknown_domains
  domain_palette <- c(domain_palette, extra_colours)
}

theme_article <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size, base_family = figure_font) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3),
      plot.subtitle = element_text(colour = "#4E5555", size = base_size),
      plot.caption = element_text(
        colour = "#4E5555", size = base_size - 1, hjust = 0
      ),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 18, 10, 10)
    )
}

save_editable_figure <- function(plot, filename, width, height) {
  pdf_path <- file.path(FIGURES_DIR, paste0(filename, ".pdf"))
  svg_path <- file.path(FIGURES_DIR, paste0(filename, ".svg"))
  png_path <- file.path(FIGURES_DIR, paste0(filename, ".png"))

  grDevices::cairo_pdf(
    filename = pdf_path,
    width = width,
    height = height,
    family = figure_font,
    bg = "white",
    onefile = FALSE
  )
  print(plot)
  grDevices::dev.off()

  svglite::svglite(
    file = svg_path,
    width = width,
    height = height,
    bg = "white",
    system_fonts = list(sans = figure_font)
  )
  print(plot)
  grDevices::dev.off()

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      filename = png_path,
      plot = plot,
      device = ragg::agg_png,
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      bg = "white"
    )
  } else {
    ggplot2::ggsave(
      filename = png_path,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      bg = "white"
    )
  }

  data.table(
    Figure = filename,
    PDF = pdf_path,
    SVG = svg_path,
    PNG = png_path,
    PDF_bytes = file.info(pdf_path)$size,
    SVG_bytes = file.info(svg_path)$size,
    PNG_bytes = file.info(png_path)$size,
    PDF_device = "grDevices::cairo_pdf",
    SVG_device = "svglite::svglite",
    Font_family = figure_font
  )
}

message("[01/06] Revisando a arquitetura por dominio...")

domain_plot <- copy(domain_table)
domain_plot[, Annotation_only_count :=
  Chemical_entity_count - Exact_structure_count]
setorder(domain_plot, Chemical_entity_count)
domain_plot[, Domain_analysis := factor(
  Domain_analysis, levels = Domain_analysis
)]

domain_resolution <- melt(
  domain_plot,
  id.vars = c("Domain_analysis", "Chemical_entity_count"),
  measure.vars = c("Annotation_only_count", "Exact_structure_count"),
  variable.name = "Resolution",
  value.name = "Entity_count"
)
domain_resolution[, Resolution := factor(
  Resolution,
  levels = c("Annotation_only_count", "Exact_structure_count"),
  labels = c("Annotation only", "Exact structure")
)]

p_domain_entities <- ggplot(
  domain_resolution,
  aes(Domain_analysis, Entity_count, fill = Resolution)
) +
  geom_col(width = 0.68) +
  geom_text(
    data = domain_plot,
    aes(Domain_analysis, Chemical_entity_count, label = Chemical_entity_count),
    inherit.aes = FALSE,
    hjust = -0.14,
    size = 3.1,
    family = figure_font
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c(
    "Annotation only" = "#CAD3D4",
    "Exact structure" = "#176B66"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "A  Chemical diversity",
    x = NULL,
    y = "Curated chemical entities",
    fill = "Identification"
  ) +
  theme_article() +
  theme(legend.position = "bottom")

p_domain_species <- ggplot(
  domain_plot,
  aes(Domain_analysis, Species_prevalence_current)
) +
  geom_col(width = 0.68, fill = "#D6803A") +
  geom_text(
    aes(label = paste0(Current_species_count, "/173")),
    hjust = -0.14,
    size = 3.1,
    family = figure_font
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.20))
  ) +
  labs(
    title = "B  Taxonomic breadth",
    x = NULL,
    y = "Studied species reporting the domain"
  ) +
  theme_article() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

figure_01 <- p_domain_entities + p_domain_species +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(
    title = "Documented chemical architecture of Arecaceae",
    caption = paste(
      "Entity richness and taxonomic breadth are distinct dimensions;",
      "exact structures are a subset of curated entities."
    ),
    theme = theme(
      plot.title = element_text(
        family = figure_font, face = "bold", size = 15, hjust = 0
      ),
      plot.caption = element_text(
        family = figure_font, colour = "#4E5555", size = 9, hjust = 0
      )
    )
  )

message("[02/06] Recuperando valores e cobertura por orgao...")

organ_plot <- copy(organ_table)
setorder(organ_plot, Article_count)
organ_plot[, Organ_primary := factor(
  Organ_primary, levels = Organ_primary
)]

p_organ_articles <- ggplot(
  organ_plot, aes(Organ_primary, Article_count)
) +
  geom_col(width = 0.68, fill = "#71B99C") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = "A  Literature effort",
    x = NULL,
    y = "Informative articles"
  ) +
  theme_article()

p_organ_entities <- ggplot(
  organ_plot, aes(Organ_primary, Chemical_entity_count)
) +
  geom_col(width = 0.68, fill = "#286C52") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "B  Observed chemical richness",
    x = NULL,
    y = "Curated chemical entities"
  ) +
  theme_article() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

figure_02 <- p_organ_articles + p_organ_entities +
  plot_layout(widths = c(1, 1.32)) +
  plot_annotation(
    title = "Sampling effort and observed chemical richness by plant organ",
    caption = paste(
      "Counts are descriptive and are not corrected for unequal numbers",
      "of articles or studied species."
    ),
    theme = theme(
      plot.title = element_text(
        family = figure_font, face = "bold", size = 15, hjust = 0
      ),
      plot.caption = element_text(
        family = figure_font, colour = "#4E5555", size = 9, hjust = 0
      )
    )
  )

message("[03/06] Restringindo a comparacao a orgaos interpretaveis...")

eligible_organs <- c(
  "Fruit", "Seed", "Leaf", "Flower/inflorescence", "Stem/shoot", "Root"
)
eligible_organs <- eligible_organs[
  eligible_organs %chin% organ_table$Organ_primary
]

top_classes_15 <- head(
  class_table[order(-Current_species_count, -Chemical_entity_count)],
  15L
)$Class_analysis

organ_grid <- CJ(
  Organ_primary = eligible_organs,
  Class_analysis = top_classes_15,
  unique = TRUE
)

organ_heat <- merge(
  organ_grid,
  organ_class[
    Organ_primary %chin% eligible_organs &
      Class_analysis %chin% top_classes_15,
    .(
      Organ_primary,
      Class_analysis,
      Species_class_count,
      Species_studied_in_organ,
      Species_prevalence_within_organ
    )
  ],
  by = c("Organ_primary", "Class_analysis"),
  all.x = TRUE
)

organ_n <- setNames(
  organ_table$Current_species_count,
  organ_table$Organ_primary
)
organ_heat[is.na(Species_class_count), Species_class_count := 0L]
organ_heat[is.na(Species_studied_in_organ),
  Species_studied_in_organ := organ_n[Organ_primary]]
organ_heat[is.na(Species_prevalence_within_organ),
  Species_prevalence_within_organ := 0]
organ_heat[, Cell_label := fifelse(
  Species_prevalence_within_organ > 0,
  percent(Species_prevalence_within_organ, accuracy = 1),
  ""
)]
organ_heat[, `:=`(
  Organ_label = factor(
    paste0(Organ_primary, " (n=", Species_studied_in_organ, " spp.)"),
    levels = paste0(
      rev(eligible_organs), " (n=", organ_n[rev(eligible_organs)], " spp.)"
    )
  ),
  Class_analysis = factor(Class_analysis, levels = top_classes_15)
)]

figure_03 <- ggplot(
  organ_heat,
  aes(Class_analysis, Organ_label, fill = Species_prevalence_within_organ)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(
    aes(label = Cell_label),
    size = 2.35,
    family = figure_font
  ) +
  scale_fill_gradientn(
    colours = c("#F4FBF7", "#C7E9D5", "#73C49C", "#16805D", "#064B37"),
    limits = c(0, 1),
    oob = squish,
    labels = percent_format(accuracy = 1),
    name = "Species prevalence\nwithin organ"
  ) +
  labs(
    title = "Chemical-class prevalence across plant organs represented by at least five species",
    subtitle = "Fifteen classes with the broadest species representation",
    caption = paste(
      "Mixed/unclear units and organs represented by fewer than five species",
      "were excluded from this comparison."
    ),
    x = NULL,
    y = NULL
  ) +
  theme_article() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1),
    legend.position = "right"
  )

message("[04/06] Inserindo nomes das classes na figura de amplitude...")

class_plot <- copy(class_table)
class_plot[, Rank_entities := frank(
  -Chemical_entity_count, ties.method = "min"
)]
class_plot[, Rank_species := frank(
  -Current_species_count, ties.method = "min"
)]
class_plot[, Joint_rank_score := Rank_entities + Rank_species]
setorder(
  class_plot,
  Joint_rank_score,
  Rank_entities,
  Rank_species,
  Class_analysis
)
class_plot[, Highlight_top20 := seq_len(.N) <= 20L]
class_context <- class_plot[Highlight_top20 == FALSE]
class_top20 <- class_plot[Highlight_top20 == TRUE]

figure_04 <- ggplot(
  mapping = aes(Current_species_count, Chemical_entity_count)
) +
  geom_point(
    data = class_context,
    aes(size = Article_count),
    colour = "#D3D8D8",
    alpha = 0.60,
    stroke = 0.25,
    show.legend = FALSE
  ) +
  geom_point(
    data = class_top20,
    aes(colour = Domain_analysis, size = Article_count),
    alpha = 0.88,
    stroke = 0.40
  ) +
  ggrepel::geom_text_repel(
    data = class_top20,
    aes(label = Class_analysis),
    family = figure_font,
    size = 3.20,
    seed = 20260813,
    max.overlaps = Inf,
    max.time = 15,
    force = 3,
    force_pull = 0.4,
    box.padding = 0.25,
    point.padding = 0.14,
    min.segment.length = 0,
    segment.colour = "#8C9393",
    segment.size = 0.25,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = domain_palette, drop = FALSE) +
  scale_size_continuous(range = c(2.2, 10)) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(
    title = "Taxonomic breadth versus observed chemical richness",
    subtitle = paste(
      "The 20 classes with the strongest joint ranks for entity richness",
      "and species breadth are highlighted and labelled."
    ),
    caption = paste(
      "The remaining 25 classes are shown in light grey for context; point",
      "area scales with the number of informative articles."
    ),
    x = "Current species represented",
    y = "Curated chemical entities",
    colour = "Chemical domain",
    size = "Articles"
  ) +
  theme_article() +
  theme(
    panel.grid.major.y = element_line(colour = "#E2E5E5", linewidth = 0.5),
    legend.position = "right"
  )

message("[05/06] Refinando a associacao entre metodos e classes...")

top_classes_richness_15 <- head(
  class_table[order(-Chemical_entity_count, Class_analysis)],
  15L
)$Class_analysis

core_methods <- c(
  "GC-based", "LC-based", "Mass spectrometry", "NMR/isolation"
)

method_grid <- CJ(
  Class_analysis = top_classes_richness_15,
  Analytical_family = core_methods,
  unique = TRUE
)

method_heat <- merge(
  method_grid,
  method_table[
    Class_analysis %chin% top_classes_richness_15 &
      Analytical_family %chin% core_methods,
    .(
      Class_analysis,
      Analytical_family,
      Percent_class_articles_using_method
    )
  ],
  by = c("Class_analysis", "Analytical_family"),
  all.x = TRUE
)

method_heat[is.na(Percent_class_articles_using_method),
  Percent_class_articles_using_method := 0]
method_heat[, `:=`(
  Class_analysis = factor(
    Class_analysis, levels = rev(top_classes_richness_15)
  ),
  Analytical_family = factor(
    Analytical_family,
    levels = core_methods,
    labels = c(
      "GC-based", "LC-based", "Mass\nspectrometry", "NMR /\nisolation"
    )
  )
)]

figure_05 <- ggplot(
  method_heat,
  aes(Analytical_family, Class_analysis,
      fill = Percent_class_articles_using_method)
) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colours = c("#F7F6FB", "#D6D1E8", "#A59ACD", "#6D5AA9", "#3E2E78"),
    limits = c(0, 1),
    oob = squish,
    labels = percent_format(accuracy = 1),
    name = "Fraction of\nclass articles"
  ) +
  labs(
    title = "Analytical-method context of the 15 richest chemical classes",
    subtitle = paste(
      "Methods are non-exclusive within articles; percentages across columns",
      "need not sum to 100%."
    ),
    x = NULL,
    y = NULL
  ) +
  theme_article() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9.5),
    legend.position = "right"
  )

message("[06/06] Corrigindo o denominador da comparacao entre generos...")

eligible_genera <- genus_table[
  !is.na(Genus_current) & Current_species_count >= 5L
][order(-Chemical_entity_count, Genus_current)]
eligible_genera <- head(eligible_genera, 15L)
top_genera <- eligible_genera$Genus_current

top_classes_breadth_15 <- head(
  class_table[order(-Current_species_count, -Chemical_entity_count)],
  15L
)$Class_analysis

genus_grid <- CJ(
  Genus_current = top_genera,
  Class_analysis = top_classes_breadth_15,
  unique = TRUE
)

genus_heat <- merge(
  genus_grid,
  genus_class[
    Genus_current %chin% top_genera &
      Class_analysis %chin% top_classes_breadth_15,
    .(
      Genus_current,
      Class_analysis,
      Species_class_count,
      Studied_species_in_genus,
      Species_prevalence_within_genus
    )
  ],
  by = c("Genus_current", "Class_analysis"),
  all.x = TRUE
)

genus_n <- setNames(
  eligible_genera$Current_species_count,
  eligible_genera$Genus_current
)
genus_heat[is.na(Species_class_count), Species_class_count := 0L]
genus_heat[is.na(Studied_species_in_genus),
  Studied_species_in_genus := genus_n[Genus_current]]
genus_heat[is.na(Species_prevalence_within_genus),
  Species_prevalence_within_genus := 0]
genus_heat[, `:=`(
  Genus_label = factor(
    paste0(Genus_current, " (n=", Studied_species_in_genus, " spp.)"),
    levels = paste0(
      rev(top_genera), " (n=", genus_n[rev(top_genera)], " spp.)"
    )
  ),
  Class_analysis = factor(
    Class_analysis, levels = top_classes_breadth_15
  )
)]

figure_06 <- ggplot(
  genus_heat,
  aes(Class_analysis, Genus_label, fill = Species_prevalence_within_genus)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  scale_fill_gradientn(
    colours = c("#F4FBF7", "#C7E9D5", "#73C49C", "#16805D", "#064B37"),
    limits = c(0, 1),
    oob = squish,
    labels = percent_format(accuracy = 1),
    name = "Proportion of studied\nspecies in genus"
  ) +
  labs(
    title = "Chemical-class representation across genera represented by at least five species",
    subtitle = paste(
      "The row label gives the number of studied species used as the",
      "denominator."
    ),
    caption = paste(
      "Colour equals the cell species count divided by the number of studied",
      "species shown in the row label."
    ),
    x = NULL,
    y = NULL
  ) +
  theme_article() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1),
    legend.position = "right"
  )

exports <- rbindlist(list(
  save_editable_figure(
    figure_01, "Figure_04C_01_domain_architecture", 13.0, 7.0
  ),
  save_editable_figure(
    figure_02, "Figure_04C_02_organ_sampling_and_richness", 13.0, 6.8
  ),
  save_editable_figure(
    figure_03, "Figure_04C_03_organ_class_incidence", 14.2, 6.2
  ),
  save_editable_figure(
    figure_04, "Figure_04C_04_class_breadth_vs_richness", 14.0, 9.5
  ),
  save_editable_figure(
    figure_05, "Figure_04C_05_method_class_association", 10.0, 8.2
  ),
  save_editable_figure(
    figure_06, "Figure_04C_06_genus_class_incidence", 14.2, 8.2
  )
), use.names = TRUE)

if (nrow(exports) != 6L ||
    any(!file.exists(c(exports$PDF, exports$SVG, exports$PNG))) ||
    any(exports$PDF_bytes <= 0 | exports$SVG_bytes <= 0 | exports$PNG_bytes <= 0)) {
  stop("A exportacao das figuras nao passou pelo controle de integridade.", call. = FALSE)
}

data.table::fwrite(
  exports,
  file.path(TABLES_DIR, "Table_QA_04C_figure_exports_v1_1_5.csv"),
  sep = ";",
  bom = TRUE,
  quote = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(LOGS_DIR, "04C_figure_export_v1_1_5_sessionInfo.txt")
)

message("ETAPA 04C — FIGURAS REVISADAS v1.1.5 CONCLUIDA")
message("PDF vetorial, SVG editavel e PNG 600 dpi salvos em:")
message(FIGURES_DIR)
