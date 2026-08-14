#!/usr/bin/env bash
set -Eeuo pipefail

# Execute este arquivo no Git Bash. O primeiro argumento, se fornecido,
# substitui a pasta-fonte padrão do projeto organizado.
SOURCE_PROJECT="${1:-/g/Meu Drive/ARTIGOS/2026/Arecaceae}"
OWNER="carloscarollo-UFMS"
REPOSITORY="arecaceae-chemical-research-landscape"
TAG="v1.0.0"
REMOTE_URL="https://github.com/${OWNER}/${REPOSITORY}.git"
REPOSITORY_URL="https://github.com/${OWNER}/${REPOSITORY}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

fail() {
  printf '\nERRO: %s\n' "$1" >&2
  exit 1
}

printf 'Pasta-fonte: %s\n' "$SOURCE_PROJECT"
printf 'Pacote GitHub: %s\n' "$REPO_DIR"

[[ -d "$SOURCE_PROJECT" ]] || fail "A pasta-fonte não foi encontrada. Informe-a como primeiro argumento."
[[ "$SOURCE_PROJECT" != "$REPO_DIR" ]] || fail "Extraia o pacote em uma subpasta; a fonte e o repositório não podem ser a mesma pasta."

# O script analítico arquivado contém 15 entradas nomeadas. A pendência
# editorial mencionava 16, mas não existe uma 16ª entrada no objeto input_paths.
declare -a REQUIRED_INPUTS=(
  "data/derived/study_design_summary.csv"
  "data/interim/04C_incidence_article_taxon_organ_entity.csv"
  "data/interim/04C_articles_250.csv"
  "data/interim/02_species_analysis_with_GBIF_GHS.csv"
  "data/reference/01_species_tree_WCVP_pruned.tre"
  "data/reference/Species_tree_all_genes_rooted_original.tre"
  "data/derived/Table_04C_02_domain_architecture.csv"
  "data/derived/Table_04C_03_class_architecture.csv"
  "data/derived/Table_04C_04_organ_architecture.csv"
  "data/derived/Table_04C_05_organ_class_incidence.csv"
  "data/derived/Table_04C_07_genus_class_incidence.csv"
  "data/derived/Table_05B_04_genus_coverage.csv"
  "data/derived/Table_05B_05_main_selection_models.csv"
  "data/derived/Table_05B_07_trait_specific_models.csv"
  "data/derived/Table_05B_08_socioeconomic_models.csv"
)

for relative_path in "${REQUIRED_INPUTS[@]}"; do
  [[ -s "$SOURCE_PROJECT/$relative_path" ]] || fail "Entrada ausente ou vazia: $SOURCE_PROJECT/$relative_path"
  mkdir -p "$REPO_DIR/$(dirname "$relative_path")"
  cp -f "$SOURCE_PROJECT/$relative_path" "$REPO_DIR/$relative_path"
done

# Copia o script realmente usado na execução arquivada.
FINAL_SCRIPT_SOURCE="$SOURCE_PROJECT/scripts/Arecaceae_single_pipeline_end.R"
[[ -s "$FINAL_SCRIPT_SOURCE" ]] || fail "Script final não encontrado: $FINAL_SCRIPT_SOURCE"
mkdir -p "$REPO_DIR/scripts"
cp -f "$FINAL_SCRIPT_SOURCE" "$REPO_DIR/scripts/Arecaceae_single_pipeline_end.R"

EXPECTED_FINAL_SCRIPT_MD5="00bfb4ae76ada6383e04ed0493368548"
OBSERVED_FINAL_SCRIPT_MD5="$(md5sum "$REPO_DIR/scripts/Arecaceae_single_pipeline_end.R" | awk '{print $1}')"
[[ "$OBSERVED_FINAL_SCRIPT_MD5" == "$EXPECTED_FINAL_SCRIPT_MD5" ]] || fail "O script final mudou. MD5 observado: $OBSERVED_FINAL_SCRIPT_MD5"

EXPECTED_STAGE05B_SCRIPT_MD5="9f38d16e009b848dcad793ce250d8ed0"
OBSERVED_STAGE05B_SCRIPT_MD5="$(md5sum "$REPO_DIR/scripts/upstream/06_model_research_selection.R" | awk '{print $1}')"
[[ "$OBSERVED_STAGE05B_SCRIPT_MD5" == "$EXPECTED_STAGE05B_SCRIPT_MD5" ]] || fail "O script Stage 05B mudou. MD5 observado: $OBSERVED_STAGE05B_SCRIPT_MD5"

# O arquivo abaixo foi produzido pela execução final de 14/08/2026, apesar do
# nome histórico mencionar apenas geração de figuras.
FINAL_SESSION_SOURCE="$SOURCE_PROJECT/results/logs/figure_generation_sessionInfo.txt"
[[ -s "$FINAL_SESSION_SOURCE" ]] || fail "SessionInfo da execução final não encontrado: $FINAL_SESSION_SOURCE"
mkdir -p "$REPO_DIR/results/logs"
cp -f "$FINAL_SESSION_SOURCE" "$REPO_DIR/results/logs/final_pipeline_sessionInfo.txt"

declare -a FINAL_MANIFESTS=(
  "figure_export_manifest.csv"
  "required_input_manifest.csv"
  "Table_iTOL_export_manifest.csv"
  "Table_pipeline_QA.csv"
)
mkdir -p "$REPO_DIR/results/manifests"
for filename in "${FINAL_MANIFESTS[@]}"; do
  source_file="$SOURCE_PROJECT/results/tables/$filename"
  [[ -s "$source_file" ]] || fail "Manifest final ausente ou vazio: $source_file"
  cp -f "$source_file" "$REPO_DIR/results/manifests/$filename"
done

# Retém apenas a matriz normalizada e os diagnósticos necessários para auditar
# os modelos de seleção. As figuras Stage 05B são regeneráveis e não entram.
STAGE05B_MATRIX_SOURCE="$SOURCE_PROJECT/data/interim/05B_species_selection_analysis.csv"
[[ -s "$STAGE05B_MATRIX_SOURCE" ]] || fail "Matriz normalizada Stage 05B ausente: $STAGE05B_MATRIX_SOURCE"
mkdir -p "$REPO_DIR/diagnostics/stage05B/data_processed" "$REPO_DIR/diagnostics/stage05B/tables"
cp -f "$STAGE05B_MATRIX_SOURCE" "$REPO_DIR/diagnostics/stage05B/data_processed/05B_species_selection_analysis.csv"

declare -a STAGE05B_DIAGNOSTICS=(
  "Table_05B_01_functional_size_axes.csv"
  "Table_05B_06_main_model_fit.csv"
  "Table_05B_09_socioeconomic_2x2_tests.csv"
  "Table_05B_10_range_sensitivity.csv"
  "Table_05B_11_continent_model_screen.csv"
  "Table_05B_12_predictor_correlations.csv"
  "Table_05B_13_input_provenance.csv"
  "Table_05B_14_taxonomic_subfamily_fallback.csv"
  "Table_QA_05B.csv"
  "Table_QA_05B_model_convergence.csv"
)
for filename in "${STAGE05B_DIAGNOSTICS[@]}"; do
  source_file="$SOURCE_PROJECT/data/derived/$filename"
  [[ -s "$source_file" ]] || fail "Diagnóstico Stage 05B ausente ou vazio: $source_file"
  cp -f "$source_file" "$REPO_DIR/diagnostics/stage05B/tables/$filename"
done

[[ -s "$REPO_DIR/diagnostics/stage05B/models/05B_selection_bias_models.rds" ]] || fail "Objeto dos modelos Stage 05B ausente no pacote."
[[ -s "$REPO_DIR/diagnostics/stage05B/logs/05B_sessionInfo.txt" ]] || fail "SessionInfo Stage 05B ausente no pacote."
[[ -s "$REPO_DIR/diagnostics/stage05B/logs/05B_execution_log.txt" ]] || fail "Log Stage 05B ausente no pacote."

# Confere se o R usado agora coincide com o R registrado na execução final.
RSCRIPT_BIN="$(command -v Rscript || true)"
if [[ -z "$RSCRIPT_BIN" ]]; then
  for candidate in /c/Program\ Files/R/R-*/bin/Rscript.exe; do
    [[ -x "$candidate" ]] && RSCRIPT_BIN="$candidate"
  done
fi
[[ -n "$RSCRIPT_BIN" ]] || fail "Rscript não foi localizado. Adicione o R ao PATH do Windows."

ARCHIVED_R_VERSION="$(sed -n '1s/^R version \([^ ]*\).*/\1/p' "$REPO_DIR/results/logs/final_pipeline_sessionInfo.txt")"
CURRENT_R_VERSION="$("$RSCRIPT_BIN" --vanilla -e 'cat(as.character(getRversion()))')"
[[ -n "$ARCHIVED_R_VERSION" ]] || fail "Não foi possível ler a versão do R no sessionInfo final."
[[ "$CURRENT_R_VERSION" == "$ARCHIVED_R_VERSION" ]] || fail "O R atual ($CURRENT_R_VERSION) difere do R da execução final ($ARCHIVED_R_VERSION)."

cd "$REPO_DIR"

# Gera o lockfile no mesmo R da execução arquivada.
"$RSCRIPT_BIN" --vanilla -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org"); renv::snapshot(project = ".", prompt = FALSE)'
[[ -s "$REPO_DIR/renv.lock" ]] || fail "O renv.lock não foi criado."

# Manifesto criptográfico de todo o depósito, excluindo apenas metadados Git.
mkdir -p results/manifests
find . -path './.git' -prune -o -type f ! -path './results/manifests/SHA256SUMS.txt' -print0 \
  | sort -z \
  | xargs -0 sha256sum > results/manifests/SHA256SUMS.txt

if [[ ! -d .git ]]; then
  git init -b main
fi
git config user.name "Carlos Alexandre Carollo"
git config user.email "carloscarollo-UFMS@users.noreply.github.com"

git add --all
if ! git diff --cached --quiet; then
  git commit -m "Archive reproducible Arecaceae analysis v1.0.0"
fi

if git rev-parse "$TAG" >/dev/null 2>&1; then
  [[ "$(git rev-list -n 1 "$TAG")" == "$(git rev-parse HEAD)" ]] || fail "A tag $TAG já existe e aponta para outro commit."
else
  git tag -a "$TAG" -m "Immutable manuscript archive $TAG"
fi

if ! git remote get-url origin >/dev/null 2>&1; then
  if git ls-remote "$REMOTE_URL" >/dev/null 2>&1; then
    git remote add origin "$REMOTE_URL"
  elif command -v gh >/dev/null 2>&1; then
    gh auth status
    gh repo create "$OWNER/$REPOSITORY" --public --description "Reproducible data and code for the Arecaceae phytochemistry evidence map" --source=. --remote=origin
  else
    printf '\nO pacote, o commit e a tag %s estão prontos.\n' "$TAG"
    printf 'Crie um repositório público e vazio com o nome %s em:\nhttps://github.com/new\n' "$REPOSITORY"
    printf 'Depois execute novamente este mesmo arquivo.\n'
    exit 2
  fi
fi

[[ "$(git remote get-url origin)" == "$REMOTE_URL" ]] || fail "O remote origin não corresponde a $REMOTE_URL"
git push -u origin main
git push origin "$TAG"

printf '\nPUBLICAÇÃO CONCLUÍDA\n'
printf 'URL estável: %s\n' "$REPOSITORY_URL"
printf 'Versão imutável: %s/tree/%s\n' "$REPOSITORY_URL" "$TAG"
printf 'R arquivado: %s\n' "$ARCHIVED_R_VERSION"
