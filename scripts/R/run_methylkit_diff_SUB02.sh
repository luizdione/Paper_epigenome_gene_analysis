#!/usr/bin/env bash
# =============================================================================
# SUB02 — Pipeline Fase 1 methylKit em loop pelas 5 comparações
# =============================================================================
# Diferenças vs SUB01:
#   - Roda 5 comparações sequenciais (C1..C5)
#   - Resultados em results/SUB02_runs/<comp_id>/fase1_methylkit/
#   - Mantém os mesmos parâmetros padrão (mincov=5, diff=15, qvalue=0.05)
#   - Adiciona também rodada conservadora (diff=25) em paralelo
#
# Pré-requisito:
#   - Arquivos FILTERED_*.bismark.cov.gz já gerados pelo
#     baixar_e_filtrar_metilacao_SUB02.py (janela ±15 kb)
#   - Cada CSV de samples (samples/samples_C*.csv) deve ter o caminho dos arquivos
#     em $DATA_DIR/<dataset>/...
#
# Uso:  bash run_methylkit_diff_SUB02.sh
# =============================================================================

set -euo pipefail
PROJ_DIR="${PROJ_DIR:-.}"
DATA_DIR="${DATA_DIR:-/data/epigenomas}"
SUB02_DIR="$PROJ_DIR/subprojeto/SUB02"
OUT_BASE="$PROJ_DIR/results/SUB02_runs"
SCRIPT_R="$PROJ_DIR/scripts/R/run_methylkit_diff.R"  # mesmo script do SUB01, parametrizado

ASSEMBLY="ARS-UCD1.3"
QVAL=0.05

declare -A DATA_SUBDIR=(
  ["C1_zenodo_NEL_vs_HF"]="ZENODO19555520"
  ["C2_GSE_NelFev_vs_AngFev_HS"]="NCBI_GSE161113"
  ["C3_GSE_NelJun_vs_AngJun"]="NCBI_GSE161113"
  ["C4_GSE_NelJun_vs_NelFev_HS"]="NCBI_GSE161113"
  ["C5_GSE_AngJun_vs_AngFev_HS"]="NCBI_GSE161113"
)

for COMP in C1_zenodo_NEL_vs_HF C2_GSE_NelFev_vs_AngFev_HS C3_GSE_NelJun_vs_AngJun C4_GSE_NelJun_vs_NelFev_HS C5_GSE_AngJun_vs_AngFev_HS; do
  META="$SUB02_DIR/samples/samples_${COMP}.csv"
  DSUB="${DATA_SUBDIR[$COMP]}"
  DSRC="$DATA_DIR/$DSUB"

  if [[ ! -f "$META" ]]; then
    echo "[SKIP] $COMP — sample sheet ausente: $META"
    continue
  fi
  if [[ ! -d "$DSRC" ]]; then
    echo "[SKIP] $COMP — data dir ausente: $DSRC"
    continue
  fi

  for DIFF in 15 25; do
    OUT="$OUT_BASE/$COMP/fase1_methylkit_diff${DIFF}"
    mkdir -p "$OUT"
    LOG="$OUT_BASE/$COMP/logs/fase1_diff${DIFF}.log"
    mkdir -p "$(dirname "$LOG")"
    echo "[RUN] $COMP  diff=${DIFF}"
    echo "  meta=$META  data=$DSRC  out=$OUT"
    Rscript "$SCRIPT_R" \
      --metadata "$META" \
      --data-dir "$DSRC" \
      --assembly "$ASSEMBLY" \
      --diff "$DIFF" --qvalue "$QVAL" --min-cov 5 \
      --out "$OUT" 2>&1 | tee "$LOG" || echo "[ERR] $COMP diff=$DIFF falhou — ver $LOG"
  done
done

echo "[FIM] Pipeline SUB02 Fase 1 concluído. Saídas em $OUT_BASE/"
