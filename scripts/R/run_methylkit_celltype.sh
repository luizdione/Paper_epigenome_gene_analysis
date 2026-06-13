#!/usr/bin/env bash
# =============================================================================
# SUB02 — Fase 1 methylKit nos contrastes resolvidos por TIPO CELULAR (C8-C14)
# =============================================================================
# Nelore (Bos indicus, treat=0, n=3) vs Holstein-Friesian (Bos taurus, treat=1, n=3)
# em cada tipo celular imune (CD4, CD8, Bcell, gdT, gran, mono, NK).
# Dados: ZENODO19555520, arquivos FILTERED_*.cov.gz (janela +/-15 kb, 73 genes).
# Mesmos parametros do restante do SUB02: mincov=5, qvalue=0.05, diff 15 e 25.
#
# Uso (dentro do container epigenoma_env):
#   bash /home/rstudio/projeto/scripts/R/run_methylkit_celltype.sh
# =============================================================================
set -uo pipefail
PROJ_DIR="${PROJ_DIR:-/home/rstudio/projeto}"
DATA_DIR="${DATA_DIR:-/data/epigenomas/ZENODO19555520}"
SUB02_DIR="$PROJ_DIR/subprojeto/SUB02"
OUT_BASE="$PROJ_DIR/results/SUB02_runs"
SCRIPT_R="$PROJ_DIR/scripts/R/run_methylkit_diff.R"
ASSEMBLY="ARS-UCD1.3"
QVAL=0.05

COMPS=(
  "C8_zenodo_NEL_vs_HF_CD4"
  "C9_zenodo_NEL_vs_HF_CD8"
  "C10_zenodo_NEL_vs_HF_Bcell"
  "C11_zenodo_NEL_vs_HF_gdT"
  "C12_zenodo_NEL_vs_HF_gran"
  "C13_zenodo_NEL_vs_HF_mono"
  "C14_zenodo_NEL_vs_HF_NK"
)

for COMP in "${COMPS[@]}"; do
  META="$SUB02_DIR/samples/samples_${COMP}.csv"
  if [[ ! -f "$META" ]]; then
    echo "[SKIP] $COMP — sample sheet ausente: $META"; continue
  fi
  for DIFF in 15 25; do
    OUT="$OUT_BASE/$COMP/fase1_methylkit_diff${DIFF}"
    mkdir -p "$OUT" "$OUT_BASE/$COMP/logs"
    LOG="$OUT_BASE/$COMP/logs/fase1_diff${DIFF}.log"
    echo "[RUN] $COMP  diff=${DIFF}"
    Rscript "$SCRIPT_R" \
      --metadata "$META" \
      --data-dir "$DATA_DIR" \
      --assembly "$ASSEMBLY" \
      --diff "$DIFF" --qvalue "$QVAL" --min-cov 5 \
      --out "$OUT" 2>&1 | tee "$LOG" || echo "[ERR] $COMP diff=$DIFF falhou — ver $LOG"
  done
done
echo "[FIM] Fase 1 por tipo celular concluida. Saidas em $OUT_BASE/"
