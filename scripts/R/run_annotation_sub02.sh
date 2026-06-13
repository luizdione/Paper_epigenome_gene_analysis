#!/usr/bin/env bash
# =============================================================================
# SUB02 — Fase 2: Anotação dos DMCs (overlap TSS / gene body / flanking)
# =============================================================================
# Para cada comparação C1-C5 e cada threshold (diff=15, 25), executa
# run_annotation_glm.R usando o JSON de coordenadas do SUB02 como fallback.
# Saídas em: results/SUB02_runs/<COMP>/fase2_annotation_diff{15,25}/
# =============================================================================

set -uo pipefail
PROJ_DIR="${PROJ_DIR:-/home/rstudio/projeto}"
SUB02_DIR="$PROJ_DIR/subprojeto/SUB02"
OUT_BASE="$PROJ_DIR/results/SUB02_runs"
SCRIPT_R="$PROJ_DIR/scripts/R/run_annotation_glm.R"
JSON="$SUB02_DIR/SUB02_target_genes_coordinates.json"

COMPS=(
  "C1_zenodo_NEL_vs_HF"
  "C2_GSE_NelFev_vs_AngFev_HS"
  "C3_GSE_NelJun_vs_AngJun"
  "C4_GSE_NelJun_vs_NelFev_HS"
  "C5_GSE_AngJun_vs_AngFev_HS"
)

for COMP in "${COMPS[@]}"; do
  META="$SUB02_DIR/samples/samples_${COMP}.csv"
  for DIFF in 15 25; do
    DMRS="$OUT_BASE/$COMP/fase1_methylkit_diff${DIFF}/DMCs_significant.csv"
    OUT="$OUT_BASE/$COMP/fase2_annotation_diff${DIFF}"
    LOG="$OUT_BASE/$COMP/logs/fase2_annotation_diff${DIFF}.log"
    mkdir -p "$OUT" "$(dirname "$LOG")"

    if [[ ! -f "$DMRS" ]]; then
      echo "[SKIP] $COMP diff=$DIFF — DMCs_significant.csv ausente: $DMRS"
      continue
    fi

    N=$(tail -n +2 "$DMRS" | wc -l)
    echo "[RUN] $COMP diff=${DIFF} — $N DMCs para anotar"

    Rscript "$SCRIPT_R" \
      --dmrs     "$DMRS" \
      --metadata "$META" \
      --json     "$JSON" \
      --out      "$OUT" 2>&1 | tee "$LOG" \
      && echo "[OK] $COMP diff=$DIFF -> $OUT" \
      || echo "[ERR] $COMP diff=$DIFF — ver $LOG"
  done
done

echo ""
echo "[FIM] Anotação SUB02 Fase 2 concluída. Resumo:"
for COMP in "${COMPS[@]}"; do
  for DIFF in 15 25; do
    ANN="$OUT_BASE/$COMP/fase2_annotation_diff${DIFF}/annotation_table.csv"
    if [[ -f "$ANN" ]]; then
      N_PROM=$(grep -c "promoter" "$ANN" 2>/dev/null || echo 0)
      N_TOTAL=$(tail -n +2 "$ANN" | wc -l)
      echo "  $COMP diff=$DIFF: $N_TOTAL DMCs anotados ($N_PROM promoter)"
    fi
  done
done
