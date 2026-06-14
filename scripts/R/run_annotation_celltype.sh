#!/usr/bin/env bash
# =============================================================================
# SUB02 — Fase 2 ANOTAÇÃO (ChIPseeker) para os contrastes por tipo celular C8-C14
# Produz fase2_annotation_diff{15,25}/annotation_table.csv — formato que a
# ponte multi-ômica (bridge_deg_dmr.py) consome.
# =============================================================================
set -uo pipefail
PROJ_DIR="${PROJ_DIR:-.}"
SUB02_DIR="$PROJ_DIR/subprojeto/SUB02"
OUT_BASE="$PROJ_DIR/results/SUB02_runs"
SCRIPT_R="$PROJ_DIR/scripts/R/run_annotation_glm.R"
JSON="$SUB02_DIR/SUB02_target_genes_coordinates.json"

COMPS=(
  "C8_zenodo_NEL_vs_HF_CD4"  "C9_zenodo_NEL_vs_HF_CD8"
  "C10_zenodo_NEL_vs_HF_Bcell" "C11_zenodo_NEL_vs_HF_gdT"
  "C12_zenodo_NEL_vs_HF_gran" "C13_zenodo_NEL_vs_HF_mono"
  "C14_zenodo_NEL_vs_HF_NK"
)

for COMP in "${COMPS[@]}"; do
  META="$SUB02_DIR/samples/samples_${COMP}.csv"
  for DIFF in 15 25; do
    DMRS="$OUT_BASE/$COMP/fase1_methylkit_diff${DIFF}/DMCs_significant.csv"
    OUT="$OUT_BASE/$COMP/fase2_annotation_diff${DIFF}"
    LOG="$OUT_BASE/$COMP/logs/fase2_annotation_diff${DIFF}.log"
    mkdir -p "$OUT" "$(dirname "$LOG")"
    if [[ ! -f "$DMRS" ]]; then echo "[SKIP] $COMP diff=$DIFF — sem DMCs"; continue; fi
    N=$(tail -n +2 "$DMRS" | wc -l)
    echo "[RUN] $COMP diff=${DIFF} — $N DMCs"
    Rscript "$SCRIPT_R" --dmrs "$DMRS" --metadata "$META" --json "$JSON" --out "$OUT" 2>&1 | tail -4 \
      && echo "[OK] $COMP diff=$DIFF" || echo "[ERR] $COMP diff=$DIFF — ver $LOG"
  done
done
echo "[FIM] Anotação por tipo celular concluída."
for COMP in "${COMPS[@]}"; do
  ANN="$OUT_BASE/$COMP/fase2_annotation_diff15/annotation_table.csv"
  if [[ -f "$ANN" ]]; then
    N_PROM=$(grep -c "promoter" "$ANN" 2>/dev/null || echo 0)
    N_TOTAL=$(($(wc -l < "$ANN") - 1))
    echo "  $COMP diff15: $N_TOTAL DMCs anotados ($N_PROM promoter)"
  fi
done
