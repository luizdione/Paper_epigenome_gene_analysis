#!/usr/bin/env bash
# =============================================================================
# SUB02 — Fase 2 (interseção DMC × genes-alvo) para os contrastes por tipo
# celular C8-C14. Usa os DMCs diff15 de cada contraste e o painel de 73 genes.
# Saída: results/SUB02_runs/<comp>/fase2_intersect/dmrs_in_target_genes.csv
# =============================================================================
set -uo pipefail
PROJ_DIR="${PROJ_DIR:-/home/rstudio/projeto}"
OUT_BASE="$PROJ_DIR/results/SUB02_runs"
JSON="$PROJ_DIR/subprojeto/SUB02/SUB02_target_genes_coordinates.json"
ISECT="$PROJ_DIR/scripts/pipeline/intersect_dmrs_genes.py"

COMPS=(
  "C8_zenodo_NEL_vs_HF_CD4"  "C9_zenodo_NEL_vs_HF_CD8"
  "C10_zenodo_NEL_vs_HF_Bcell" "C11_zenodo_NEL_vs_HF_gdT"
  "C12_zenodo_NEL_vs_HF_gran" "C13_zenodo_NEL_vs_HF_mono"
  "C14_zenodo_NEL_vs_HF_NK"
)

for COMP in "${COMPS[@]}"; do
  DMC="$OUT_BASE/$COMP/fase1_methylkit_diff15/DMCs_significant.csv"
  OUT="$OUT_BASE/$COMP/fase2_intersect/dmrs_in_target_genes.csv"
  mkdir -p "$(dirname "$OUT")"
  if [[ ! -f "$DMC" ]]; then echo "[SKIP] $COMP — sem DMCs: $DMC"; continue; fi
  echo "[F2] $COMP"
  python3 "$ISECT" --dmrs "$DMC" --json "$JSON" --out "$OUT" 2>&1 | tail -3
  n=$(($(wc -l < "$OUT") - 1))
  genes=$(tail -n +2 "$OUT" | awk -F, '{print $NF}' | sort -u | wc -l)
  echo "    -> $n DMCs intersectadas | linhas de genes: ver arquivo"
done
echo "[FIM] Fase 2 por tipo celular concluida."
