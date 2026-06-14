#!/bin/bash
# Fase 3 — TSS Methylation Profile com janela ±10 kb
# Apenas amostras de Junho (run04 canônico: Angus=treat=1, Nelore=treat=0)
# Executar dentro do container: docker exec epigenoma_env bash ./scripts/pipeline/run_fase3_10kb.sh

set -e

DATA="/data/epigenomas/NCBI_GSE161113"
JSON="/data/genes_alvo/all_target_genes_coordinates.json"
SCRIPT="./scripts/pipeline/tss_methylation_profile.py"
OUT="./results/SUB01_GSE161113_run04/fase3_deeptools_10kb"
LOG="./results/SUB01_GSE161113_run04/logs/fase3_10kb.log"

mkdir -p "$OUT"
mkdir -p "$(dirname $LOG)"

echo "[$(date)] Iniciando Fase 3 TSS ±10 kb" | tee "$LOG"

python3 "$SCRIPT" \
  --methylation \
    "$DATA/FILTERED_GSM4889760_CpG_context_A01-Jun.bismark.cov.gz" \
    "$DATA/FILTERED_GSM4889762_CpG_context_A02-Jun.bismark.cov.gz" \
    "$DATA/FILTERED_GSM4889764_CpG_context_A03-Jun.bismark.cov.gz" \
    "$DATA/FILTERED_GSM4889766_CpG_context_A04-Jun.bismark.cov.gz" \
    "$DATA/FILTERED_GSM4889768_CpG_context_A05-Jun.bismark.cov.gz" \
    "$DATA/FILTERED_GSM4889770_CpG_context_N01-Jun.bismark.cov.gz" \
    "$DATA/FILTERED_GSM4889772_CpG_context_N02-Jun.bismark.cov.gz" \
    "$DATA/FILTERED_GSM4889774_CpG_context_N03-Jun.bismark.cov.gz" \
    "$DATA/FILTERED_GSM4889776_CpG_context_N04-Jun.bismark.cov.gz" \
    "$DATA/FILTERED_GSM4889778_CpG_context_N05-Jun.bismark.cov.gz" \
  --breeds Angus,Angus,Angus,Angus,Angus,Nelore,Nelore,Nelore,Nelore,Nelore \
  --json "$JSON" \
  --window 10000 \
  --bins 40 \
  --out "$OUT" \
  2>&1 | tee -a "$LOG"

echo "[$(date)] Fase 3 concluída. Outputs em: $OUT" | tee -a "$LOG"
ls -lh "$OUT"
