#!/bin/bash
# =============================================================================
# SUB02 — Fase 3: TSS Methylation Profile ±10 kb (adaptado do SUB01)
# Comparações: C1 (Nelore vs HF, ZENODO), C2 (NelFev vs AngFev), C3 (NelJun vs AngJun)
# =============================================================================
set -uo pipefail

SCRIPT="/home/rstudio/projeto/scripts/pipeline/tss_methylation_profile.py"
JSON="/home/rstudio/projeto/subprojeto/SUB02/SUB02_target_genes_coordinates.json"
OUT_BASE="/home/rstudio/projeto/results/SUB02_runs"
GSE="/data/epigenomas/NCBI_GSE161113"
ZEN="/data/epigenomas/ZENODO19555520"

# ─── C1: Nelore vs Holstein-Friesian (sangue lisado, ZENODO) ──────────────────
echo "[$(date)] === C1: Nelore vs HF (ZENODO LysedBlood) ==="
mkdir -p "$OUT_BASE/C1_zenodo_NEL_vs_HF/fase3_tss_10kb"
python3 "$SCRIPT" \
  --methylation \
    "$ZEN/FILTERED_Nelore_ERS4410324_LysedBlood.bismark.cov.gz" \
    "$ZEN/FILTERED_Nelore_ERS4410325_LysedBlood.bismark.cov.gz" \
    "$ZEN/FILTERED_HF_ERS5889248_LysedBlood.bismark.cov.gz" \
    "$ZEN/FILTERED_HF_ERS5889249_LysedBlood.bismark.cov.gz" \
    "$ZEN/FILTERED_HF_ERS5889250_LysedBlood.bismark.cov.gz" \
  --breeds Nelore,Nelore,HF,HF,HF \
  --json "$JSON" \
  --window 10000 \
  --bins 40 \
  --out "$OUT_BASE/C1_zenodo_NEL_vs_HF/fase3_tss_10kb" \
  2>&1 | tee "$OUT_BASE/C1_zenodo_NEL_vs_HF/logs/fase3_tss_10kb.log" \
  && echo "[OK] C1 Fase 3 concluída" \
  || echo "[ERR] C1 Fase 3 falhou"

# ─── C2: Nelore-Fev vs Angus-Fev (heat stress, GSE161113) ────────────────────
echo "[$(date)] === C2: NelFev vs AngFev HS (GSE161113) ==="
mkdir -p "$OUT_BASE/C2_GSE_NelFev_vs_AngFev_HS/fase3_tss_10kb"
python3 "$SCRIPT" \
  --methylation \
    "$GSE/FILTERED_GSM4889769_CpG_context_N01-Feb.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889771_CpG_context_N02-Feb.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889773_CpG_context_N03-Feb.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889775_CpG_context_N04-Feb.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889777_CpG_context_N05-Feb.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889759_CpG_context_A01-Feb.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889761_CpG_context_A02-Feb.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889763_CpG_context_A03-Feb.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889765_CpG_context_A04-Feb.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889767_CpG_context_A05-Feb.bismark.cov.gz" \
  --breeds Nelore,Nelore,Nelore,Nelore,Nelore,Angus,Angus,Angus,Angus,Angus \
  --json "$JSON" \
  --window 10000 \
  --bins 40 \
  --out "$OUT_BASE/C2_GSE_NelFev_vs_AngFev_HS/fase3_tss_10kb" \
  2>&1 | tee "$OUT_BASE/C2_GSE_NelFev_vs_AngFev_HS/logs/fase3_tss_10kb.log" \
  && echo "[OK] C2 Fase 3 concluída" \
  || echo "[ERR] C2 Fase 3 falhou"

# ─── C3: Nelore-Jun vs Angus-Jun (termoneutral, GSE161113) ───────────────────
echo "[$(date)] === C3: NelJun vs AngJun (GSE161113) ==="
mkdir -p "$OUT_BASE/C3_GSE_NelJun_vs_AngJun/fase3_tss_10kb"
python3 "$SCRIPT" \
  --methylation \
    "$GSE/FILTERED_GSM4889770_CpG_context_N01-Jun.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889772_CpG_context_N02-Jun.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889774_CpG_context_N03-Jun.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889776_CpG_context_N04-Jun.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889778_CpG_context_N05-Jun.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889760_CpG_context_A01-Jun.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889762_CpG_context_A02-Jun.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889764_CpG_context_A03-Jun.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889766_CpG_context_A04-Jun.bismark.cov.gz" \
    "$GSE/FILTERED_GSM4889768_CpG_context_A05-Jun.bismark.cov.gz" \
  --breeds Nelore,Nelore,Nelore,Nelore,Nelore,Angus,Angus,Angus,Angus,Angus \
  --json "$JSON" \
  --window 10000 \
  --bins 40 \
  --out "$OUT_BASE/C3_GSE_NelJun_vs_AngJun/fase3_tss_10kb" \
  2>&1 | tee "$OUT_BASE/C3_GSE_NelJun_vs_AngJun/logs/fase3_tss_10kb.log" \
  && echo "[OK] C3 Fase 3 concluída" \
  || echo "[ERR] C3 Fase 3 falhou"

echo "[FIM] SUB02 Fase 3 concluída. Outputs em $OUT_BASE/<COMP>/fase3_tss_10kb/"
