#!/bin/bash
# Script de recuperação pós-redownload
# Rodar APÓS baixar_cromossomos.py terminar de baixar Chr4, 16, 22, 25, 29, X
#
# Passos:
#   1. Gerar índices .tbi para os VCFs re-baixados
#   2. Re-rodar rebuild_fasta_missing_genes.py para preencher os 8 genes
#
# Uso (no container Docker):
#   docker exec epigenoma_env bash ./scripts/rebuild_after_redownload.sh

set -e
cd .

echo "=== [1/3] Verificando integridade dos VCFs re-baixados ==="
python3 scripts/check_vcf_integrity_fast.py /data/1000bg_vcf --quick

echo ""
echo "=== [2/3] Gerando índices tabix .tbi ==="
for ch in 4 16 22 25 29 X; do
    vcf="/data/1000bg_vcf/Chr${ch}-Run8-TAUIND-public.vcf.gz"
    tbi="${vcf}.tbi"
    if [ ! -f "$vcf" ]; then
        echo "  [SKIP] $vcf não existe"
        continue
    fi
    if [ -f "$tbi" ] && [ "$tbi" -nt "$vcf" ]; then
        echo "  [OK] $tbi já atualizado"
        continue
    fi
    echo "  [INDEX] tabix -p vcf $vcf (~5-15 min para arquivo grande)..."
    tabix -p vcf "$vcf"
    echo "  [OK] $tbi gerado"
done

echo ""
echo "=== [3/3] Re-extraindo FASTA dos 8 genes faltantes ==="
python3 scripts/rebuild_fasta_missing_genes.py --sub SUB01 --vcf-dir /data/1000bg_vcf --window 500

echo ""
echo "=== Validando resultado final ==="
python3 scripts/fasta_quality_report.py SUB01

echo ""
echo "=== CONCLUÍDO ==="
echo "FASTA atualizado em: subprojeto/SUB01/SUB01_reference_sequences.fa"
echo "Próximo passo: re-rodar pipeline para usar novas sequências (Fases 4 HOMER e 6 MOFA2)"
