#!/usr/bin/env python3
"""Debug completo da Fase 3: verifica overlap cromossômico e formato de leitura."""
import gzip, json, sys
from collections import defaultdict

DATA = "/data/epigenomas/NCBI_GSE161113"
JSON = "/data/genes_alvo/all_target_genes_coordinates.json"

# Amostras Jun
SAMPLES = {
    "Angus_A01": f"{DATA}/FILTERED_GSM4889760_CpG_context_A01-Jun.bismark.cov.gz",
    "Angus_A02": f"{DATA}/FILTERED_GSM4889762_CpG_context_A02-Jun.bismark.cov.gz",
    "Nelore_N01": f"{DATA}/FILTERED_GSM4889770_CpG_context_N01-Jun.bismark.cov.gz",
    "Nelore_N02": f"{DATA}/FILTERED_GSM4889772_CpG_context_N02-Jun.bismark.cov.gz",
}

print("=== Formato das colunas (primeiras 3 linhas de cada) ===")
for label, path in SAMPLES.items():
    print(f"\n{label}:")
    with gzip.open(path, 'rt') as fh:
        for i, line in enumerate(fh):
            if i >= 3: break
            cols = line.strip().split('\t')
            print(f"  cols: {cols[:6]} → col[3]={cols[3]}, col[4]={cols[4]}, col[5]={cols[5]}")

print("\n\n=== Cromossomos em cada amostra ===")
for label, path in SAMPLES.items():
    chroms = defaultdict(int)
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if cols: chroms[cols[0]] += 1
    print(f"\n{label}: {dict(sorted(chroms.items()))}")

print("\n\n=== Cromossomos no JSON (genes-alvo) ===")
with open(JSON) as fh:
    data = json.load(fh)
gene_chroms = {}
for name, info in data.get("genes", {}).items():
    gene_chroms[name] = info.get("chromosome","?")
chrom_counts = defaultdict(int)
for chrom in gene_chroms.values():
    chrom_counts[chrom] += 1
print("Cromossomos únicos:", dict(sorted(chrom_counts.items())))
print("Exemplo de genes:", list(gene_chroms.items())[:5])
