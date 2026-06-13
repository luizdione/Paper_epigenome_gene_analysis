#!/usr/bin/env python3
"""
SUB02 — Filtragem de arquivos de cobertura Bismark para janelas dos 73 genes-alvo
Diferenças vs SUB01:
  - Janela flanking ±15 kb (era ±10 kb)
  - JSON aponta para SUB02_target_genes_coordinates.json
  - Mantém compatibilidade com layout /data/epigenomas/{NCBI_GSE161113,ZENODO19555520}/

Uso:
  python3 baixar_e_filtrar_metilacao_SUB02.py [--dry-run]

O script varre BASE_DIR, encontra arquivos .gz que NÃO começam com FILTERED_,
recorta para os intervalos flanking_start..flanking_end de cada gene SUB02 e
escreve FILTERED_<file> no mesmo diretório.
"""

import json
import gzip
import os
import platform
import sys
from pathlib import Path

# --- CONFIGURAÇÕES ---
_HERE = Path(__file__).resolve().parent
PATH_JSON = str(_HERE.parent.parent / "subprojeto" / "SUB02" / "SUB02_target_genes_coordinates.json")

# Auto-detecção Linux/Windows (igual SUB01)
BASE_DIR = "/data/epigenomas" if platform.system() == "Linux" else r"E:\data\epigenomas"

DRY_RUN = "--dry-run" in sys.argv


def carregar_alvos():
    with open(PATH_JSON, 'r', encoding='utf-8') as f:
        config = json.load(f)
    genes = config['genes']
    flank_label = config.get('_meta', {}).get('flanking_window_label', '?')
    print(f"[INFO] SUB02 — {len(genes)} genes, janela flanking {flank_label}")
    regioes = {}
    for g, info in genes.items():
        chrom = str(info['chromosome'])
        start, end = info['flanking_start'], info['flanking_end']
        regioes.setdefault(chrom, []).append((start, end, g))
    # Ordenar para busca eficiente
    for c in regioes:
        regioes[c].sort()
    return regioes


def processar_pastas():
    regioes = carregar_alvos()
    if not Path(BASE_DIR).exists():
        print(f"[ERRO] BASE_DIR não existe: {BASE_DIR}")
        print("       Confira se o disco com /data/epigenomas (ou E:\\data\\epigenomas) está montado.")
        sys.exit(1)

    n_files = 0
    n_total_lines = 0
    for root, dirs, files in os.walk(BASE_DIR):
        for file in files:
            if not file.endswith(".gz"):
                continue
            if file.startswith("FILTERED_"):
                continue
            path_in = os.path.join(root, file)
            path_out = os.path.join(root, "FILTERED_" + file)

            if os.path.exists(path_out):
                print(f"  [skip] {file} (FILTERED_ já existe)")
                continue
            n_files += 1

            if DRY_RUN:
                print(f"  [dry-run] processaria {file}")
                continue

            print(f"  [filtrar] {path_in}")
            try:
                with gzip.open(path_in, 'rt') as f_in, gzip.open(path_out, 'wt') as f_out:
                    count = 0
                    for line in f_in:
                        parts = line.split()
                        if len(parts) < 2:
                            continue
                        chr_file = parts[0].replace("chr", "")
                        try:
                            pos = int(parts[1])
                        except ValueError:
                            continue
                        if chr_file in regioes:
                            for s, e, _gene in regioes[chr_file]:
                                if s <= pos <= e:
                                    f_out.write(line)
                                    count += 1
                                    break
                                if pos < s:
                                    break
                    print(f"    -> {count} CpGs salvos")
                    n_total_lines += count
            except Exception as ex:
                print(f"    -> ERRO: {ex}")

    print(f"\n[FIM] Arquivos processados: {n_files}, total CpGs filtrados: {n_total_lines}")


if __name__ == "__main__":
    processar_pastas()
