import json
import gzip
import os
import platform
from pathlib import Path

# --- CONFIGURAÇÕES ---
# Caminho absoluto para o JSON de genes-alvo (relativo a este script)
_HERE = Path(__file__).resolve().parent
PATH_JSON = str(_HERE / "genes_alvo" / "all_target_genes_coordinates.json")
# Fallback: tenta subprojeto/SUB001/ se genes_alvo não existir
if not Path(PATH_JSON).exists():
    _fallback = _HERE / "subprojeto" / "SUB001" / "SUB001_target_genes_coordinates.json"
    if _fallback.exists():
        PATH_JSON = str(_fallback)
# Fix #2 (2026-04-30): auto-detecção Linux/Windows para compatibilidade Docker
# No Docker (WSL/Linux) o SSD E:/ está montado em /data; no Windows usa o path nativo.
# Fix #3 (2026-05-07): diretório correto é "epigenomas" (com 's'), não "epigenoma"
BASE_DIR = "/data/epigenomas" if platform.system() == "Linux" else r"E:\data\epigenomas"


def carregar_alvos():
    with open(PATH_JSON, 'r', encoding='utf-8') as f:
        config = json.load(f)
    genes = config['genes']
    regioes = {}
    for g, info in genes.items():
        chrom = str(info['chromosome'])
        start, end = info['flanking_start'], info['flanking_end']
        if chrom not in regioes: regioes[chrom] = []
        regioes[chrom].append((start, end))
    return regioes


def processar_pastas():
    regioes = carregar_alvos()

    # Varre as subpastas (ex: NCBI_GSE161113, NCBI_GSE231711)
    for root, dirs, files in os.walk(BASE_DIR):
        for file in files:
            # Filtra apenas arquivos .gz que NAO sao os ja filtrados
            if file.endswith(".gz") and not file.startswith("FILTERED_"):
                path_in = os.path.join(root, file)
                path_out = os.path.join(root, "FILTERED_" + file)

                # Se o filtrado ja existe, pula para economizar tempo
                if os.path.exists(path_out):
                    print(f"Pulando {file}, filtro ja existe.")
                    continue

                print(f"Filtrando: {path_in}")
                try:
                    with gzip.open(path_in, 'rt') as f_in, gzip.open(path_out, 'wt') as f_out:
                        count = 0
                        for line in f_in:
                            parts = line.split()
                            if len(parts) < 2: continue

                            chr_file = parts[0].replace("chr", "")
                            try:
                                pos = int(parts[1])
                            except ValueError:
                                continue

                            if chr_file in regioes:
                                for s, e in regioes[chr_file]:
                                    if s <= pos <= e:
                                        f_out.write(line)
                                        count += 1
                                        break
                        print(f"  -> Sucesso! {count} CpGs salvos em {path_out}")
                except Exception as e:
                    print(f"  -> Erro em {file}: {e}")