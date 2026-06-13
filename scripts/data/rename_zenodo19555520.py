#!/usr/bin/env python3
"""
rename_zenodo19555520.py
------------------------
Renomeia os arquivos Bismark coverage do dataset ZENODO19555520
de IDs internos do dataset para o padrão raça + accession ENA.

Padrão ATUAL  : <ID_INTERNO>_<CELULA>_RRBS_trimmed_bismark_bt2_SORTED_CowChrOnly.bismark.cov.gz
Padrão DESTINO: <RACA>_<ENA_SAMPLE_ACCESSION>_<CELULA>.bismark.cov.gz

Exemplo:
  HF3457_Bcell_RRBS_trimmed_bismark_bt2_SORTED_CowChrOnly.bismark.cov.gz
  →  HF_ERS4410317_Bcell.bismark.cov.gz

Fonte do mapeamento: file3_tables10.csv (Additional File 3 — Tabela S10)

Uso:
  python rename_zenodo19555520.py [--data-dir CAMINHO] [--dry-run]

  --data-dir   Pasta com os arquivos .gz  (padrão: diretório do script)
  --dry-run    Mostra o que seria feito sem renomear nada
"""

import os
import re
import argparse

# ---------------------------------------------------------------------------
# Mapeamento: ID interno do dataset → (grupo_raca, ENA_sample_accession)
# Origem: file3_tables10.csv — Table S10 (Powell et al. 2023, Genome Biology; dataset público)
# ---------------------------------------------------------------------------
SAMPLE_ID_TO_ENA = {
    "HF3457":  ("HF", "ERS4410317"),
    "HF3458":  ("HF", "ERS4410318"),
    "HF3471":  ("HF", "ERS4410319"),
    "HF3649":  ("HF", "ERS5889248"),
    "HF3666":  ("HF", "ERS5889249"),
    "HF3683":  ("HF", "ERS5889250"),
    "HF3465":  ("HF", "ERS5889247"),
    "ND230":   ("ND", "ERS4410320"),
    "ND231":   ("ND", "ERS4410321"),
    "ND234":   ("ND", "ERS4410322"),
    "Nelore1": ("Nelore", "ERS4410323"),
    "Nelore2": ("Nelore", "ERS4410324"),
    "Nelore3": ("Nelore", "ERS4410325"),
}

# Tipos celulares reconhecidos pela plataforma Shiny
CELL_TYPES = {"Bcell", "CD4", "CD8", "gdT", "NK", "mono", "gran", "LysedBlood"}

# Regex para extrair: ID_manuscrito + tipo_celular a partir do nome atual
# Exemplos de nomes capturados:
#   HF3457_Bcell_RRBS_trimmed_bismark_bt2_SORTED_CowChrOnly.bismark.cov.gz
#   ND230_LysedBlood_trimmed_bismark_bt2_SORTED_CowChrOnly.bismark.cov.gz
#   Nelore1_CD4_RRBS_trimmed_bismark_bt2_SORTED_CowChrOnly.bismark.cov.gz
FNAME_RE = re.compile(
    r"^(?P<manuscript_id>HF\d+|ND\d+|Nelore\d+)"
    r"_(?P<cell_type>Bcell|CD4|CD8|gdT|NK|mono|gran|LysedBlood)"
    r"_.*\.bismark\.cov\.gz$",
    re.IGNORECASE
)


def build_new_name(manuscript_id: str, cell_type: str) -> str | None:
    """Retorna o novo nome do arquivo ou None se o ID não estiver no mapeamento."""
    entry = MANUSCRIPT_TO_ENA.get(manuscript_id)
    if entry is None:
        return None
    raca, ena_acc = entry
    return f"{raca}_{ena_acc}_{cell_type}.bismark.cov.gz"


def rename_files(data_dir: str, dry_run: bool = False) -> None:
    print(f"\n{'[DRY-RUN] ' if dry_run else ''}Processando arquivos em: {data_dir}\n")

    files = sorted(f for f in os.listdir(data_dir) if f.endswith(".bismark.cov.gz"))

    if not files:
        print("Nenhum arquivo .bismark.cov.gz encontrado na pasta.")
        return

    renamed = 0
    skipped = 0
    errors  = 0

    for fname in files:
        m = FNAME_RE.match(fname)
        if not m:
            print(f"  ⚠  IGNORADO (padrão não reconhecido): {fname}")
            skipped += 1
            continue

        manuscript_id = m.group("manuscript_id")
        # Preserva capitalização original do tipo celular
        cell_type = next(
            (c for c in CELL_TYPES if c.lower() == m.group("cell_type").lower()),
            m.group("cell_type")
        )

        new_name = build_new_name(manuscript_id, cell_type)
        if new_name is None:
            print(f"  ✗  ERRO — ID '{manuscript_id}' não encontrado no mapeamento: {fname}")
            errors += 1
            continue

        src = os.path.join(data_dir, fname)
        dst = os.path.join(data_dir, new_name)

        if os.path.exists(dst) and src != dst:
            print(f"  ⚠  DESTINO JÁ EXISTE — pulando: {new_name}")
            skipped += 1
            continue

        action = "RENOMEAR" if not dry_run else "SIMULAÇÃO"
        print(f"  {action}: {fname}")
        print(f"       →  {new_name}")

        if not dry_run:
            os.rename(src, dst)

        renamed += 1

    print(f"\n{'=' * 60}")
    print(f"Renomeados : {renamed}")
    print(f"Ignorados  : {skipped}")
    print(f"Erros      : {errors}")
    if dry_run:
        print("(Nenhum arquivo foi alterado — modo dry-run)")


def main():
    parser = argparse.ArgumentParser(
        description="Renomeia arquivos Bismark do ZENODO19555520 para o padrão do app Shiny."
    )
    parser.add_argument(
        "--data-dir",
        default=r"E:\data\epigenomas\ZENODO19555520",
        help="Pasta com os arquivos .bismark.cov.gz (padrão: E:\\data\\epigenomas\\ZENODO19555520)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Apenas mostra o que seria feito, sem renomear"
    )
    args = parser.parse_args()

    if not os.path.isdir(args.data_dir):
        print(f"ERRO: diretório não encontrado: {args.data_dir}")
        return

    rename_files(args.data_dir, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
