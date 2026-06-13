#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gerar_samples_csv.py — Gerador interativo de samples.csv para o pipeline Epigenoma Bovino
===========================================================================================
Varre uma pasta de dados, interpreta automaticamente os nomes dos arquivos
(GSE161113, ZENODO19555520 ou qualquer padrão customizado), exibe uma tabela
editável e deixa o usuário:

    1. Revisar/corrigir raça, espécie, condição e tipo celular por arquivo
    2. Marcar individualmente quais arquivos entram na análise
    3. Definir quais amostras são CONTROLE e quais são TESTE
    4. Exportar o samples.csv pronto para o run_pipeline.py

Padrões reconhecidos automaticamente:
    ┌─────────────────────────────────────────────────────────────────────┐
    │ GSE161113 (NCBI/GEO — Heat stress Angus vs Nelore)                 │
    │   A01_Feb_CpG_filtered.txt.gz  → Angus, Bos taurus, Fevereiro      │
    │   A03_Jun_CpG_filtered.txt.gz  → Angus, Bos taurus, Junho          │
    │   N02_Feb_CpG_filtered.txt.gz  → Nelore, Bos indicus, Fevereiro    │
    │   FILTERED_A01_...             → prefixo FILTERED ignorado          │
    ├─────────────────────────────────────────────────────────────────────┤
    │ ZENODO19555520 (Immune Epigenome HF vs ND)                          │
    │   HF_ERS4410317_CD4.bismark.cov.gz  → Holstein Friesian, CD4       │
    │   ND_ERS4410320_Bcell.bismark.cov.gz → N'Dama, Bcell               │
    ├─────────────────────────────────────────────────────────────────────┤
    │ Genérico (qualquer .gz / .cov.gz / .txt.gz)                         │
    │   → campos preenchidos com "?" para edição manual                   │
    └─────────────────────────────────────────────────────────────────────┘

Uso:
    # Modo interativo (padrão — recomendado)
    python gerar_samples_csv.py --data-dir E:/data/epigenomas/NCBI_GSE161113

    # Modo não-interativo: gera CSV sem perguntas (útil em pipeline)
    python gerar_samples_csv.py --data-dir /data/epigenomas/NCBI_GSE161113 --auto --out samples.csv

    # Apenas imprime a tabela parseada sem salvar
    python gerar_samples_csv.py --data-dir /data/epigenomas/NCBI_GSE161113 --dry-run

Requisitos: Python 3.8+ (sem dependências externas)

Autor: Luiz Dione Barbosa De-Melo — GenOmics Epigenoma 2026
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path
from typing import Optional

# ─────────────────────────────────────────────────────────────────────────────
# Tabelas de conhecimento
# ─────────────────────────────────────────────────────────────────────────────

# GSE161113: prefixo A = Angus (taurus), N = Nelore (indicus)
GSE_PREFIX = {
    "A": {
        "breed":   "Angus",
        "species": "Bos taurus",
        "species_code": "taurus",
    },
    "N": {
        "breed":   "Nelore",
        "species": "Bos indicus",
        "species_code": "indicus",
    },
}

# GSE161113: mês de coleta codificado no nome
GSE_MONTH = {
    "Feb": {"month": "Fevereiro", "season": "verao_sul"},   # Fevereiro = pico do calor
    "Jun": {"month": "Junho",     "season": "inverno_sul"},  # Junho = período frio
}

# ZENODO: código de raça
ZENODO_BREED = {
    "HF": {
        "breed":        "Holstein Friesian",
        "species":      "Bos taurus",
        "species_code": "taurus",
    },
    "ND": {
        "breed":        "N'Dama",
        "species":      "Bos taurus africanus",
        "species_code": "taurus",
    },
    "Nelore": {
        "breed":        "Nelore",
        "species":      "Bos indicus",
        "species_code": "indicus",
    },
}

# Tipos celulares reconhecidos
CELL_TYPES = {
    "Bcell": "Célula B",
    "CD4":   "CD4 T cell",
    "CD8":   "CD8 T cell",
    "gdT":   "γδ T cell",
    "NK":    "NK cell",
    "mono":  "Monócito",
    "gran":  "Granulócito",
    "LysedBlood": "Sangue lisado (bulk)",
}

# Extensões de arquivo aceitas
ACCEPTED_EXT = re.compile(
    r"\.(gz|cov\.gz|txt\.gz|bismark\.cov\.gz|bed\.gz|tsv\.gz)$",
    re.IGNORECASE
)


# ─────────────────────────────────────────────────────────────────────────────
# Parsers de nome de arquivo
# ─────────────────────────────────────────────────────────────────────────────

def parse_gse161113(filename: str) -> Optional[dict]:
    """
    Padrão GSE161113 — detecta A01–A05 (Angus) e N01–N05 (Nelore)
    com mês Feb (fevereiro/calor) ou Jun (junho/controle).

    Aceita qualquer variação de prefixo:
      FILTERED_A01_Feb_CpG.txt.gz
      FILTERED_N03_Jun_bismark.cov.gz
      A01_Feb_CpG_report.txt.gz
      GSM4905205_A01_Feb_CpG_report.txt.gz   ← prefixo GEO
      SRR12345_N02_Jun.cov.gz                ← prefixo SRA
      FILTERED_A01_bismark.cov.gz            ← sem mês
    """
    # Busca o padrão [AN][01-05] em qualquer posição do nome (search, não match)
    # (?<![A-Za-z]) evita capturar letras no meio de palavras
    m = re.search(
        r"(?<![A-Za-z])([AN])(\d{1,2})(?:_(Feb|Jun))?(?=[_\.\-]|$)",
        filename, re.IGNORECASE
    )
    if not m:
        return None

    prefix = m.group(1).upper()
    num    = m.group(2).zfill(2)   # garante zero à esquerda: 1 → "01"
    month  = (m.group(3) or "").capitalize()   # "Feb", "Jun" ou ""

    if prefix not in GSE_PREFIX:
        return None

    prefix = m.group(1).upper()
    num    = m.group(2)
    month  = (m.group(3) or "").capitalize()  # "Feb" ou "Jun" ou ""

    if prefix not in GSE_PREFIX:
        return None

    info = GSE_PREFIX[prefix].copy()
    month_info = GSE_MONTH.get(month, {})

    sample_id = f"{info['breed'].lower()}_{prefix}{num}"
    if month:
        sample_id += f"_{month.lower()}"

    # Heurística de tratamento: Junho = controle (frio), Fevereiro = estresse calor
    if month == "Feb":
        treatment = "heat_stress"
    elif month == "Jun":
        treatment = "control"
    else:
        treatment = "?"

    return {
        "sampleID":    sample_id,
        "breed":       info["breed"],
        "species":     info["species"],
        "species_code": info["species_code"],
        "tissue_cell": "sangue_total",   # GSE161113 = sangue completo
        "month":       month_info.get("month", month or "?"),
        "treatment":   treatment,
        "file":        filename,
        "dataset":     "GSE161113",
        "_parser":     "gse161113",
    }


def parse_zenodo(filename: str) -> Optional[dict]:
    """
    Padrão ZENODO19555520:
      <RACA>_<ENA_ID>_<TIPO_CELULA>.<ext>
      Ex: HF_ERS4410317_CD4.bismark.cov.gz
          ND_ERS4410320_Bcell.bismark.cov.gz
          Nelore_ERS4410323_CD4.bismark.cov.gz
    """
    base = re.sub(r"\.(bismark\.cov|cov|txt|bed|tsv)?\.gz$", "", filename, flags=re.IGNORECASE)
    m = re.match(
        r"^(HF|ND|Nelore)_([A-Za-z0-9]+)_(Bcell|CD4|CD8|gdT|NK|mono|gran|LysedBlood)",
        base, re.IGNORECASE
    )
    if not m:
        return None

    breed_code = m.group(1)
    # normaliza capitalização para chave do dict
    breed_key = next((k for k in ZENODO_BREED if k.lower() == breed_code.lower()), None)
    if not breed_key:
        return None

    ena_id    = m.group(2)
    cell_code = m.group(3)
    cell_full = CELL_TYPES.get(cell_code, cell_code)
    info      = ZENODO_BREED[breed_key]

    sample_id = f"{breed_key}_{ena_id}_{cell_code}"

    return {
        "sampleID":    sample_id,
        "breed":       info["breed"],
        "species":     info["species"],
        "species_code": info["species_code"],
        "tissue_cell": cell_full,
        "month":       "N/A",
        "treatment":   "?",   # usuário define qual grupo é controle
        "file":        filename,
        "dataset":     "ZENODO19555520",
        "_parser":     "zenodo",
        "_ena_id":     ena_id,
        "_cell_code":  cell_code,
    }


def parse_generic(filename: str) -> dict:
    """Fallback: extrai o que puder do nome do arquivo."""
    stem = re.sub(r"(\.(bismark\.cov|cov|txt|bed|tsv))?\.gz$", "", filename, re.IGNORECASE)
    stem = re.sub(r"^FILTERED_", "", stem, flags=re.IGNORECASE)

    # Tenta detectar pistas no nome
    breed = "?"
    species = "?"
    species_code = "?"
    cell = "?"
    treatment = "?"

    stem_lower = stem.lower()
    if "holstein" in stem_lower or "hf" in stem_lower:
        breed, species, species_code = "Holstein Friesian", "Bos taurus", "taurus"
    elif "angus" in stem_lower:
        breed, species, species_code = "Angus", "Bos taurus", "taurus"
    elif "nelore" in stem_lower:
        breed, species, species_code = "Nelore", "Bos indicus", "indicus"
    elif "gir" in stem_lower:
        breed, species, species_code = "Gir", "Bos indicus", "indicus"
    elif "girolando" in stem_lower:
        breed, species, species_code = "Girolando", "Híbrido", "hybrid"
    elif "nd" in stem_lower or "ndama" in stem_lower:
        breed, species, species_code = "N'Dama", "Bos taurus africanus", "taurus"

    for ct in CELL_TYPES:
        if ct.lower() in stem_lower:
            cell = CELL_TYPES[ct]
            break

    if "control" in stem_lower or "ctrl" in stem_lower or "jun" in stem_lower:
        treatment = "control"
    elif "heat" in stem_lower or "stress" in stem_lower or "feb" in stem_lower:
        treatment = "heat_stress"

    return {
        "sampleID":    stem,
        "breed":       breed,
        "species":     species,
        "species_code": species_code,
        "tissue_cell": cell,
        "month":       "?",
        "treatment":   treatment,
        "file":        filename,
        "dataset":     "custom",
        "_parser":     "generic",
    }


def parse_filename(filename: str) -> dict:
    """Tenta parsers em ordem de especificidade."""
    return (
        parse_gse161113(filename)
        or parse_zenodo(filename)
        or parse_generic(filename)
    )


# ─────────────────────────────────────────────────────────────────────────────
# Varredura de pasta
# ─────────────────────────────────────────────────────────────────────────────

def scan_directory(data_dir: Path) -> list[dict]:
    """Varre data_dir recursivamente e parseia todos os arquivos aceitos."""
    if not data_dir.exists():
        print(f"\n✗ Pasta não encontrada: {data_dir}")
        sys.exit(1)

    files = sorted(
        f for f in data_dir.rglob("*")
        if f.is_file() and ACCEPTED_EXT.search(f.name)
    )

    if not files:
        print(f"\n✗ Nenhum arquivo .gz encontrado em: {data_dir}")
        sys.exit(1)

    results = []
    for f in files:
        # Caminho relativo ao data_dir para a coluna 'file'
        rel = str(f.relative_to(data_dir))
        row = parse_filename(f.name)
        row["file"] = rel
        results.append(row)

    return results


# ─────────────────────────────────────────────────────────────────────────────
# Exibição de tabela no terminal
# ─────────────────────────────────────────────────────────────────────────────

COL_W = {
    "idx":         4,
    "sampleID":   22,
    "breed":      16,
    "species_code": 8,
    "tissue_cell": 14,
    "treatment":  12,
    "group":       7,
    "file":       40,
}

def trunc(s: str, n: int) -> str:
    return (s[:n-1] + "…") if len(s) > n else s.ljust(n)


def print_table(rows: list[dict], selected: list[bool], groups: list[str]) -> None:
    """Imprime tabela formatada com status de seleção e grupo."""
    header = (
        f"{'#':>{COL_W['idx']}}  "
        f"{'✓':<2}  "
        f"{'sampleID':<{COL_W['sampleID']}}  "
        f"{'breed':<{COL_W['breed']}}  "
        f"{'sp':<{COL_W['species_code']}}  "
        f"{'tissue/cell':<{COL_W['tissue_cell']}}  "
        f"{'treatment':<{COL_W['treatment']}}  "
        f"{'group':<{COL_W['group']}}  "
        f"file"
    )
    sep = "─" * len(header)
    print()
    print(sep)
    print(header)
    print(sep)
    for i, (row, sel, grp) in enumerate(zip(rows, selected, groups), 1):
        mark = "✓" if sel else " "
        grp_disp = {"C": "CTRL", "T": "TESTE", "": "─"}.get(grp, grp)
        print(
            f"{i:>{COL_W['idx']}}  "
            f"{mark:<2}  "
            f"{trunc(row['sampleID'], COL_W['sampleID']):<{COL_W['sampleID']}}  "
            f"{trunc(row['breed'], COL_W['breed']):<{COL_W['breed']}}  "
            f"{trunc(row['species_code'], COL_W['species_code']):<{COL_W['species_code']}}  "
            f"{trunc(row['tissue_cell'], COL_W['tissue_cell']):<{COL_W['tissue_cell']}}  "
            f"{trunc(row['treatment'], COL_W['treatment']):<{COL_W['treatment']}}  "
            f"{grp_disp:<{COL_W['group']}}  "
            f"{trunc(row['file'], COL_W['file'])}"
        )
    print(sep)
    n_sel  = sum(selected)
    n_ctrl = groups.count("C")
    n_test = groups.count("T")
    print(f"  {n_sel} arquivo(s) selecionado(s)  |  Controle: {n_ctrl}  |  Teste: {n_test}")
    print()


# ─────────────────────────────────────────────────────────────────────────────
# Helpers de input
# ─────────────────────────────────────────────────────────────────────────────

def ask(prompt: str, default: str = "") -> str:
    resp = input(prompt).strip()
    return resp if resp else default


def parse_range(s: str, max_idx: int) -> list[int]:
    """
    Interpreta '1,3,5-8,all' → [0,2,4,5,6,7] (0-indexed).
    Retorna lista de índices válidos.
    """
    if s.strip().lower() in ("all", "todos", "*"):
        return list(range(max_idx))
    indices = set()
    for part in re.split(r"[,\s]+", s):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            try:
                indices.update(range(int(a) - 1, int(b)))
            except ValueError:
                pass
        else:
            try:
                indices.add(int(part) - 1)
            except ValueError:
                pass
    return sorted(i for i in indices if 0 <= i < max_idx)


# ─────────────────────────────────────────────────────────────────────────────
# Menu principal interativo
# ─────────────────────────────────────────────────────────────────────────────

MENU_HELP = """
  Comandos disponíveis:
  ─────────────────────────────────────────────────────────────
  s <nums>   Selecionar amostras  (ex: s 1,3,5-8  |  s all)
  d <nums>   Desmarcar amostras   (ex: d 2,4)
  c <nums>   Marcar como CONTROLE (ex: c 1,3)
  t <nums>   Marcar como TESTE    (ex: t 2,4,6-8)
  e <num>    Editar campos de uma amostra manualmente
  r          Resetar seleção e grupos
  p          Visualizar tabela novamente
  ok         Salvar samples.csv e sair
  q          Sair sem salvar
  ─────────────────────────────────────────────────────────────
"""


def edit_row(row: dict) -> dict:
    """Permite ao usuário editar os campos de uma amostra individualmente."""
    print()
    print(f"  Editando: {row['file']}")
    print("  (Enter = manter valor atual)")
    print()

    campos = [
        ("sampleID",    "Sample ID"),
        ("breed",       "Raça (ex: Angus, Nelore, Holstein, Gir, Girolando)"),
        ("species",     "Espécie (ex: Bos taurus, Bos indicus, Híbrido)"),
        ("species_code","Código espécie (taurus/indicus/hybrid)"),
        ("tissue_cell", "Tecido/célula (ex: sangue_total, CD4, CD8, Bcell, fígado)"),
        ("treatment",   "Condição (ex: control, heat_stress, lactating, dry, etc.)"),
    ]

    for key, label in campos:
        val = row.get(key, "?")
        novo = ask(f"  {label} [{val}]: ", default=val)
        row[key] = novo

    return row


def interactive_session(rows: list[dict], out_path: Path) -> None:
    n = len(rows)
    selected = [True] * n   # começa com todos selecionados
    groups   = [""] * n     # "" = sem grupo, "C" = controle, "T" = teste

    # Tenta atribuição automática de grupo a partir do campo treatment
    for i, row in enumerate(rows):
        t = row.get("treatment", "").lower()
        if t in ("control", "controle", "ctrl"):
            groups[i] = "C"
        elif t in ("heat_stress", "heat", "stress", "teste", "case", "treated"):
            groups[i] = "T"

    print()
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║   GenOmics Epigenoma — Gerador de samples.csv               ║")
    print("╚══════════════════════════════════════════════════════════════╝")

    # Detecta parsers usados
    parsers_usados = set(r.get("_parser") for r in rows)
    datasets = set(r.get("dataset") for r in rows)
    print(f"\n  {n} arquivo(s) encontrado(s) | Datasets: {', '.join(datasets)}")
    print(f"  Parsers: {', '.join(parsers_usados)}")

    # Avisa se algum campo ficou como "?"
    n_duvida = sum(
        1 for r in rows
        if "?" in (r.get("breed",""), r.get("treatment",""), r.get("species_code",""))
    )
    if n_duvida:
        print(f"\n  ⚠  {n_duvida} amostra(s) com campos indefinidos (?) — edite com 'e <num>'")

    print_table(rows, selected, groups)
    print(MENU_HELP)

    while True:
        cmd = ask("  > ").strip()
        if not cmd:
            continue

        parts = cmd.split(None, 1)
        op    = parts[0].lower()
        arg   = parts[1] if len(parts) > 1 else ""

        if op == "q":
            print("  Saindo sem salvar.")
            sys.exit(0)

        elif op == "ok":
            break

        elif op == "p":
            print_table(rows, selected, groups)

        elif op == "r":
            selected = [True] * n
            groups   = [""] * n
            print("  Seleção e grupos resetados.")
            print_table(rows, selected, groups)

        elif op == "s":
            idx_list = parse_range(arg, n)
            if not idx_list:
                print("  ⚠ Nenhum índice válido. Ex: s 1,3,5-8  ou  s all")
            else:
                for i in idx_list:
                    selected[i] = True
                print(f"  ✓ {len(idx_list)} amostra(s) selecionada(s).")
                print_table(rows, selected, groups)

        elif op == "d":
            idx_list = parse_range(arg, n)
            for i in idx_list:
                selected[i] = False
                groups[i]   = ""
            print(f"  ✓ {len(idx_list)} amostra(s) desmarcada(s).")
            print_table(rows, selected, groups)

        elif op == "c":
            idx_list = parse_range(arg, n)
            for i in idx_list:
                selected[i] = True
                groups[i]   = "C"
            print(f"  ✓ {len(idx_list)} amostra(s) marcada(s) como CONTROLE.")
            print_table(rows, selected, groups)

        elif op == "t":
            idx_list = parse_range(arg, n)
            for i in idx_list:
                selected[i] = True
                groups[i]   = "T"
            print(f"  ✓ {len(idx_list)} amostra(s) marcada(s) como TESTE.")
            print_table(rows, selected, groups)

        elif op == "e":
            try:
                idx = int(arg.strip()) - 1
                if 0 <= idx < n:
                    rows[idx] = edit_row(rows[idx])
                    print(f"  ✓ Amostra {idx+1} atualizada.")
                    print_table(rows, selected, groups)
                else:
                    print(f"  ✗ Índice fora do intervalo (1–{n}).")
            except ValueError:
                print("  ✗ Use: e <número>  (ex: e 3)")

        elif op in ("h", "help", "?", "ajuda"):
            print(MENU_HELP)

        else:
            print(f"  ✗ Comando desconhecido: '{op}'. Digite 'h' para ajuda.")

    # ── Validação antes de salvar ──────────────────────────────────────────
    sel_rows   = [(rows[i], groups[i]) for i in range(n) if selected[i]]
    n_ctrl     = sum(1 for _, g in sel_rows if g == "C")
    n_test     = sum(1 for _, g in sel_rows if g == "T")
    n_sem_grp  = sum(1 for _, g in sel_rows if g == "")

    print()
    print("─" * 60)
    print(f"  Amostras selecionadas : {len(sel_rows)}")
    print(f"  Controle (C)          : {n_ctrl}")
    print(f"  Teste    (T)          : {n_test}")
    if n_sem_grp:
        print(f"  ⚠  Sem grupo definido  : {n_sem_grp} (serão marcadas como 'undefined')")

    if n_ctrl == 0 or n_test == 0:
        print()
        print("  ⚠  ATENÇÃO: você precisa de pelo menos 1 controle e 1 teste")
        print("     para a análise diferencial methylKit funcionar.")
        confirm = ask("  Salvar assim mesmo? [s/N]: ", "n")
        if confirm.lower() != "s":
            print("  Voltando ao menu... (use 'c' e 't' para definir grupos)")
            interactive_session(rows, out_path)
            return

    save_csv(sel_rows, out_path)


def save_csv(sel_rows: list[tuple[dict, str]], out_path: Path) -> None:
    """Salva o samples.csv no formato esperado pelo run_pipeline.py."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Monta treatment numérico: CONTROLE = 0, TESTE = 1
    final_rows = []
    for row, grp in sel_rows:
        treatment_code = 0 if grp == "C" else (1 if grp == "T" else "undefined")
        final_rows.append({
            "sampleID":    row["sampleID"],
            "breed":       row["breed"],
            "species":     row["species"],
            "species_code": row["species_code"],
            "tissue_cell": row["tissue_cell"],
            "treatment":   treatment_code,
            "condition":   row.get("treatment", "?"),   # rótulo legível
            "month":       row.get("month", "?"),
            "dataset":     row.get("dataset", "?"),
            "file":        row["file"],
        })

    fieldnames = [
        "sampleID", "breed", "species", "species_code",
        "tissue_cell", "treatment", "condition", "month", "dataset", "file"
    ]

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(final_rows)

    print()
    print(f"  ✓ samples.csv salvo em: {out_path}")
    print(f"  {len(final_rows)} amostras exportadas.")
    print()
    print("  Colunas geradas:")
    print("    sampleID     — identificador único da amostra")
    print("    breed        — raça (Angus, Nelore, Holstein, ...)")
    print("    species      — espécie completa")
    print("    species_code — taurus | indicus | hybrid")
    print("    tissue_cell  — tecido ou tipo celular")
    print("    treatment    — 0 = Controle | 1 = Teste (para methylKit)")
    print("    condition    — rótulo legível da condição")
    print("    month        — mês de coleta (se disponível)")
    print("    dataset      — origem dos dados")
    print("    file         — caminho relativo ao --data-dir")
    print()
    print("  Próximo passo:")
    print(f"    python scripts/run_pipeline.py \\")
    print(f"        --metadata {out_path} \\")
    print(f"        --data-dir <pasta_dos_dados> \\")
    print(f"        --out results/run_$(date +%Y%m%d)")
    print()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--data-dir", required=True,
        help="Pasta com os arquivos de metilação (.gz)"
    )
    ap.add_argument(
        "--out", default="samples.csv",
        help="Caminho de saída do samples.csv (padrão: samples.csv)"
    )
    ap.add_argument(
        "--auto", action="store_true",
        help="Modo não-interativo: aceita autoparse sem confirmação"
    )
    ap.add_argument(
        "--dry-run", action="store_true",
        help="Apenas imprime a tabela parseada sem salvar nada"
    )
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    out_path = Path(args.out)

    rows = scan_directory(data_dir)

    if args.dry_run:
        selected = [True] * len(rows)
        groups   = [""] * len(rows)
        for i, row in enumerate(rows):
            t = row.get("treatment", "").lower()
            if t in ("control", "controle", "ctrl"):
                groups[i] = "C"
            elif t not in ("?", ""):
                groups[i] = "T"
        print_table(rows, selected, groups)
        print("  [dry-run] Nenhum arquivo salvo.")
        return

    if args.auto:
        selected = [True] * len(rows)
        groups   = []
        for row in rows:
            t = row.get("treatment", "").lower()
            if t in ("control", "controle", "ctrl"):
                groups.append("C")
            elif t not in ("?", ""):
                groups.append("T")
            else:
                groups.append("")
        sel_rows = [(rows[i], groups[i]) for i in range(len(rows)) if selected[i]]
        save_csv(sel_rows, out_path)
        return

    interactive_session(rows, out_path)


if __name__ == "__main__":
    main()
