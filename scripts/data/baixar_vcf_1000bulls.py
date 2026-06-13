#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
baixar_vcf_1000bulls.py — Download interativo de VCFs para o projeto Epigenoma Bovino
=======================================================================================
Baixa variantes (SNPs/Indels) de duas fontes:

    MODO COMPLETO  → VCFs por cromossomo do 1000 Bull Genomes Run 9 (EBI/FTP)
    MODO ALVO      → Regiões dos genes-alvo via Ensembl REST, por subprojeto

Ao rodar sem argumentos, o script entra em modo interativo:
    1. Pergunta: completo ou alvo?
    2. Se alvo → lista subprojetos do JSON e pergunta qual(is)
    3. Exibe resumo (genes, regiões, tamanho estimado) e pede confirmação
    4. Só então inicia o download

Flanking padrão: usa `flanking_start` / `flanking_end` já presentes no JSON
(±10 kb em relação ao gene), com opção de expandir via --flanking-extra.

Saída padrão: E:\\data\\epigenomas\\vcfs\\<subprojeto>\\<gene>_variants.vcf
              E:\\data\\epigenomas\\vcfs\\completo\\Run9.all.chr<N>.vcf.gz

Uso direto (sem interação):
    python baixar_vcf_1000bulls.py --mode targeted --subproject thermotolerance
    python baixar_vcf_1000bulls.py --mode targeted --subproject all
    python baixar_vcf_1000bulls.py --mode full --chromosomes 1 2 3
    python baixar_vcf_1000bulls.py --mode check

Requisitos:
    pip install requests tqdm

Autor: Luiz Dione Barbosa De-Melo — GenOmics Epigenoma 2026
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path
from typing import Optional

# ─────────────────────────────────────────────────────────────────────────────
# Configurações — Caminhos padrão no SSD
# ─────────────────────────────────────────────────────────────────────────────

HERE = Path(__file__).resolve().parent

# SSD principal do projeto (Windows)
SSD_BASE          = Path(r"E:\data\epigenomas")
VCF_DIR_TARGETED  = SSD_BASE / "vcfs" / "targeted"
VCF_DIR_FULL      = SSD_BASE / "vcfs" / "completo"

# JSON de genes-alvo (relativo ao script ou absoluto)
DEFAULT_JSON = HERE / "all_target_genes_coordinates.json"

# ─────────────────────────────────────────────────────────────────────────────
# Fontes remotas
# ─────────────────────────────────────────────────────────────────────────────

BULLS_FTP_BASE   = "https://ftp.ebi.ac.uk/pub/databases/1000bullgenomes/Variants/Run9/SNPs_and_Indels"
BULLS_VCF_PAT    = "Run9.all.chr{chrom}.vcf.gz"
BULLS_TBI_PAT    = "Run9.all.chr{chrom}.vcf.gz.tbi"

ENSEMBL_REST     = "https://rest.ensembl.org"
ENSEMBL_SPECIES  = "bos_taurus"

ALL_CHROMOSOMES  = [str(i) for i in range(1, 30)] + ["X"]

# Descrições amigáveis dos módulos
MODULE_LABELS = {
    "milk_production"      : "Produção leiteira",
    "thermo_reproduction"  : "Termotolerância / reprodução",
    "thermotolerance"      : "Termotolerância (core)",
    "thermotolerance_extra": "Termotolerância (extra)",
    "weight_gain"          : "Ganho de peso",
}

# ─────────────────────────────────────────────────────────────────────────────
# Utilitários
# ─────────────────────────────────────────────────────────────────────────────

def try_tqdm():
    try:
        from tqdm import tqdm
        return tqdm
    except ImportError:
        return None


def try_requests():
    try:
        import requests
        return requests
    except ImportError:
        return None


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def human_size(n: int) -> str:
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f} {u}"
        n /= 1024
    return f"{n:.1f} PB"


def sep(char="─", n=60):
    print(char * n)


def ask(prompt: str, valid: list[str] | None = None, default: str = "") -> str:
    """Lê input do usuário com validação opcional."""
    while True:
        resp = input(prompt).strip()
        if not resp and default:
            return default
        if valid is None or resp.lower() in [v.lower() for v in valid]:
            return resp
        print(f"  Opção inválida. Escolha entre: {', '.join(valid)}")


# ─────────────────────────────────────────────────────────────────────────────
# Carrega JSON e extrai subprojetos
# ─────────────────────────────────────────────────────────────────────────────

def load_genes(json_path: Path) -> dict:
    if not json_path.exists():
        print(f"\n✗ JSON não encontrado: {json_path}")
        print("  Coloque all_target_genes_coordinates.json na pasta do script.")
        sys.exit(1)
    with open(json_path, encoding="utf-8") as f:
        return json.load(f)


def get_subprojects(genes: dict) -> dict[str, list[str]]:
    """Retorna {modulo: [lista de genes]} ordenado."""
    modules: dict[str, list[str]] = {}
    for gene, info in genes.items():
        mod = info.get("module", "sem_modulo")
        modules.setdefault(mod, []).append(gene)
    return dict(sorted(modules.items()))


def genes_for_subproject(all_genes: dict, subproject: str) -> dict:
    """Filtra genes de um subprojeto específico (ou todos se 'all')."""
    if subproject == "all":
        return all_genes
    return {g: i for g, i in all_genes.items() if i.get("module") == subproject}


# ─────────────────────────────────────────────────────────────────────────────
# Download com retry + resume
# ─────────────────────────────────────────────────────────────────────────────

def download_file(
    url: str,
    dest: Path,
    resume: bool = True,
    max_retries: int = 5,
    chunk_size: int = 1024 * 1024,
) -> bool:
    tqdm = try_tqdm()
    part = dest.with_suffix(dest.suffix + ".part")
    existing = part.stat().st_size if (resume and part.exists()) else 0

    headers: dict[str, str] = {}
    if existing > 0:
        headers["Range"] = f"bytes={existing}-"
        print(f"    [resume] continuando de {human_size(existing)}")

    for attempt in range(1, max_retries + 1):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=60) as resp:
                total = int(resp.headers.get("Content-Length", 0)) + existing
                mode  = "ab" if existing > 0 else "wb"

                with open(part, mode) as fh:
                    pbar = (tqdm(total=total, initial=existing, unit="B",
                                 unit_scale=True, unit_divisor=1024,
                                 desc=dest.name[:40], leave=False)
                            if tqdm else None)
                    downloaded = existing
                    while True:
                        chunk = resp.read(chunk_size)
                        if not chunk:
                            break
                        fh.write(chunk)
                        downloaded += len(chunk)
                        if pbar:
                            pbar.update(len(chunk))
                    if pbar:
                        pbar.close()

            part.rename(dest)
            print(f"    ✓ {dest.name} ({human_size(dest.stat().st_size)})")
            return True

        except (urllib.error.URLError, OSError) as e:
            wait = 2 ** attempt
            print(f"    [tentativa {attempt}/{max_retries}] {e} — aguardando {wait}s...")
            if attempt < max_retries:
                time.sleep(wait)
                existing = part.stat().st_size if part.exists() else 0
                if existing > 0:
                    headers["Range"] = f"bytes={existing}-"
            else:
                print(f"    ✗ Falha definitiva: {dest.name}")
                return False
    return False


# ─────────────────────────────────────────────────────────────────────────────
# Modo CHECK
# ─────────────────────────────────────────────────────────────────────────────

def mode_check() -> None:
    print("\n[CHECK] Testando acesso ao servidor 1000 Bull Genomes (Run 9)...")
    print(f"  URL base: {BULLS_FTP_BASE}\n")
    ok_list = []
    for chrom in ["1", "2", "X"]:
        url = f"{BULLS_FTP_BASE}/{BULLS_VCF_PAT.format(chrom=chrom)}"
        try:
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=15) as r:
                size = int(r.headers.get("Content-Length", 0))
                print(f"  ✓ chr{chrom:>2} disponível — {human_size(size)}")
                ok_list.append(chrom)
        except urllib.error.HTTPError as e:
            print(f"  ✗ chr{chrom:>2} — HTTP {e.code}")
        except urllib.error.URLError as e:
            print(f"  ✗ chr{chrom:>2} — {e.reason}")

    print()
    if ok_list:
        print("✓ Servidor acessível. Modo --mode full disponível.")
        print("  Estimativa: ~5–50 GB por cromossomo (Run 9 completo).")
    else:
        print("✗ Servidor inacessível agora. Use --mode targeted (Ensembl REST).")


# ─────────────────────────────────────────────────────────────────────────────
# Modo FULL — cromossomos inteiros do 1000 Bull Genomes
# ─────────────────────────────────────────────────────────────────────────────

def mode_full(chromosomes: list[str], out_dir: Path, resume: bool) -> None:
    ensure_dir(out_dir)
    print(f"\n[FULL] {len(chromosomes)} cromossomo(s) → {out_dir}")
    sep()
    ok = fail = 0
    for chrom in chromosomes:
        vcf_name = BULLS_VCF_PAT.format(chrom=chrom)
        tbi_name = BULLS_TBI_PAT.format(chrom=chrom)
        vcf_dest = out_dir / vcf_name
        tbi_dest = out_dir / tbi_name
        print(f"\n  Cromossomo {chrom}")
        if vcf_dest.exists():
            print(f"    ○ Já existe — pulando (delete o arquivo para forçar re-download)")
            ok += 1
        else:
            if download_file(f"{BULLS_FTP_BASE}/{vcf_name}", vcf_dest, resume=resume):
                ok += 1
            else:
                fail += 1
                continue
        if not tbi_dest.exists():
            download_file(f"{BULLS_FTP_BASE}/{tbi_name}", tbi_dest,
                          resume=False)
    sep()
    print(f"[FULL] Concluído: {ok} OK | {fail} falha(s)")
    if fail:
        print("  → Re-execute para retomar downloads interrompidos.")


# ─────────────────────────────────────────────────────────────────────────────
# Modo TARGETED — regiões dos genes-alvo via Ensembl REST
# ─────────────────────────────────────────────────────────────────────────────

def mode_targeted(
    genes_subset: dict,
    subproject_name: str,
    out_dir: Path,
    flanking_extra: int = 0,
) -> None:
    """
    Para cada gene em genes_subset, baixa variantes da região via Ensembl REST.

    Flanking aplicado:
        start = flanking_start (já inclui ±10 kb do JSON) - flanking_extra
        end   = flanking_end   (já inclui ±10 kb do JSON) + flanking_extra

    Gera:
        <out_dir>/<subprojeto>/<gene>_variants.vcf   — por gene
        <out_dir>/<subprojeto>/all_variants.vcf       — consolidado
    """
    requests = try_requests()
    if requests is None:
        print("\nERRO: módulo 'requests' não instalado.")
        print("  Execute: pip install requests")
        sys.exit(1)

    sub_dir = ensure_dir(out_dir / subproject_name)
    total   = len(genes_subset)

    print(f"\n[TARGETED] Subprojeto: {subproject_name}  |  {total} genes")
    print(f"  Flanking JSON (±10 kb) + extra {flanking_extra:,} pb")
    print(f"  Fonte: Ensembl REST ({ENSEMBL_SPECIES})")
    print(f"  Saída: {sub_dir}")
    sep()

    session = requests.Session()
    session.headers.update({"Content-Type": "text/x-vcf"})

    ok = fail = skip = 0
    all_lines: list[str] = []
    header_captured = False

    for idx, (gene, info) in enumerate(genes_subset.items(), 1):
        chrom        = info.get("chromosome")
        f_start      = info.get("flanking_start") or info.get("start")
        f_end        = info.get("flanking_end")   or info.get("end")
        raw_start    = info.get("start", f_start)
        raw_end      = info.get("end",   f_end)

        if not all([chrom, f_start, f_end]):
            print(f"  [{idx:02d}/{total}] {gene:<18} ○ coordenadas incompletas — pulando")
            skip += 1
            continue

        # Janela final: flanking do JSON (±10 kb) + extra solicitado pelo usuário
        q_start = max(1, int(f_start) - flanking_extra)
        q_end   = int(f_end) + flanking_extra
        region_size = q_end - q_start

        out_file = sub_dir / f"{gene}_variants.vcf"
        if out_file.exists() and out_file.stat().st_size > 0:
            n_vars = sum(1 for l in out_file.read_text().splitlines()
                         if l and not l.startswith("#"))
            print(f"  [{idx:02d}/{total}] {gene:<18} ○ já existe ({n_vars} variantes)")
            ok += 1
            for l in out_file.read_text().splitlines():
                if not l.startswith("#"):
                    all_lines.append(l)
            continue

        url = (f"{ENSEMBL_REST}/overlap/region/{ENSEMBL_SPECIES}"
               f"/{chrom}:{q_start}-{q_end}?feature=variation")

        for attempt in range(1, 4):
            try:
                resp = session.get(url, timeout=40)

                if resp.status_code == 200:
                    content = resp.text
                    lines   = content.splitlines()
                    n_vars  = sum(1 for l in lines if l and not l.startswith("#"))

                    if not header_captured:
                        for l in lines:
                            if l.startswith("#"):
                                all_lines.append(l)
                        header_captured = True

                    out_file.write_text(content, encoding="utf-8")
                    all_lines.extend(l for l in lines if l and not l.startswith("#"))

                    gene_region = (f"chr{chrom}:{int(raw_start):,}-{int(raw_end):,}"
                                   f"  [±10kb+{flanking_extra}pb → {region_size/1000:.0f}kb total]")
                    print(f"  [{idx:02d}/{total}] {gene:<18} ✓  {n_vars:>6} variantes  {gene_region}")
                    ok += 1
                    break

                elif resp.status_code == 429:
                    wait = int(resp.headers.get("Retry-After", 60))
                    print(f"  [{idx:02d}/{total}] {gene:<18} ⏸  rate limit — aguardando {wait}s...")
                    time.sleep(wait)

                elif resp.status_code == 400:
                    print(f"  [{idx:02d}/{total}] {gene:<18} ✗  região inválida para Ensembl (HTTP 400)")
                    fail += 1
                    break

                else:
                    print(f"  [{idx:02d}/{total}] {gene:<18} ✗  HTTP {resp.status_code} (tentativa {attempt}/3)")
                    time.sleep(2 ** attempt)

            except Exception as e:
                print(f"  [{idx:02d}/{total}] {gene:<18} ✗  {e} (tentativa {attempt}/3)")
                time.sleep(2 ** attempt)
        else:
            fail += 1

        time.sleep(0.4)  # respeita rate limit da Ensembl (~15 req/s)

    # Arquivo consolidado do subprojeto
    consolidated = sub_dir / "all_variants.vcf"
    consolidated.write_text("\n".join(all_lines), encoding="utf-8")
    n_total = sum(1 for l in all_lines if l and not l.startswith("#"))

    sep()
    print(f"[TARGETED] Concluído: {ok} OK | {fail} falha(s) | {skip} pulado(s)")
    print(f"  VCF consolidado: {consolidated}")
    print(f"  Total de variantes: {n_total:,}")
    print()
    print("  Próximo passo — ativar Fase 6 (MOFA2):")
    print(f"  python scripts/run_pipeline.py --metadata samples.csv \\")
    print(f"      --data-dir /data/NCBI_GSE161113 --out results/run_mofa2 --skip 4")


# ─────────────────────────────────────────────────────────────────────────────
# Modo INTERATIVO — pergunta ao usuário antes de baixar
# ─────────────────────────────────────────────────────────────────────────────

def interactive_mode(json_path: Path, out_dir_targeted: Path, out_dir_full: Path) -> None:
    data      = load_genes(json_path)
    all_genes = data.get("genes", {})
    subprojects = get_subprojects(all_genes)

    print()
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║   GenOmics Epigenoma — Download de VCFs (modo interativo)   ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    print()
    print("  Escolha o tipo de download:")
    print()
    print("  [1] ALVO   — variantes nas regiões dos genes-alvo (Ensembl REST)")
    print("               Leve (~10–500 MB total). Recomendado para MOFA2.")
    print()
    print("  [2] COMPLETO — VCFs cromossômicos do 1000 Bull Genomes Run 9")
    print("               Pesado (~5–50 GB por cromossomo).")
    print()
    choice = ask("  Sua escolha [1/2]: ", valid=["1", "2"])

    # ── MODO COMPLETO ─────────────────────────────────────────────────────────
    if choice == "2":
        print()
        print("  Cromossomos disponíveis: 1–29, X")
        print("  Digite os cromossomos separados por espaço,")
        print("  ou 'all' para todos (atenção: centenas de GB!).")
        print()
        raw = input("  Cromossomos [ex: 1 2 3 X] ou 'all': ").strip()
        if raw.lower() == "all":
            chroms = ALL_CHROMOSOMES
        else:
            chroms = [c.strip() for c in raw.split() if c.strip()]
            invalid = [c for c in chroms if c not in ALL_CHROMOSOMES]
            if invalid:
                print(f"\n  ✗ Cromossomos inválidos: {invalid}")
                print(f"  Válidos: {', '.join(ALL_CHROMOSOMES)}")
                sys.exit(1)

        print()
        sep()
        print(f"  RESUMO — Download COMPLETO")
        sep()
        print(f"  Cromossomos : {', '.join(chroms)}")
        print(f"  Fonte       : {BULLS_FTP_BASE}")
        print(f"  Destino     : {out_dir_full}")
        print(f"  Estimativa  : ~{len(chroms) * 15} GB (varia por cromossomo)")
        sep()
        confirm = ask("\n  Confirmar download? [s/N]: ", valid=["s", "n", "S", "N"], default="n")
        if confirm.lower() != "s":
            print("\n  Download cancelado.")
            sys.exit(0)
        mode_full(chroms, out_dir_full, resume=True)

    # ── MODO ALVO ─────────────────────────────────────────────────────────────
    else:
        print()
        print("  Subprojetos disponíveis no JSON:")
        print()
        sp_keys = list(subprojects.keys())
        for i, mod in enumerate(sp_keys, 1):
            label  = MODULE_LABELS.get(mod, mod)
            n_genes = len(subprojects[mod])
            genes_preview = ", ".join(subprojects[mod][:4])
            if n_genes > 4:
                genes_preview += f" ... (+{n_genes - 4})"
            print(f"  [{i}] {label:<35} ({n_genes} genes)")
            print(f"       {genes_preview}")
            print()
        print(f"  [0] TODOS os subprojetos ({len(all_genes)} genes no total)")
        print()

        valid_choices = ["0"] + [str(i) for i in range(1, len(sp_keys) + 1)]
        raw = ask(f"  Escolha o subprojeto [{'/'.join(valid_choices)}]: ",
                  valid=valid_choices)

        if raw == "0":
            selected_mod   = "all"
            selected_genes = all_genes
        else:
            selected_mod   = sp_keys[int(raw) - 1]
            selected_genes = genes_for_subproject(all_genes, selected_mod)

        # Flanking extra
        print()
        print("  Flanking já embutido no JSON: ±10.000 pb por gene.")
        print("  Deseja adicionar flanking extra além dos 10 kb? (0 = não)")
        raw_flank = input("  Flanking extra em pb [padrão 0]: ").strip()
        flanking_extra = int(raw_flank) if raw_flank.isdigit() else 0

        # Resumo antes de confirmar
        print()
        sep()
        print(f"  RESUMO — Download ALVO")
        sep()
        if selected_mod == "all":
            print(f"  Subprojeto  : TODOS")
        else:
            print(f"  Subprojeto  : {MODULE_LABELS.get(selected_mod, selected_mod)}")
        print(f"  Genes       : {len(selected_genes)}")
        print()
        for gene, info in selected_genes.items():
            chrom   = info.get("chromosome", "?")
            f_start = info.get("flanking_start") or info.get("start", 0)
            f_end   = info.get("flanking_end")   or info.get("end", 0)
            q_start = max(1, int(f_start) - flanking_extra)
            q_end   = int(f_end) + flanking_extra
            kb      = (q_end - q_start) / 1000
            print(f"    {gene:<18} chr{chrom}  {int(q_start):>12,} – {int(q_end):>12,}  ({kb:.0f} kb)")
        print()
        print(f"  Flanking    : ±10.000 pb (JSON) + {flanking_extra:,} pb extra")
        print(f"  Fonte       : Ensembl REST API ({ENSEMBL_SPECIES})")
        print(f"  Destino     : {out_dir_targeted / selected_mod}")
        sep()
        confirm = ask("\n  Confirmar download? [s/N]: ", valid=["s", "n", "S", "N"], default="n")
        if confirm.lower() != "s":
            print("\n  Download cancelado.")
            sys.exit(0)

        mode_targeted(selected_genes, selected_mod, out_dir_targeted,
                      flanking_extra=flanking_extra)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--mode", choices=["interactive", "targeted", "full", "check"],
        default="interactive",
        help="Modo de operação (padrão: interactive)"
    )
    ap.add_argument(
        "--json", default=str(DEFAULT_JSON),
        help="Caminho para all_target_genes_coordinates.json"
    )
    ap.add_argument(
        "--subproject", default=None,
        help="Subprojeto para modo targeted (nome do módulo ou 'all')"
    )
    ap.add_argument(
        "--chromosomes", nargs="+", default=None,
        help="Cromossomos para modo full (ex: 1 2 3 X)"
    )
    ap.add_argument(
        "--flanking-extra", type=int, default=0,
        help="pb extras além dos ±10 kb já embutidos no JSON (modo targeted)"
    )
    ap.add_argument(
        "--out-targeted", default=str(VCF_DIR_TARGETED),
        help=f"Diretório de saída para modo targeted (padrão: {VCF_DIR_TARGETED})"
    )
    ap.add_argument(
        "--out-full", default=str(VCF_DIR_FULL),
        help=f"Diretório de saída para modo full (padrão: {VCF_DIR_FULL})"
    )
    ap.add_argument(
        "--no-resume", action="store_true",
        help="Forçar re-download mesmo se arquivo já existe"
    )
    args = ap.parse_args()

    json_path      = Path(args.json)
    out_targeted   = Path(args.out_targeted)
    out_full       = Path(args.out_full)
    resume         = not args.no_resume

    if args.mode == "check":
        mode_check()
        return

    if args.mode == "interactive":
        interactive_mode(json_path, out_targeted, out_full)
        return

    if args.mode == "full":
        chroms = args.chromosomes or ALL_CHROMOSOMES
        invalid = [c for c in chroms if c not in ALL_CHROMOSOMES]
        if invalid:
            print(f"✗ Cromossomos inválidos: {invalid}")
            sys.exit(1)
        print(f"\n⚠  Modo full: {len(chroms)} cromossomo(s) — pode ser centenas de GB!")
        confirm = ask("  Confirmar? [s/N]: ", valid=["s", "n", "S", "N"], default="n")
        if confirm.lower() != "s":
            sys.exit(0)
        mode_full(chroms, out_full, resume=resume)
        return

    if args.mode == "targeted":
        data      = load_genes(json_path)
        all_genes = data.get("genes", {})
        sub       = args.subproject or "all"
        subset    = genes_for_subproject(all_genes, sub)
        if not subset:
            print(f"✗ Subprojeto '{sub}' não encontrado. Disponíveis:")
            for k in get_subprojects(all_genes):
                print(f"  {k}")
            sys.exit(1)
        mode_targeted(subset, sub, out_targeted,
                      flanking_extra=args.flanking_extra)


if __name__ == "__main__":
    main()
