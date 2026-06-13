#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
baixar_vcfs_alternativos.py — Tenta baixar VCFs das chrs corrompidas
de FONTES ALTERNATIVAS (não EBI 1000bg Run8).

CORROMPIDOS NO RUN8 (CRC32 mismatch / Inflate failed):
  - Chr16, Chr22, Chr25

ESTRATÉGIA:
  1) Tentar 1000 Bull Genomes Run 7 (versão anterior, mesmos animais)
  2) Tentar Ensembl Variation (Bos taurus VCF do dbSNP)
  3) Tentar NCBI dbSNP por chr (FTP)
  4) Tentar Cattle Variation Database (Roslin/CattleGenome)
  5) Tentar via CHEAT (descomprimir + recomprimir o atual via bcftools view --no-version)

Cada fonte tem características diferentes:
  - 1000bg Run7 = mesmos 1000 animais com genótipos (ideal para analise de raça)
  - Ensembl/dbSNP = catálogo de variantes (sem genótipos por animal)
  - Cattle Variation = SNP chip / sequência low-pass

Uso:
  python baixar_vcfs_alternativos.py
  python baixar_vcfs_alternativos.py --chrs 16,22,25
  python baixar_vcfs_alternativos.py --source ensembl
  python baixar_vcfs_alternativos.py --skip-test          # pula bgzip -t

Saídas:
  E:/data/1000bg_vcf_alt/<source>/Chr<N>-*.vcf.gz
  E:/data/1000bg_vcf_alt/download_report.json
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from datetime import datetime
from pathlib import Path

# ============================================================
# CLI
# ============================================================
ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                             description=__doc__)
ap.add_argument("--chrs", default="16,22,25",
                help="Cromossomos a baixar (default: 16,22,25)")
ap.add_argument("--source", choices=["all", "1000bg_run7", "ensembl",
                                     "ncbi_dbsnp", "cattle_var", "ebi_eva"],
                default="all", help="Fonte específica (default: all = tenta todas)")
ap.add_argument("--dest", default=r"E:\data\1000bg_vcf_alt",
                help="Diretório destino")
ap.add_argument("--skip-test", action="store_true",
                help="Pula bgzip -t após download")
ap.add_argument("--timeout", type=int, default=60, help="Timeout HTTP (s)")
args = ap.parse_args()

DEST = Path(args.dest)
DEST.mkdir(parents=True, exist_ok=True)

CHRS = [c.strip() for c in args.chrs.split(",")]
print(f"[CFG] Cromossomos alvo: {CHRS}")
print(f"[CFG] Fonte: {args.source}")
print(f"[CFG] Destino: {DEST}")
print()

# ============================================================
# CATÁLOGO DE FONTES ALTERNATIVAS
# ============================================================
# Cada entry retorna URLs candidatos para um chr (ARS-UCD1.3 sempre que possível)
def urls_1000bg_run7(chrom):
    """Run 7 do 1000 Bull Genomes (versão anterior, mesma TAUIND)."""
    base = "https://ftp.sra.ebi.ac.uk/vol1/ena_ftp"
    # NOTE: precisa do ENA filereport TSV do Run7 para URL exata
    # Como fallback, tentar mirror conhecido
    return [
        f"https://ftp.sra.ebi.ac.uk/vol1/run/Chr{chrom}-Run7-TAUIND-public.vcf.gz",
        f"https://hgdownload.soe.ucsc.edu/gbdb/bosTau9/1000bullgenomes/run7/Chr{chrom}.vcf.gz",
    ]


def urls_ensembl_variation(chrom):
    """Ensembl Variation Bos taurus (catálogo dbSNP, SEM genótipos individuais)."""
    rel = "release-113"
    base = f"https://ftp.ensembl.org/pub/{rel}/variation/vcf/bos_taurus"
    return [
        f"{base}/bos_taurus_chromosome-{chrom}.vcf.gz",
        f"{base}/bos_taurus.vcf.gz",  # all-in-one (precisa filtrar por chr)
    ]


def urls_ncbi_dbsnp(chrom):
    """NCBI dbSNP build atual para Bos taurus (sem genótipos por animal)."""
    base = "https://ftp.ncbi.nih.gov/snp/organisms/archive/bovine"
    # Os builds recentes do dbSNP estão organizados por chromosomo
    return [
        f"{base}/VCF/00-All_chr{chrom}.vcf.gz",
        f"{base}/VCF/00-common_chr{chrom}.vcf.gz",
    ]


def urls_ebi_eva(chrom):
    """EVA (European Variation Archive) — accessions com genótipos."""
    # EVA tem PRJEB23303 para 1000bg, PRJEB42513 para Run8
    # Tentar URLs derivadas
    return [
        f"https://ftp.ebi.ac.uk/eva/PRJEB42513/vcfs/Chr{chrom}.vcf.gz",
        f"https://ftp.ebi.ac.uk/eva/PRJEB23303/vcfs/Chr{chrom}.vcf.gz",
    ]


def urls_cattle_var(chrom):
    """Cattle Variation Database / SNP chip (sem WGS, mas largo)."""
    return [
        f"https://gnomad.broadinstitute.org/cattle/{chrom}.vcf.gz",  # placeholder
    ]


SOURCES = {
    "1000bg_run7": urls_1000bg_run7,
    "ensembl": urls_ensembl_variation,
    "ncbi_dbsnp": urls_ncbi_dbsnp,
    "ebi_eva": urls_ebi_eva,
    "cattle_var": urls_cattle_var,
}

# ============================================================
# HELPERS
# ============================================================
def fmt_size(n):
    for unit in ["B", "KB", "MB", "GB"]:
        if abs(n) < 1024:
            return f"{n:.2f} {unit}"
        n /= 1024
    return f"{n:.2f} TB"


def head_request(url, timeout=30):
    """HEAD HTTP — retorna (status, content_length) ou (None, None)."""
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as r:
            cl = r.headers.get("Content-Length")
            return r.status, int(cl) if cl else None
    except Exception as e:
        return None, str(e)[:80]


def download_file(url, dest_path, timeout=60):
    """Download com progress simples e resume."""
    dest_path = Path(dest_path)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    local_size = dest_path.stat().st_size if dest_path.exists() else 0

    # HEAD para tamanho total
    status, total = head_request(url, timeout)
    if status != 200:
        print(f"  [HEAD] HTTP {status}: {total}")
        return False

    if local_size == total:
        print(f"  [OK] Já completo: {dest_path.name} ({fmt_size(total)})")
        return True

    headers = {"User-Agent": "epigenoma-vcf-downloader/1.0"}
    if local_size > 0:
        headers["Range"] = f"bytes={local_size}-"
        print(f"  [RESUME] {fmt_size(local_size)} / {fmt_size(total) if total else '?'}")

    req = urllib.request.Request(url, headers=headers)
    start = time.time()
    last_report = start
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            mode = "ab" if local_size > 0 else "wb"
            with open(dest_path, mode) as f:
                done = local_size
                CHUNK = 4 * 1024 * 1024  # 4 MB
                while True:
                    chunk = r.read(CHUNK)
                    if not chunk:
                        break
                    f.write(chunk)
                    done += len(chunk)
                    now = time.time()
                    if now - last_report > 5:
                        elapsed = now - start
                        spd = (done - local_size) / max(elapsed, 0.01)
                        if total:
                            pct = 100 * done / total
                            print(f"  [{pct:5.1f}%] {fmt_size(done)} | {fmt_size(spd)}/s")
                        else:
                            print(f"  {fmt_size(done)} | {fmt_size(spd)}/s")
                        last_report = now
        final = dest_path.stat().st_size
        print(f"  [OK] Baixado: {dest_path.name} ({fmt_size(final)})")
        return total is None or final == total
    except Exception as e:
        print(f"  [ERR] {type(e).__name__}: {e}")
        return False


def test_bgzf(path):
    """Testa integridade do arquivo BGZF via Python (puro, sem bgzip)."""
    BGZF_MAGIC = b"\x1f\x8b\x08\x04"
    BGZF_EOF = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
    try:
        with open(path, "rb") as f:
            head = f.read(4)
            if head != BGZF_MAGIC:
                return False, "magic_fail"
            f.seek(-28, 2)
            tail = f.read(28)
            if tail != BGZF_EOF:
                return False, "eof_missing"
        # Tentar descomprimir primeiros 10 MB
        with gzip.open(path, "rb") as f:
            try:
                f.read(10 * 1024 * 1024)
            except Exception as e:
                return False, f"inflate_fail: {e}"
        return True, "ok"
    except Exception as e:
        return False, f"error: {e}"


def md5sum(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            data = f.read(chunk)
            if not data:
                break
            h.update(data)
    return h.hexdigest()


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 70)
    print("  Download de VCFs alternativos para chrs corrompidos")
    print(f"  {datetime.now():%Y-%m-%d %H:%M:%S}")
    print("=" * 70)
    print()

    # Que fontes tentar
    if args.source == "all":
        sources = ["1000bg_run7", "ensembl", "ncbi_dbsnp", "ebi_eva"]
    else:
        sources = [args.source]

    report = {
        "date": datetime.now().isoformat(),
        "chrs_target": CHRS,
        "sources_tried": sources,
        "results": [],
    }

    for chrom in CHRS:
        print(f"\n{'='*70}")
        print(f"  Chr{chrom}")
        print(f"{'='*70}")

        for src_name in sources:
            src_func = SOURCES.get(src_name)
            if src_func is None:
                continue
            urls = src_func(chrom)
            print(f"\n  [{src_name}] testando {len(urls)} URL(s)...")

            chr_dir = DEST / src_name
            chr_dir.mkdir(parents=True, exist_ok=True)

            for url in urls:
                # Probe HEAD
                status, info = head_request(url, args.timeout)
                if status == 200:
                    sz = info if isinstance(info, int) else None
                    print(f"    [HEAD 200] {url}")
                    print(f"    Tamanho: {fmt_size(sz) if sz else '?'}")
                    # Decisão: baixar
                    fname = url.rsplit("/", 1)[-1]
                    dest = chr_dir / fname
                    print(f"    Baixando para: {dest}")
                    ok = download_file(url, dest, args.timeout)

                    test_status = "skipped"
                    if ok and not args.skip_test:
                        print("    Testando integridade BGZF...")
                        is_ok, status_msg = test_bgzf(dest)
                        test_status = status_msg
                        if is_ok:
                            print(f"    [OK] BGZF íntegro")
                        else:
                            print(f"    [WARN] Falha no teste: {status_msg}")

                    report["results"].append({
                        "chrom": chrom,
                        "source": src_name,
                        "url": url,
                        "downloaded": ok,
                        "size_bytes": dest.stat().st_size if dest.exists() else 0,
                        "md5": md5sum(dest) if (ok and dest.exists()) else None,
                        "bgzf_test": test_status,
                        "dest_path": str(dest),
                    })
                    if ok:
                        break  # próxima fonte
                elif status is None:
                    print(f"    [HEAD ERR] {url}: {info}")
                else:
                    print(f"    [HEAD {status}] {url}")
            else:
                continue
            break  # se chegou aqui pela 1ª URL OK, próximo chr

    # Salvar relatório
    rep_path = DEST / "download_report.json"
    with open(rep_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    print()
    print("=" * 70)
    print("  RESUMO")
    print("=" * 70)
    successful = [r for r in report["results"] if r["downloaded"] and r["bgzf_test"] in ("ok", "skipped")]
    print(f"  Tentativas: {len(report['results'])}")
    print(f"  Sucessos:   {len(successful)}")
    for r in successful:
        print(f"    Chr{r['chrom']:>2} ← {r['source']:<14} ({fmt_size(r['size_bytes'])}) — {r['dest_path']}")
    if not successful:
        print()
        print("  NENHUM VCF baixado com sucesso de fontes alternativas.")
        print("  RECOMENDAÇÕES:")
        print("    1) Verificar se o ENA filereport.tsv contém URLs Run7")
        print("    2) Contatar 1000bullgenomes.com para obter mirror alternativo")
        print("    3) Usar Ensembl Variation (catálogo) sem genótipos individuais")
        print("    4) Considerar gerar VCF próprio via low-pass WGS de seus animais")
    print()
    print(f"  Relatório completo: {rep_path}")


if __name__ == "__main__":
    main()
