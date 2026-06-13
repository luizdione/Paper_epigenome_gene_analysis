#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
baixar_ensembl_variation.py — Baixa Ensembl Variation Bos taurus (1.03 GB)
e extrai por cromossomo INTEIRO (catálogo reutilizável por TODOS os subprojetos).

ESTRATÉGIA EM 2 NÍVEIS (one-time download, multi-subproject use):

  NÍVEL 1 (este script): Baixa o catálogo completo Ensembl uma única vez,
                         splita por cromossomo (todos os 30 chrs + X + MT).
                         → Salva em E:\\data\\ensembl_variation\\bos_taurus_v113\\
                         → ~1.5 GB total (todos os 31 chrs splitados + .tbi)
                         → Reutilizado por TODOS os subprojetos futuros (SUB01, SUB02, ...)

  NÍVEL 2 (snps_in_cpg.py): Cada subprojeto consome só os chrs dos seus genes-alvo
                            via tabix query (rápido, sem re-baixar).
                            → tabix Chr16-Ensembl.vcf.gz chr16:25032085-25040832

CONTEXTO:
  Os VCFs Chr16/22/25 do 1000 Bull Genomes Run8 (TAUIND) estão corrompidos
  com erros BGZF (CRC32 mismatch / Inflate failed) APÓS 2 re-downloads.
  Ensembl Variation = fallback robusto, catálogo de variantes (baseado em dbSNP).

LIMITAÇÕES:
  • Ensembl Variation NÃO contém GENÓTIPOS por animal (só REF/ALT, AF global)
  • Para análise de divergência Angus vs Nelore precisaríamos dos genótipos
  • Mas serve para: identificar SNPs comuns nos genes-alvo, mapear localização,
    cruzar com DMCs e motifs, frequência alélica global de população cattle

URL:
  https://ftp.ensembl.org/pub/release-113/variation/vcf/bos_taurus/bos_taurus.vcf.gz
  Tamanho: 1.034.206.231 bytes (~1 GB compactado, ~5 GB descompactado)

USO:
  python baixar_ensembl_variation.py                       # baixa + splita TODOS os chrs (default)
  python baixar_ensembl_variation.py --chrs 16,22,25       # só chrs específicos (legado)
  python baixar_ensembl_variation.py --release 110         # versão alternativa
  python baixar_ensembl_variation.py --skip-download       # só processa
"""
from __future__ import annotations
import argparse
import gzip
import os
import sys
import time
import urllib.request
from datetime import datetime, timedelta
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--release", type=int, default=113, help="Release Ensembl (default: 113)")
ap.add_argument("--dest", default=r"E:\data\ensembl_variation",
                help="Diretório destino reutilizável (default: E:\\data\\ensembl_variation)")
ap.add_argument("--chrs", default="all",
                help="Chrs a extrair: 'all' (todos 1-29, X, MT) OU lista '16,22,25' (default: all)")
ap.add_argument("--skip-download", action="store_true",
                help="Pula download, só extrai do arquivo já presente")
args = ap.parse_args()

# Subdiretório por release para versionamento
DEST = Path(args.dest) / f"bos_taurus_v{args.release}"
DEST.mkdir(parents=True, exist_ok=True)

URL = f"https://ftp.ensembl.org/pub/release-{args.release}/variation/vcf/bos_taurus/bos_taurus.vcf.gz"
GZ_PATH = DEST / f"bos_taurus_release-{args.release}.vcf.gz"

# Resolver chrs: "all" expande para 1-29, X, MT
if args.chrs.lower() == "all":
    CHRS_SET = set([str(i) for i in range(1, 30)] + ["X", "Y", "MT"])
else:
    CHRS_SET = set(c.strip() for c in args.chrs.split(","))


def fmt(n):
    for u in ["B", "KB", "MB", "GB"]:
        if abs(n) < 1024:
            return f"{n:.2f} {u}"
        n /= 1024
    return f"{n:.2f} TB"


def download_with_resume(url, dest):
    """Download com resume e progress."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    local_size = dest.stat().st_size if dest.exists() else 0

    # HEAD para tamanho total
    try:
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=30) as r:
            total = int(r.headers.get("Content-Length", 0))
    except Exception as e:
        print(f"[ERR] HEAD falhou: {e}")
        return False

    print(f"[INFO] Tamanho total: {fmt(total)}")
    if local_size == total and total > 0:
        print(f"[OK] Já completo: {dest.name}")
        return True
    if local_size > total:
        print(f"[WARN] Local maior que remoto, deletando")
        dest.unlink()
        local_size = 0

    headers = {"User-Agent": "epigenoma-ensembl/1.0"}
    if local_size > 0:
        headers["Range"] = f"bytes={local_size}-"
        print(f"[RESUME] {fmt(local_size)} / {fmt(total)} ({100*local_size/total:.1f}%)")

    req = urllib.request.Request(url, headers=headers)
    start = time.time()
    last_p = start
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            mode = "ab" if local_size > 0 else "wb"
            with open(dest, mode) as f:
                done = local_size
                CHUNK = 4 * 1024 * 1024
                while True:
                    data = r.read(CHUNK)
                    if not data:
                        break
                    f.write(data)
                    done += len(data)
                    now = time.time()
                    if now - last_p > 5:
                        elapsed = now - start
                        spd = (done - local_size) / max(elapsed, 0.01)
                        pct = 100 * done / total if total else 0
                        eta_s = (total - done) / max(spd, 1)
                        print(f"  [{pct:5.1f}%] {fmt(done)}/{fmt(total)} | "
                              f"{fmt(spd)}/s | ETA {timedelta(seconds=int(eta_s))}")
                        last_p = now
        print(f"[OK] Baixado: {fmt(dest.stat().st_size)}")
        return True
    except Exception as e:
        print(f"[ERR] Download: {e}")
        return False


def extract_chr_variants(gz_path, chrs_set, out_dir):
    """Lê o VCF combinado e splita em arquivos por cromossomo (lazy: só os pedidos)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Detectar quais já estão prontos (skip se .gz + .tbi)
    todo = []
    for chrom in chrs_set:
        gz = out_dir / f"Chr{chrom}-Ensembl-bos_taurus.vcf.gz"
        tbi = out_dir / f"Chr{chrom}-Ensembl-bos_taurus.vcf.gz.tbi"
        if gz.exists() and tbi.exists() and gz.stat().st_size > 1000:
            print(f"  [SKIP] Chr{chrom}: já existe ({fmt(gz.stat().st_size)} + .tbi)")
        else:
            todo.append(chrom)

    if not todo:
        print(f"[INFO] Todos os chrs solicitados já estão extraídos.")
        return {c: 0 for c in chrs_set}

    print(f"[INFO] Chrs a extrair: {sorted(todo)}")

    # Abrir arquivos de saída por chr
    out_files = {}
    for chrom in todo:
        out_path = out_dir / f"Chr{chrom}-Ensembl-bos_taurus.vcf"
        out_files[chrom] = open(out_path, "w", encoding="utf-8")

    # Header escrito apenas uma vez por arquivo
    print(f"[INFO] Lendo {gz_path} (~1 GB compactado, ~5-10 min de processamento)")
    n_total = 0
    counters = {c: 0 for c in todo}
    header_lines = []
    start_time = time.time()
    last_report = start_time

    with gzip.open(gz_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line)
                continue
            n_total += 1
            chrom = line.split("\t", 1)[0]
            chrom_norm = chrom.replace("chr", "").replace("Chr", "")
            if chrom_norm in counters:
                f_out = out_files[chrom_norm]
                if counters[chrom_norm] == 0:
                    f_out.writelines(header_lines)
                f_out.write(line)
                counters[chrom_norm] += 1

            now = time.time()
            if now - last_report > 30:
                elapsed = now - start_time
                rate = n_total / max(elapsed, 0.01)
                # Estatística de progresso por chr (só dos top 5 mais lidos)
                top = sorted(counters.items(), key=lambda x: -x[1])[:5]
                top_str = ", ".join(f"Chr{c}={n:,}" for c, n in top)
                print(f"  Lidos {n_total:,} | {timedelta(seconds=int(elapsed))} "
                      f"({rate:,.0f}/s) | top: {top_str}")
                last_report = now

    for f_out in out_files.values():
        f_out.close()

    print()
    print(f"[INFO] Total de variantes lidas: {n_total:,}")
    for chrom in sorted(counters.keys()):
        n = counters[chrom]
        if n > 0:
            print(f"  Chr{chrom}: {n:,} variantes extraídas")
        else:
            print(f"  Chr{chrom}: 0 variantes (chr não presente no VCF Ensembl)")
            # Remover arquivo vazio
            try:
                (out_dir / f"Chr{chrom}-Ensembl-bos_taurus.vcf").unlink(missing_ok=True)
            except Exception:
                pass

    return counters


def compress_and_index(vcf_path):
    """Comprime via Python (gzip) e tenta gerar .tbi via tabix se disponível."""
    import shutil
    gz_path = vcf_path + ".gz"

    print(f"  Comprimindo {Path(vcf_path).name} (gzip puro Python)...")
    with open(vcf_path, "rb") as f_in:
        with gzip.open(gz_path, "wb", compresslevel=6) as f_out:
            shutil.copyfileobj(f_in, f_out, length=4 * 1024 * 1024)

    os.remove(vcf_path)
    print(f"  [OK] Comprimido: {fmt(Path(gz_path).stat().st_size)}")

    # Tentar tabix se disponível
    if shutil.which("tabix"):
        try:
            import subprocess
            r = subprocess.run(["tabix", "-p", "vcf", gz_path],
                               capture_output=True, timeout=300)
            if r.returncode == 0:
                print(f"  [OK] Índice .tbi gerado")
            else:
                print(f"  [WARN] tabix falhou: {r.stderr.decode()[:200]}")
        except Exception as e:
            print(f"  [WARN] tabix erro: {e}")


def main():
    print("=" * 70)
    print(f"  Ensembl Variation Bos taurus (release {args.release})")
    print(f"  {datetime.now():%Y-%m-%d %H:%M:%S}")
    print("=" * 70)
    print(f"  URL:    {URL}")
    print(f"  Local:  {GZ_PATH}")
    print(f"  Chrs:   {sorted(CHRS_SET)}")
    print()

    # 1) Download
    if not args.skip_download:
        print("[1/3] Baixando bos_taurus.vcf.gz (~1 GB)...")
        if not download_with_resume(URL, GZ_PATH):
            sys.exit(1)
        print()
    else:
        print("[1/3] Pulando download (--skip-download)")
        if not GZ_PATH.exists():
            print(f"[ERR] Arquivo não existe: {GZ_PATH}")
            sys.exit(1)

    # 2) Extrair chrs
    print(f"[2/3] Extraindo Chr{','.join(sorted(CHRS_SET))}...")
    counters = extract_chr_variants(GZ_PATH, CHRS_SET, DEST)
    print()

    # 3) Comprimir + indexar
    print("[3/3] Comprimindo + indexando arquivos extraídos...")
    for chrom in sorted(CHRS_SET):
        vcf_file = DEST / f"Chr{chrom}-Ensembl-bos_taurus.vcf"
        if vcf_file.exists() and counters[chrom] > 0:
            compress_and_index(str(vcf_file))
    print()

    print("=" * 70)
    print("  CONCLUÍDO")
    print("=" * 70)
    for chrom in sorted(CHRS_SET):
        gz = DEST / f"Chr{chrom}-Ensembl-bos_taurus.vcf.gz"
        tbi = DEST / f"Chr{chrom}-Ensembl-bos_taurus.vcf.gz.tbi"
        if gz.exists():
            print(f"  Chr{chrom}: {fmt(gz.stat().st_size)} "
                  f"({'.tbi OK' if tbi.exists() else 'sem .tbi'})")
    print()
    print("  ⚠️  IMPORTANTE: Ensembl Variation NÃO contém genótipos por animal.")
    print("     Use estes VCFs para: identificar SNPs nos genes-alvo, AF global,")
    print("     cruzar com DMCs/motifs. NÃO serve para análise de divergência")
    print("     Angus vs Nelore (precisa dos VCFs do 1000bg ou similar).")


if __name__ == "__main__":
    main()
