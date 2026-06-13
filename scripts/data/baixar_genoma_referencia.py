#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
baixar_genoma_referencia.py — Baixa o genoma de referência ARS-UCD1.3 (Bos taurus)
do Ensembl, descomprime e gera índice samtools faidx.

OBJETIVO:
  Ter o genoma de referência COMPLETO (não só sites variantes do 1000bg) para
  permitir extração de sequências reais nos genes-alvo. Isso vai melhorar a
  qualidade do FASTA do subprojeto de ~5% known para 100% known.

DESTINO:
  E:\\data\\genome_reference\\ARS-UCD1.3\\
    ├── Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.gz   (~800 MB)
    ├── Bos_taurus.ARS-UCD1.3.dna.toplevel.fa      (~3 GB descompactado)
    └── Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.fai  (índice samtools)

USO:
  python baixar_genoma_referencia.py
  python baixar_genoma_referencia.py --skip-decompress  # só baixa
  python baixar_genoma_referencia.py --release 113
"""
from __future__ import annotations
import argparse
import gzip
import os
import shutil
import subprocess
import sys
import time
from datetime import timedelta
from pathlib import Path
import urllib.request

# ============================================================
# CLI
# ============================================================
ap = argparse.ArgumentParser(description="Baixar genoma ARS-UCD1.3 do Ensembl")
ap.add_argument("--dest", default=r"E:\data\genome_reference\ARS-UCD1.3",
                help="Diretório destino (default: E:\\data\\genome_reference\\ARS-UCD1.3)")
ap.add_argument("--release", type=int, default=113,
                help="Release do Ensembl (default: 113)")
ap.add_argument("--skip-decompress", action="store_true",
                help="Não descomprimir após download")
ap.add_argument("--force", action="store_true",
                help="Re-baixar mesmo se já existe")
args = ap.parse_args()

DEST = Path(args.dest)
RELEASE = args.release
URL = f"https://ftp.ensembl.org/pub/release-{RELEASE}/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.gz"
GZ_PATH = DEST / "Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.gz"
FA_PATH = DEST / "Bos_taurus.ARS-UCD1.3.dna.toplevel.fa"
FAI_PATH = DEST / "Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.fai"

# ============================================================
# HELPERS
# ============================================================
def fmt_size(b):
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if abs(b) < 1024:
            return f"{b:.2f} {unit}"
        b /= 1024
    return f"{b:.2f} PB"


def download_with_progress(url, dest_path):
    """Baixa URL com barra de progresso e suporte a resume."""
    dest_path = Path(dest_path)
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    # Resume: se já existe, ler tamanho atual
    local_size = dest_path.stat().st_size if dest_path.exists() else 0

    # HEAD para descobrir tamanho total
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            total_size = int(resp.headers.get("Content-Length", 0))
    except Exception as e:
        print(f"[WARN] HEAD falhou: {e}. Continuando sem total.")
        total_size = 0

    if local_size == total_size and total_size > 0:
        print(f"[OK] Já completo: {dest_path.name} ({fmt_size(total_size)})")
        return True
    if local_size > total_size and total_size > 0:
        print(f"[WARN] Local maior que remoto, deletando: {dest_path}")
        dest_path.unlink()
        local_size = 0

    # GET com Range se resume
    headers = {"User-Agent": "epigenoma-genome-downloader/1.0"}
    if local_size > 0:
        headers["Range"] = f"bytes={local_size}-"
        print(f"[RESUME] Continuando de {fmt_size(local_size)}")

    req = urllib.request.Request(url, headers=headers)
    start = time.time()
    last_print = start
    with urllib.request.urlopen(req, timeout=60) as resp:
        mode = "ab" if local_size > 0 else "wb"
        with open(dest_path, mode) as f:
            done = local_size
            CHUNK = 1024 * 1024  # 1 MB
            while True:
                chunk = resp.read(CHUNK)
                if not chunk:
                    break
                f.write(chunk)
                done += len(chunk)
                now = time.time()
                if now - last_print > 2:
                    elapsed = now - start
                    speed = (done - local_size) / max(elapsed, 0.01)
                    if total_size > 0:
                        pct = 100 * done / total_size
                        eta_s = (total_size - done) / max(speed, 1)
                        print(f"  [{pct:5.1f}%] {fmt_size(done)}/{fmt_size(total_size)} | "
                              f"{fmt_size(speed)}/s | ETA {timedelta(seconds=int(eta_s))}")
                    else:
                        print(f"  {fmt_size(done)} | {fmt_size(speed)}/s")
                    last_print = now

    final_size = dest_path.stat().st_size
    print(f"[OK] Baixado: {dest_path.name} ({fmt_size(final_size)})")
    return total_size == 0 or final_size == total_size


def decompress_gz(gz_path, fa_path):
    """Descomprime .fa.gz mantendo o original."""
    if fa_path.exists():
        print(f"[OK] Já descompactado: {fa_path.name} ({fmt_size(fa_path.stat().st_size)})")
        return True
    print(f"[INFO] Descomprimindo {gz_path.name} -> {fa_path.name}...")
    start = time.time()
    with gzip.open(gz_path, "rb") as src:
        with open(fa_path, "wb") as dst:
            shutil.copyfileobj(src, dst, length=16 * 1024 * 1024)  # 16 MB chunks
    elapsed = time.time() - start
    print(f"[OK] Descomprimido em {timedelta(seconds=int(elapsed))} -> "
          f"{fmt_size(fa_path.stat().st_size)}")
    return True


def index_fasta(fa_path, fai_path):
    """Gera índice samtools faidx (necessário para extração rápida via pysam)."""
    if fai_path.exists() and fai_path.stat().st_mtime >= fa_path.stat().st_mtime:
        print(f"[OK] Índice atualizado: {fai_path.name}")
        return True

    # Tentar samtools primeiro
    if shutil.which("samtools"):
        print(f"[INFO] Indexando via samtools faidx...")
        r = subprocess.run(["samtools", "faidx", str(fa_path)], capture_output=True)
        if r.returncode == 0:
            print(f"[OK] Índice gerado: {fai_path.name}")
            return True

    # Fallback via pysam
    try:
        import pysam
        print(f"[INFO] Indexando via pysam.faidx...")
        pysam.faidx(str(fa_path))
        print(f"[OK] Índice gerado: {fai_path.name}")
        return True
    except ImportError:
        pass
    except Exception as e:
        print(f"[WARN] pysam falhou: {e}")

    # Fallback manual em Python puro (mais lento)
    print(f"[INFO] Indexando manualmente (puro Python, pode demorar)...")
    with open(fa_path, "rb") as f:
        with open(fai_path, "w") as fai:
            offset = 0
            current_name = None
            current_seq_start = None
            current_seq_len = 0
            current_line_blen = None
            current_line_len = None
            for line in f:
                line_len = len(line)
                line_str = line.decode("ascii", errors="replace")
                if line_str.startswith(">"):
                    # Fechar entrada anterior
                    if current_name:
                        fai.write(f"{current_name}\t{current_seq_len}\t{current_seq_start}\t"
                                  f"{current_line_blen}\t{current_line_len}\n")
                    current_name = line_str[1:].strip().split()[0]
                    current_seq_len = 0
                    current_line_blen = None
                    current_line_len = None
                    offset += line_len
                    current_seq_start = offset
                else:
                    # Linha de sequência
                    seq_len = line_len - 1 if line.endswith(b"\n") else line_len
                    if current_line_blen is None:
                        current_line_blen = seq_len
                        current_line_len = line_len
                    current_seq_len += seq_len
                    offset += line_len
            # Fechar última entrada
            if current_name:
                fai.write(f"{current_name}\t{current_seq_len}\t{current_seq_start}\t"
                          f"{current_line_blen}\t{current_line_len}\n")
    print(f"[OK] Índice gerado manualmente: {fai_path.name}")
    return True


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 70)
    print(f"  DOWNLOAD GENOMA REFERÊNCIA — Bos taurus ARS-UCD1.3")
    print(f"  Ensembl release: {RELEASE}")
    print(f"  Destino: {DEST}")
    print("=" * 70)
    print()

    DEST.mkdir(parents=True, exist_ok=True)

    if args.force and GZ_PATH.exists():
        print(f"[FORCE] Removendo {GZ_PATH}")
        GZ_PATH.unlink()

    # ── 1) Download do .fa.gz ──────────────────────────────
    print(f"[1/3] Baixando de {URL}")
    print(f"      -> {GZ_PATH}")
    print(f"      Tamanho esperado: ~800 MB compactado")
    print()
    if not download_with_progress(URL, GZ_PATH):
        print(f"[ERRO] Download falhou")
        sys.exit(1)
    print()

    # ── 2) Descompressão ──────────────────────────────────
    if args.skip_decompress:
        print(f"[2/3] Pulando descompressão (--skip-decompress)")
    else:
        print(f"[2/3] Descomprimindo {GZ_PATH.name}")
        print(f"      -> {FA_PATH}")
        print(f"      Tamanho esperado: ~3 GB descompactado")
        print()
        if not decompress_gz(GZ_PATH, FA_PATH):
            print(f"[ERRO] Descompressão falhou")
            sys.exit(1)
        print()

        # ── 3) Indexação samtools faidx ───────────────────────
        print(f"[3/3] Indexando FASTA")
        print(f"      -> {FAI_PATH}")
        print()
        if not index_fasta(FA_PATH, FAI_PATH):
            print(f"[WARN] Indexação falhou — índice pode ser gerado depois com samtools faidx")
        print()

    # ── Resumo final ──────────────────────────────────────
    print("=" * 70)
    print(f"  CONCLUÍDO")
    print("=" * 70)
    print(f"  GZ:   {GZ_PATH}  ({fmt_size(GZ_PATH.stat().st_size)})")
    if FA_PATH.exists():
        print(f"  FA:   {FA_PATH}  ({fmt_size(FA_PATH.stat().st_size)})")
    if FAI_PATH.exists():
        print(f"  FAI:  {FAI_PATH}  ({fmt_size(FAI_PATH.stat().st_size)})")
    print()
    print(f"  Próximo passo: rodar extract_sequences_from_vcf.py com --reference {FA_PATH}")
    print(f"  (após adaptação do script para usar genoma + integrar variantes)")


if __name__ == "__main__":
    main()
