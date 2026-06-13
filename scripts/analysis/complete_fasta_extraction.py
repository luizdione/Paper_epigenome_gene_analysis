#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Completar extração de sequências FASTA — apenas genes faltantes.
Appenda ao arquivo existente sem re-processar o que já foi escrito.

Uso:
  # Modo automático (detecta subprojeto ativo do JSON):
  python3 scripts/complete_fasta_extraction.py

  # Modo explícito:
  python3 scripts/complete_fasta_extraction.py --sub SUB01 \
      --vcf-dir /data/1000bg_vcf --min-density 0.0

Notas:
  - Paths parametrizados por --sub (default: lê subprojeto/subprojeto_config.json)
  - MIN_DENSITY default = 0.0 (aceita tudo, mesmo sequências all-N) para evitar
    falso skip de genes em chr X / chrs com baixa cobertura no 1000 Bull Genomes
"""
from __future__ import annotations
import argparse
import datetime
import gzip
import json
import re
import subprocess
import sys
import zlib
from pathlib import Path

# pysam path
for _pip_path in [
    '/usr/local/lib/python3.10/dist-packages',
    '/usr/lib/python3/dist-packages',
    '/home/rstudio/.local/lib/python3.10/site-packages',
]:
    if _pip_path not in sys.path:
        sys.path.insert(0, _pip_path)

# ─── CLI ──────────────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser(description="Append sequências FASTA faltantes (genes do subprojeto)")
ap.add_argument("--sub", default=None, help="ID do subprojeto (ex: SUB01). Default: lê config.")
ap.add_argument("--base", default="/home/rstudio/projeto", help="Diretório base do projeto")
ap.add_argument("--vcf-dir", default="/data/1000bg_vcf", help="Diretório com VCFs do 1000 Bull Genomes")
ap.add_argument("--window", type=int, default=500, help="Flanking ±N pb ao redor do gene")
ap.add_argument("--min-density", type=float, default=0.0,
                help="Densidade mínima de bases ACGT no FASTA (0.0 = aceitar tudo)")
args = ap.parse_args()

BASE = Path(args.base)

# Detectar subprojeto ativo do config se não passado
if args.sub is None:
    cfg_path = BASE / "subprojeto" / "subprojeto_config.json"
    try:
        cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
        args.sub = cfg.get("_meta", {}).get("id") or "SUB01"
    except Exception:
        args.sub = "SUB01"

SUB = args.sub
FA_PATH   = BASE / "subprojeto" / SUB / f"{SUB}_reference_sequences.fa"
JSON_PATH = BASE / "subprojeto" / SUB / f"{SUB}_target_genes_coordinates.json"
VCF_DIR   = Path(args.vcf_dir)
WINDOW    = args.window
MIN_DENSITY = args.min_density

print(f"[CFG] Subprojeto: {SUB}")
print(f"[CFG] FASTA: {FA_PATH}")
print(f"[CFG] JSON:  {JSON_PATH}")
print(f"[CFG] VCF dir: {VCF_DIR}")
print(f"[CFG] Window: ±{WINDOW} pb | min_density={MIN_DENSITY}")
print()

# ─── Carregar genes já escritos ───────────────────────────────────────────────
written_genes = set()
if FA_PATH.exists():
    with open(FA_PATH) as f:
        for line in f:
            if line.startswith(">"):
                gene = line.strip()[1:].split("_chr")[0]
                written_genes.add(gene)
print(f"Genes já escritos: {len(written_genes)} → {sorted(written_genes)[:5]}...")

# ─── Carregar JSON ────────────────────────────────────────────────────────────
with open(JSON_PATH) as f:
    data = json.load(f)
genes = data.get("genes", {})
if isinstance(genes, list):
    genes = {g.get("gene_name", f"gene_{i}"): g for i, g in enumerate(genes)}

missing = {k: v for k, v in genes.items() if k not in written_genes}
print(f"Genes faltantes: {len(missing)}")
for g, info in sorted(missing.items(), key=lambda x: str(x[1].get("chromosome", ""))):
    print(f"  {g} chr{info.get('chromosome')}")

# ─── Helpers ──────────────────────────────────────────────────────────────────
def norm_chrom(c):
    c = str(c)
    for pfx in ["Chr", "chr", "CHR"]:
        if c.startswith(pfx):
            c = c[len(pfx):]
    return c

def find_vcf(chrom_norm):
    patterns = [
        re.compile(rf'[Cc]hr{re.escape(chrom_norm)}[-_\.]', re.IGNORECASE),
        re.compile(rf'[Cc]hr0*{re.escape(chrom_norm)}[-_\.]', re.IGNORECASE),
    ]
    for vcf in sorted(VCF_DIR.glob("*.vcf.gz")):
        for pat in patterns:
            if pat.search(vcf.name):
                return vcf
    return None

def ensure_index(vcf_path):
    tbi = Path(str(vcf_path) + ".tbi")
    if tbi.exists():
        return True
    try:
        r = subprocess.run(
            ["tabix", "-p", "vcf", str(vcf_path)],
            capture_output=True, timeout=120
        )
        return r.returncode == 0
    except Exception:
        return False

def detect_prefix(vcf_path):
    try:
        import pysam
        tbx = pysam.TabixFile(str(vcf_path))
        contigs = list(tbx.contigs)
        tbx.close()
        if contigs:
            c = contigs[0]
            if c.startswith("Chr"): return "Chr"
            if c.startswith("chr"): return "chr"
        return ""
    except Exception:
        return ""

def extract_tabix(vcf_path, regions, prefix):
    import pysam
    positions = {}
    try:
        tbx = pysam.TabixFile(str(vcf_path))
        vcf_chroms = set(tbx.contigs)
        for reg in regions:
            chrom_norm = norm_chrom(reg["chrom"])
            candidates = [prefix + chrom_norm, "Chr" + chrom_norm, chrom_norm]
            chrom_vcf = next((c for c in candidates if c in vcf_chroms), None)
            if not chrom_vcf:
                continue
            start = reg["start"] - 1  # tabix 0-based
            end   = reg["end"]
            try:
                for rec in tbx.fetch(chrom_vcf, max(0, start), end):
                    parts = rec.split("\t")
                    if len(parts) < 4:
                        continue
                    try:
                        pos = int(parts[1])
                    except ValueError:
                        continue
                    ref = parts[3].strip()
                    if len(ref) == 1 and ref.upper() in "ACGT":
                        positions[pos] = ref.upper()
            except Exception:
                pass
        tbx.close()
    except Exception as e:
        print(f"  [WARN] tabix error: {e}", file=sys.stderr)
    return positions

def extract_linear(vcf_path, regions):
    positions = {}
    intervals = [(r["start"], r["end"]) for r in regions]
    min_pos = min(s for s, e in intervals)
    max_pos = max(e for s, e in intervals)
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    print(f"  → [linear] {vcf_path.name} ({min_pos:,}–{max_pos:,})", flush=True)
    try:
        with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t", 5)
                if len(parts) < 4:
                    continue
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue
                if pos < min_pos or pos > max_pos:
                    continue
                ref = parts[3].strip()
                if len(ref) != 1 or ref.upper() not in "ACGT":
                    continue
                for s, e in intervals:
                    if s <= pos <= e:
                        positions[pos] = ref.upper()
                        break
    except (OSError, EOFError, zlib.error, Exception) as ex:
        print(f"  [WARN] Erro linear: {ex}", file=sys.stderr)
    print(f"  → {len(positions):,} posições coletadas.", flush=True)
    return positions

def reconstruct(positions, start, end):
    seq = []
    for pos in range(start, end + 1):
        seq.append(positions.get(pos, "N"))
    return "".join(seq)

# ─── Verificar pysam ──────────────────────────────────────────────────────────
has_pysam = False
try:
    import pysam
    has_pysam = True
    print("[INFO] pysam disponível → tabix ativado")
except ImportError:
    print("[INFO] pysam não disponível → fallback linear")

# ─── Agrupar por cromossomo ───────────────────────────────────────────────────
by_chrom = {}
for gene_name, info in missing.items():
    chrom = norm_chrom(info.get("chromosome", ""))
    start = info.get("start", 0) - WINDOW
    end   = info.get("end", 0) + WINDOW
    label = f"{gene_name}_chr{chrom}:{start}-{end}"
    entry = {"gene": gene_name, "chrom": chrom, "start": start, "end": end, "label": label}
    by_chrom.setdefault(chrom, []).append(entry)

print(f"\nProcessando {len(by_chrom)} cromossomos em append mode...\n")
written = 0
skipped = 0

with open(FA_PATH, "a", encoding="ascii") as fa:
    for chrom_norm, regions in sorted(by_chrom.items(), key=lambda x: x[0].zfill(3)):
        vcf = find_vcf(chrom_norm)
        if vcf is None:
            print(f"[WARN] Sem VCF para chr{chrom_norm} → skip {len(regions)} genes", file=sys.stderr)
            skipped += len(regions)
            continue

        t0 = datetime.datetime.now()
        print(f"[{t0:%H:%M:%S}] Chr{chrom_norm} ({len(regions)} genes) | {vcf.name}", flush=True)

        positions = {}
        if has_pysam:
            indexed = ensure_index(vcf)
            if indexed:
                prefix = detect_prefix(vcf)
                positions = extract_tabix(vcf, regions, prefix)
                print(f"  → [tabix] {len(positions):,} posições REF", flush=True)
            else:
                positions = extract_linear(vcf, regions)
        else:
            positions = extract_linear(vcf, regions)

        for reg in regions:
            seq = reconstruct(positions, reg["start"], reg["end"])
            known = sum(seq.count(b) for b in "ACGT")
            density = known / max(len(seq), 1)
            if density < MIN_DENSITY:
                print(f"  [skip] {reg['gene']}: densidade {density:.2%} < {MIN_DENSITY:.2%}")
                skipped += 1
                continue
            fa.write(f">{reg['label']}\n")
            for i in range(0, len(seq), 60):
                fa.write(seq[i:i+60] + "\n")
            written += 1
            elapsed = (datetime.datetime.now() - t0).total_seconds()
            print(f"  ✓ {reg['gene']} ({len(seq):,} bp, {density:.2%} known) [{elapsed:.1f}s]", flush=True)

print(f"\n=== CONCLUÍDO ===")
print(f"  Escritas: {written}")
print(f"  Ignoradas: {skipped}")
print(f"  FASTA: {FA_PATH}")

# Contagem final
with open(FA_PATH) as f:
    total = sum(1 for l in f if l.startswith(">"))
print(f"  Total no arquivo: {total} sequências")
