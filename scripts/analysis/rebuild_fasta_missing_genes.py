#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reconstrução cirúrgica do FASTA do subprojeto:
1) Remove entradas com 100% de N (extrações falhas anteriores)
2) Re-extrai os genes faltantes usando pysam.TabixFile (query indexada)
3) Reindexa via samtools faidx (ou cria .fai manualmente)

Uso:
  python3 scripts/rebuild_fasta_missing_genes.py --sub SUB01 \
      --vcf-dir /data/1000bg_vcf --window 500
"""
from __future__ import annotations
import argparse
import datetime
import json
import re
import subprocess
import sys
from pathlib import Path

import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--sub", default=None)
ap.add_argument("--base", default=".")
ap.add_argument("--vcf-dir", default="/data/1000bg_vcf")
ap.add_argument("--window", type=int, default=500)
ap.add_argument("--min-density", type=float, default=0.0)
args = ap.parse_args()

BASE = Path(args.base)
if args.sub is None:
    cfg = json.loads((BASE / "subprojeto" / "subprojeto_config.json").read_text(encoding="utf-8"))
    args.sub = cfg.get("_meta", {}).get("id", "SUB01")
SUB = args.sub
FA   = BASE / "subprojeto" / SUB / f"{SUB}_reference_sequences.fa"
JSON = BASE / "subprojeto" / SUB / f"{SUB}_target_genes_coordinates.json"
VCF_DIR = Path(args.vcf_dir)

print(f"[CFG] SUB={SUB}  FA={FA}  VCF={VCF_DIR}")

# ── 1) Carregar FASTA atual e identificar entradas all-N ──────────────────────
print("\n[1/3] Lendo FASTA atual e detectando entradas all-N...")
keep = {}      # gene_name → (header, seq) para genes com pelo menos 1 base ACGT
all_n = []     # genes com 100% N (para re-extrair)
current_header = None
current_seq = []
def flush():
    if not current_header:
        return
    name = current_header[1:].split("_chr")[0]
    seq = "".join(current_seq)
    n_known = sum(1 for c in seq.upper() if c in "ACGT")
    if n_known == 0:
        all_n.append(name)
    else:
        keep[name] = (current_header, seq)

if FA.exists():
    with open(FA) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                flush()
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        flush()

print(f"  Genes com sequência válida: {len(keep)}")
print(f"  Genes all-N (a re-extrair): {len(all_n)} → {sorted(all_n)}")

# ── 2) Identificar genes faltantes (não estão nem em keep nem em all_n) ──────
data = json.loads(JSON.read_text(encoding="utf-8"))
target_genes = data.get("genes", {})
present = set(keep.keys()) | set(all_n)
missing = [g for g in target_genes.keys() if g not in present]
print(f"  Genes ausentes do FASTA: {len(missing)} → {missing}")

# Lista final de genes a (re)extrair
to_extract = sorted(set(all_n) | set(missing))
print(f"\n[2/3] Genes a (re)extrair via tabix: {len(to_extract)}")

# ── 3) Extrair via pysam tabix ────────────────────────────────────────────────
def norm_chrom(c):
    return re.sub(r"^[Cc][Hh][Rr]", "", str(c))

def find_vcf(chrom_norm):
    for vcf in sorted(VCF_DIR.glob("*.vcf.gz")):
        m = re.search(r"[Cc]hr0*(\w+?)[-_\.]", vcf.name)
        if m and norm_chrom(m.group(1)).upper() == chrom_norm.upper():
            return vcf
    return None

def detect_prefix(vcf_path):
    try:
        tbx = pysam.TabixFile(str(vcf_path))
        contigs = list(tbx.contigs)
        tbx.close()
        if contigs and contigs[0].lower().startswith("chr"):
            return contigs[0][:3]  # "Chr" ou "chr"
        return ""
    except Exception:
        return ""

def extract_via_tabix(vcf_path, chrom_norm, start, end, prefix):
    """Retorna dict {pos: REF} usando query tabix."""
    positions = {}
    try:
        tbx = pysam.TabixFile(str(vcf_path))
        chrom_vcf_candidates = [
            prefix + chrom_norm,
            "Chr" + chrom_norm,
            "chr" + chrom_norm,
            chrom_norm,
        ]
        chrom_vcf = next((c for c in chrom_vcf_candidates if c in tbx.contigs), None)
        if chrom_vcf is None:
            print(f"  [WARN] chr{chrom_norm} não encontrado em {vcf_path.name}; contigs disponíveis: {list(tbx.contigs)[:5]}")
            tbx.close()
            return positions
        for rec in tbx.fetch(chrom_vcf, max(0, start - 1), end):
            parts = rec.split("\t", 5)
            if len(parts) < 4:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            ref = parts[3].strip()
            if len(ref) == 1 and ref.upper() in "ACGT":
                positions[pos] = ref.upper()
        tbx.close()
    except Exception as e:
        print(f"  [ERROR] tabix em {vcf_path.name}: {e}", file=sys.stderr)
    return positions

extracted = {}  # gene → (header, seq)
for gene in to_extract:
    info = target_genes.get(gene, {})
    chrom = norm_chrom(info.get("chromosome", ""))
    s = int(info.get("start", 0)) - args.window
    e = int(info.get("end", 0)) + args.window
    label = f"{gene}_chr{chrom}:{s}-{e}"

    vcf = find_vcf(chrom)
    if vcf is None:
        print(f"  [SKIP] {gene}: VCF para chr{chrom} não encontrado")
        continue
    prefix = detect_prefix(vcf)
    t0 = datetime.datetime.now()
    print(f"  [{t0:%H:%M:%S}] {gene} chr{chrom}:{s:,}-{e:,} ({vcf.name})")
    positions = extract_via_tabix(vcf, chrom, s, e, prefix)
    seq = "".join(positions.get(p, "N") for p in range(s, e + 1))
    n_known = sum(1 for c in seq if c in "ACGT")
    density = n_known / max(len(seq), 1)
    elapsed = (datetime.datetime.now() - t0).total_seconds()
    print(f"    → {len(seq):,} bp, {n_known:,} known ({density:.2%}) [{elapsed:.1f}s]")
    if density < args.min_density:
        print(f"    [SKIP] densidade < {args.min_density:.2%}")
        continue
    extracted[gene] = (f">{label}", seq)

# ── 3) Reescrever FASTA ──────────────────────────────────────────────────────
print(f"\n[3/3] Reescrevendo FASTA: keep={len(keep)} + new={len(extracted)}")
all_entries = {**keep, **extracted}
# Ordenar por chromosomo+start conforme aparece no JSON
order = list(target_genes.keys())
ordered = [g for g in order if g in all_entries]
extras = [g for g in all_entries if g not in order]
final_keys = ordered + extras

with open(FA, "w", encoding="ascii") as f:
    for gene in final_keys:
        header, seq = all_entries[gene]
        f.write(header + "\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

# ── Reindexar com samtools faidx ─────────────────────────────────────────────
print("\n[reindex] Tentando samtools faidx...")
try:
    r = subprocess.run(["samtools", "faidx", str(FA)], capture_output=True, timeout=120)
    if r.returncode == 0:
        print(f"  ✓ {FA}.fai gerado via samtools")
    else:
        print(f"  [WARN] samtools faidx retornou {r.returncode}: {r.stderr.decode()[:200]}")
        # Fallback: gerar .fai manualmente via pysam
        pysam.faidx(str(FA))
        print(f"  ✓ {FA}.fai gerado via pysam.faidx")
except FileNotFoundError:
    pysam.faidx(str(FA))
    print(f"  ✓ {FA}.fai gerado via pysam.faidx (samtools não disponível)")
except Exception as e:
    print(f"  [WARN] reindex falhou: {e}")

# ── Resumo final ─────────────────────────────────────────────────────────────
total = sum(1 for line in open(FA) if line.startswith(">"))
print(f"\n=== CONCLUÍDO ===")
print(f"  Total no FASTA: {total} sequências (esperado: {len(target_genes)})")
print(f"  FASTA: {FA}")
print(f"  FAI:   {FA}.fai")
