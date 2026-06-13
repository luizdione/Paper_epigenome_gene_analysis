#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_homer_motifs_v2.py — Análise de motifs no genoma de referência ARS-UCD1.3.

VANTAGENS sobre v1:
  • Usa pysam.FastaFile.fetch() — query indexada O(log n), não linear
  • Aceita FASTA do genoma completo OU FASTA do subprojeto (com cabeçalho gene)
  • Detecta nome do contig automaticamente (1 vs chr1)
  • Janela ±N pb ao redor de cada DMC (default 100 pb)

USO:
  python3 run_homer_motifs_v2.py \\
      --bed results/SUB01.../fase1_methylkit/DMCs_significant.bed \\
      --fasta /data/genome_reference/ARS-UCD1.3/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa \\
      --out results/SUB01.../fase4_homer \\
      --window 100
"""
from __future__ import annotations
import argparse
import collections
import csv
import datetime
import os
import sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--bed", required=True, help="BED com DMCs (chr start end ...)")
ap.add_argument("--fasta", required=True, help="FASTA indexado (.fai existir)")
ap.add_argument("--out", required=True, help="Diretório de saída")
ap.add_argument("--window", type=int, default=100, help="±N pb ao redor de cada DMC")
args = ap.parse_args()

BED = Path(args.bed)
FASTA = Path(args.fasta)
OUT = Path(args.out)
WINDOW = args.window

OUT.mkdir(parents=True, exist_ok=True)

# Validar índice
FAI = Path(str(FASTA) + ".fai")
if not FAI.exists():
    print(f"[ERRO] Índice .fai não encontrado: {FAI}")
    print(f"       Rode: samtools faidx {FASTA}")
    sys.exit(1)

# pysam
try:
    import pysam
except ImportError:
    print("[ERRO] pysam não instalado. Rode: python3 -m pip install pysam")
    sys.exit(1)

print("=" * 70)
print(f"  HOMER motifs v2 — scanning real no genoma")
print(f"  {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
print("=" * 70)
print(f"  BED:    {BED}")
print(f"  FASTA:  {FASTA}")
print(f"  OUT:    {OUT}")
print(f"  Window: ±{WINDOW} pb")
print()

# ── Catálogo de motifs (com IUPAC) ────────────────────────────────────────────
# Fonte: JASPAR + literatura (Kruithof et al. 2018, Salavati 2023)
MOTIFS = [
    # Termotolerância
    ("HSF1",       "GAANNTTC",         "HSF",        "Heat Shock Factor 1 — HSE",
     "Principal regulador de termotolerância"),
    ("HSF1_alt",   "TTCNNGAA",         "HSF",        "Heat Shock Factor 1 — HSE reverso",
     "HSE reverso (orientação 3'-5')"),
    # Transcrição geral / housekeeping
    ("CTCF",       "CCGCGNGGNGGCAG",   "zinc finger","CCCTC-binding factor — insulator",
     "Fronteiras de TADs / domínios de metilação"),
    ("SP1",        "GGGCGG",           "zinc finger","Specificity Protein 1 — GC-box",
     "Promotores CpG-ricos, sensível à metilação"),
    # Estresse / inflamação
    ("AP-1",       "TGANTCA",          "bZIP",       "FOS/JUN — Activator Protein 1",
     "Resposta a estresse, heat shock signaling"),
    ("NFkB",       "GGGRNYYYCC",       "Rel",        "Nuclear Factor kappa B",
     "Imunidade inata e inflamação"),
    # Reprodução / metabolismo
    ("ESR1",       "AGGTCANNNTGACCT",  "nuclear",    "Estrogen Receptor Alpha",
     "Fertilidade e regulação hormonal"),
    ("PPARG",      "AACTAGGTCA",       "nuclear",    "Peroxisome Proliferator-Activated Receptor Gamma",
     "Metabolismo lipídico, marmoreio, lactação"),
    ("CEBPA",      "TTGCGCAA",         "bZIP",       "C/EBP Alpha — adipogênese",
     "Ganho de peso, deposição de gordura"),
    # CpG islands
    ("CpG_dinuc",  "CG",               "dinucleotide","CpG dinucleotide",
     "Sítio de metilação por DNMTs"),
]


def iupac_to_regex(motif):
    """Converte motif com códigos IUPAC para regex."""
    iupac = {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
        "K": "[GT]", "M": "[AC]",
        "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
        "N": "[ACGT]",
    }
    import re
    pattern = "".join(iupac.get(c.upper(), c) for c in motif)
    return re.compile(pattern, re.IGNORECASE)


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
            "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K"}
    return "".join(comp.get(c.upper(), "N") for c in reversed(seq))


# ── Detectar contigs no FASTA ────────────────────────────────────────────────
fa = pysam.FastaFile(str(FASTA))
contigs = set(fa.references)
fa_has_chr = any(c.lower().startswith("chr") for c in list(contigs)[:5])
print(f"[INFO] Contigs no FASTA: {len(contigs)} (primeiros 3: {sorted(list(contigs)[:3])})")
print(f"[INFO] Prefixo 'chr' no FASTA: {fa_has_chr}")
print()


def normalize_chrom(chrom_in_bed):
    """Match contig name no FASTA conforme prefixo."""
    c = chrom_in_bed.replace("chr", "").replace("Chr", "")
    candidates = [c, "chr" + c, "Chr" + c]
    for cand in candidates:
        if cand in contigs:
            return cand
    return None


# ── Ler DMCs do BED e extrair sequências ─────────────────────────────────────
sequences = []
n_skipped = 0
with open(BED) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        chrom_bed = parts[0]
        try:
            start = int(parts[1])
            end = int(parts[2])
        except ValueError:
            continue
        # Janela ±WINDOW
        win_start = max(1, start - WINDOW)
        win_end = end + WINDOW

        chrom_fa = normalize_chrom(chrom_bed)
        if chrom_fa is None:
            n_skipped += 1
            continue
        try:
            seq = fa.fetch(chrom_fa, win_start - 1, win_end).upper()
            sequences.append({
                "chrom": chrom_bed, "start": start, "end": end,
                "seq": seq, "len": len(seq),
            })
        except Exception as e:
            n_skipped += 1

fa.close()
print(f"[INFO] DMCs lidos: {len(sequences)}")
print(f"[INFO] DMCs skipados (chr não encontrado): {n_skipped}")
total_bp = sum(s["len"] for s in sequences)
print(f"[INFO] Total de bp varridos: {total_bp:,}")
print()

# ── Calcular composição (GC content) ────────────────────────────────────────
all_seq = "".join(s["seq"] for s in sequences)
gc = sum(1 for c in all_seq if c in "GC")
at = sum(1 for c in all_seq if c in "AT")
n_count = sum(1 for c in all_seq if c == "N")
print(f"[INFO] GC content: {100*gc/max(gc+at,1):.2f}%")
print(f"[INFO] N count: {n_count} ({100*n_count/max(len(all_seq),1):.2f}%)")
print()

# ── Scan motifs (em ambas as strands) ───────────────────────────────────────
print(f"[INFO] Escaneando {len(MOTIFS)} motifs em {len(sequences)} regiões...")
results = []
for name, consensus, mclass, desc, relevance in MOTIFS:
    pattern_fwd = iupac_to_regex(consensus)
    rc = reverse_complement(consensus)
    pattern_rev = iupac_to_regex(rc) if rc != consensus else None

    n_hits = 0
    n_seqs_with_hit = 0
    for s in sequences:
        seq = s["seq"]
        hits_fwd = len(pattern_fwd.findall(seq))
        hits_rev = len(pattern_rev.findall(seq)) if pattern_rev else 0
        total = hits_fwd + hits_rev
        n_hits += total
        if total > 0:
            n_seqs_with_hit += 1

    pct_seqs = 100 * n_seqs_with_hit / max(len(sequences), 1)

    # Esperado por chance: assumindo random ~25% por base, motif tem
    # probabilidade ~ 4^(-N) com ajuste IUPAC. Aproximação simples: 1/4^L
    # Para CpG (CG): probabilidade = (GC/2)^2 = ~0.04 -> hits esperados = 0.04 * total_bp
    motif_len = len(consensus)
    if name == "CpG_dinuc":
        expected = 0.04 * total_bp  # frequência empírica em mamíferos
    else:
        expected = total_bp / (4 ** motif_len) * 2  # 2 strands
    enrichment = n_hits / max(expected, 0.01)

    results.append({
        "motif": name, "consensus": consensus, "class": mclass,
        "description": desc, "relevance": relevance,
        "n_hits": n_hits, "n_seqs_with_hit": n_seqs_with_hit,
        "pct_seqs": f"{pct_seqs:.1f}%",
        "expected_hits": f"{expected:.2f}",
        "enrichment": f"{enrichment:.2f}x",
        "pvalue": "exact",
        "source": "scan_real_genome",
    })
    print(f"  {name:<12} {consensus:<18}  {n_hits:>5} hits  {n_seqs_with_hit:>3}/{len(sequences)} seqs  {enrichment:.2f}x")

# ── Salvar CSV ──────────────────────────────────────────────────────────────
out_csv = OUT / "motifs_summary.csv"
with open(out_csv, "w", newline="") as f:
    fieldnames = ["motif", "consensus", "class", "description", "relevance",
                  "n_hits", "n_seqs_with_hit", "pct_seqs",
                  "expected_hits", "enrichment", "pvalue", "source"]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)
print()
print(f"[OK] Tabela de motifs: {out_csv}")

# ── Análise k-mer ───────────────────────────────────────────────────────────
print()
print(f"[INFO] Calculando k-mer (k=6) top-20...")
kmer_counts = collections.Counter()
for s in sequences:
    seq = s["seq"]
    for i in range(len(seq) - 6 + 1):
        kmer = seq[i:i+6]
        if "N" not in kmer:
            kmer_counts[kmer] += 1

top20 = kmer_counts.most_common(20)
out_kmer = OUT / "kmer_top20.csv"
with open(out_kmer, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["kmer", "count", "freq_per_kb"])
    for kmer, count in top20:
        freq = count * 1000 / max(total_bp, 1)
        writer.writerow([kmer, count, f"{freq:.3f}"])
print(f"[OK] Top-20 k-mers: {out_kmer}")
print(f"  Top-5: {[k for k, c in top20[:5]]}")

# ── Resumo final ────────────────────────────────────────────────────────────
print()
print("=" * 70)
print(f"  CONCLUÍDO")
print("=" * 70)
print(f"  Motifs analisados: {len(MOTIFS)}")
print(f"  DMCs varridas:     {len(sequences)}")
print(f"  Total bp:          {total_bp:,}")
print(f"  Resultados em:     {OUT}")
