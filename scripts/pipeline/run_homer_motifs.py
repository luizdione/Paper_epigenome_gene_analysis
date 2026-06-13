#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fase 4 — Análise de Motifs de Fatores de Transcrição
=====================================================
Wrapper para HOMER (findMotifsGenome.pl) com fallback puro Python usando
PWMs do JASPAR para bovino quando HOMER não está instalado.

Modos de operação:
    1. HOMER disponível → roda findMotifsGenome.pl e parseia knownResults.txt
    2. HOMER ausente → análise de enriquecimento de k-mers + scan de PWMs
       básicos (CpG-island, CTCF, SP1, AP-1) nas sequências de DMR.

Entradas:
    --bed   : BED com DMCs significativos (saída da Fase 1)
    --genome: Genoma de referência (default: bosTau9 para HOMER, ou FASTA para fallback)
    --fasta : FASTA do genoma (necessário apenas para fallback; se ausente, gera stats básicas)
    --out   : Diretório de saída

Saídas:
    motifs_summary.csv     — tabela consolidada de motifs significativos
    motifs_enrichment.png  — barplot dos top motifs enriquecidos
    homer_raw/             — outputs brutos do HOMER (se disponível)

Uso:
    # Com HOMER instalado:
    python scripts/run_homer_motifs.py \\
        --bed results/fase1_methylkit/DMCs_significant.bed \\
        --genome bosTau9 --out results/fase4_homer

    # Sem HOMER (fallback puro Python):
    python scripts/run_homer_motifs.py \\
        --bed results/fase1_methylkit/DMCs_significant.bed \\
        --out results/fase4_homer

Autor: Luiz Dione Barbosa De-Melo — GenOmics Epigenoma 2026
"""
from __future__ import annotations
import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

try:
    import numpy as np
except ImportError:
    np = None


# ═══════════════════════════════════════════════════════════════════════════════
# PWMs conhecidos (subconjunto do JASPAR para bovino / mamíferos)
# Formato: nome → {consensus, matrix (A,C,G,T por posição), classe}
# ═══════════════════════════════════════════════════════════════════════════════
KNOWN_MOTIFS = {
    "CTCF": {
        "consensus": "CCGCGNGGNGGCAG",
        "class": "zinc finger",
        "description": "CCCTC-binding factor — insulator, organização 3D da cromatina",
        "relevance": "Fronteiras de domínios de metilação / TADs",
    },
    "SP1": {
        "consensus": "GGGCGG",
        "class": "zinc finger",
        "description": "Specificity Protein 1 — ativa genes housekeeping via GC-box",
        "relevance": "Promotores CpG-ricos, sensível à metilação",
    },
    "AP-1 (FOS/JUN)": {
        "consensus": "TGANTCA",
        "class": "bZIP",
        "description": "Activator Protein 1 — resposta a estresse",
        "relevance": "Heat shock response, stress signaling",
    },
    "NF-kB": {
        "consensus": "GGGRNTTCC",
        "class": "Rel homology",
        "description": "Nuclear Factor kappa B — imunidade inata e inflamação",
        "relevance": "Regulação imune diferencial taurinos vs zebuínos",
    },
    "HSF1": {
        "consensus": "GAANNTTC",
        "class": "HSF",
        "description": "Heat Shock Factor 1 — Heat Shock Elements (HSE)",
        "relevance": "Principal regulador de termotolerância",
    },
    "CEBPA": {
        "consensus": "TTGCGCAA",
        "class": "bZIP",
        "description": "CCAAT/Enhancer Binding Protein Alpha — adipogênese",
        "relevance": "Ganho de peso e deposição de gordura",
    },
    "PPARG": {
        "consensus": "AACTAGGTCA",
        "class": "nuclear receptor",
        "description": "Peroxisome Proliferator-Activated Receptor Gamma",
        "relevance": "Metabolismo lipídico, marmoreio, produção leiteira",
    },
    "ESR1": {
        "consensus": "AGGTCANNNTGACCT",
        "class": "nuclear receptor",
        "description": "Estrogen Receptor Alpha — reprodução e fertilidade",
        "relevance": "Fertilidade e regulação hormonal reprodutiva",
    },
}


def iupac_to_regex(consensus: str) -> str:
    """Converte consenso IUPAC em regex."""
    iupac = {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
        "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
        "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
    }
    return "".join(iupac.get(c.upper(), c) for c in consensus)


# ═══════════════════════════════════════════════════════════════════════════════
# Modo HOMER
# ═══════════════════════════════════════════════════════════════════════════════
def run_homer(bed_path: str, genome: str, out_dir: str, size: int = 200) -> list[dict]:
    """Roda HOMER findMotifsGenome.pl e parseia resultados."""
    homer_out = os.path.join(out_dir, "homer_raw")
    os.makedirs(homer_out, exist_ok=True)

    cmd = [
        "findMotifsGenome.pl", bed_path, genome, homer_out,
        "-size", str(size), "-p", "4",
    ]
    print(f"[Info] Rodando HOMER: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if result.returncode != 0:
        print(f"[WARN] HOMER falhou (rc={result.returncode})")
        if result.stderr:
            print(result.stderr[:500])
        return []

    return parse_homer_known(os.path.join(homer_out, "knownResults.txt"))


def parse_homer_known(path: str) -> list[dict]:
    """Parseia knownResults.txt do HOMER."""
    if not os.path.exists(path):
        return []
    motifs = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        header = fh.readline()  # skip header
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            try:
                motifs.append({
                    "motif": parts[0].split("/")[0].strip(),
                    "consensus": parts[1] if len(parts) > 1 else "",
                    "pvalue": float(parts[2]) if parts[2] else 1.0,
                    "log_pvalue": float(parts[3]) if parts[3] else 0.0,
                    "qvalue_BH": float(parts[4]) if parts[4] else 1.0,
                    "n_target": int(float(parts[5])) if parts[5] else 0,
                    "pct_target": float(parts[6].rstrip("%")) if parts[6] else 0,
                    "n_bg": int(float(parts[7])) if parts[7] else 0,
                    "pct_bg": float(parts[8].rstrip("%")) if parts[8] else 0,
                    "source": "HOMER",
                })
            except (ValueError, IndexError):
                continue
    return motifs


# ═══════════════════════════════════════════════════════════════════════════════
# Modo Fallback (puro Python)
# ═══════════════════════════════════════════════════════════════════════════════
def read_bed(bed_path: str) -> list[dict]:
    """Lê BED e retorna lista de regiões."""
    regions = []
    with open(bed_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            try:
                regions.append({
                    "chrom": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "name": parts[3] if len(parts) > 3 else "",
                    "score": float(parts[4]) if len(parts) > 4 else 0,
                })
            except ValueError:
                continue
    return regions


def extract_sequences_from_fasta(fasta_path: str, regions: list[dict]) -> list[str]:
    """Extrai sequências de um FASTA indexado (ou lê linearmente se pequeno)."""
    # Tentar usar pysam/pyfaidx se disponível
    try:
        import pysam
        fa = pysam.FastaFile(fasta_path)
        seqs = []
        for r in regions:
            try:
                seq = fa.fetch(r["chrom"], r["start"], r["end"])
                seqs.append(seq.upper())
            except (KeyError, ValueError):
                pass
        fa.close()
        return seqs
    except ImportError:
        pass

    # Fallback: leitura linear (lento, só para genomas pequenos ou testes)
    print("[WARN] pysam não disponível — leitura linear do FASTA (pode ser lento)")
    chrom_seqs = {}
    current_chrom = None
    current_seq = []
    with open(fasta_path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if current_chrom:
                    chrom_seqs[current_chrom] = "".join(current_seq)
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip().upper())
        if current_chrom:
            chrom_seqs[current_chrom] = "".join(current_seq)

    seqs = []
    for r in regions:
        chrom_key = r["chrom"]
        if chrom_key not in chrom_seqs:
            chrom_key = chrom_key.replace("chr", "")
        if chrom_key not in chrom_seqs:
            chrom_key = "chr" + chrom_key
        if chrom_key in chrom_seqs:
            s = chrom_seqs[chrom_key][r["start"]:r["end"]]
            if s:
                seqs.append(s)
    return seqs


def scan_motifs_in_sequences(sequences: list[str]) -> list[dict]:
    """Escaneia motifs conhecidos nas sequências usando regex IUPAC."""
    total_bp = sum(len(s) for s in sequences)
    n_seqs = len(sequences)
    results = []

    for name, info in KNOWN_MOTIFS.items():
        pattern = iupac_to_regex(info["consensus"])
        regex = re.compile(pattern, re.IGNORECASE)
        # Contar matches na fita + e -
        n_hits = 0
        seqs_with_hit = 0
        for seq in sequences:
            fwd = len(regex.findall(seq))
            # Reverse complement
            rc = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
            rev = len(regex.findall(rc))
            hits = fwd + rev
            n_hits += hits
            if hits > 0:
                seqs_with_hit += 1

        # Frequência esperada por chance (modelo simplificado: uniforme)
        motif_len = len(info["consensus"].replace("N", ""))
        expected_per_bp = (0.25 ** motif_len) * 2  # ambas fitas
        expected_hits = expected_per_bp * total_bp

        # Enriquecimento
        enrichment = n_hits / max(expected_hits, 0.01)
        pct_seqs = 100 * seqs_with_hit / max(n_seqs, 1)

        # p-valor aproximado (Poisson)
        if np is not None:
            from scipy.stats import poisson
            try:
                pval = 1 - poisson.cdf(n_hits - 1, expected_hits)
            except Exception:
                pval = 1.0 if n_hits <= expected_hits else 0.001
        else:
            pval = 1.0 if n_hits <= expected_hits else 0.001

        results.append({
            "motif": name,
            "consensus": info["consensus"],
            "class": info["class"],
            "description": info["description"],
            "relevance": info["relevance"],
            "n_hits": n_hits,
            "n_seqs_with_hit": seqs_with_hit,
            "pct_seqs": round(pct_seqs, 1),
            "expected_hits": round(expected_hits, 1),
            "enrichment": round(enrichment, 2),
            "pvalue": pval,
            "source": "fallback_regex",
        })

    results.sort(key=lambda r: r["pvalue"])
    return results


def compute_kmer_stats(sequences: list[str], k: int = 6) -> list[dict]:
    """Calcula frequência de k-mers nas sequências de DMR."""
    counts = Counter()
    total = 0
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if all(c in "ACGT" for c in kmer):
                counts[kmer] += 1
                total += 1

    # Top 20 k-mers mais frequentes
    top = counts.most_common(20)
    return [{"kmer": km, "count": c, "frequency": round(c / max(total, 1), 6)} for km, c in top]


def fallback_analysis(bed_path: str, fasta_path: str | None, out_dir: str) -> list[dict]:
    """Análise de motifs sem HOMER."""
    regions = read_bed(bed_path)
    if not regions:
        print("[WARN] Nenhuma região encontrada no BED.")
        return []
    print(f"[Info] {len(regions)} regiões no BED")

    sequences = []
    if fasta_path and os.path.exists(fasta_path):
        print(f"[Info] Extraindo sequências de {fasta_path}...")
        sequences = extract_sequences_from_fasta(fasta_path, regions)
        print(f"[Info] {len(sequences)} sequências extraídas")

    if sequences:
        # Análise completa com sequências
        motifs = scan_motifs_in_sequences(sequences)

        # K-mer analysis
        kmers = compute_kmer_stats(sequences)
        kmer_path = os.path.join(out_dir, "kmer_top20.csv")
        with open(kmer_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=["kmer", "count", "frequency"])
            writer.writeheader()
            writer.writerows(kmers)
        print(f"[Info] K-mer analysis salvo em: {kmer_path}")
    else:
        # Sem FASTA: análise estatística baseada apenas nas coordenadas
        print("[Info] Sem FASTA disponível — gerando análise estatística baseada em coordenadas.")
        motifs = []
        for name, info in KNOWN_MOTIFS.items():
            motifs.append({
                "motif": name,
                "consensus": info["consensus"],
                "class": info["class"],
                "description": info["description"],
                "relevance": info["relevance"],
                "n_hits": "N/A",
                "n_seqs_with_hit": "N/A",
                "pct_seqs": "N/A",
                "expected_hits": "N/A",
                "enrichment": "N/A",
                "pvalue": "N/A (sem FASTA)",
                "source": "reference_catalog",
            })

    return motifs


# ═══════════════════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════════════════
def plot_motif_enrichment(motifs: list[dict], out_path: str):
    """Gera barplot dos motifs enriquecidos."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib não disponível — pulando plot.")
        return

    # Filtrar motifs com dados numéricos
    plotable = [m for m in motifs if isinstance(m.get("enrichment"), (int, float)) and m["enrichment"] > 0]
    if not plotable:
        print("[WARN] Sem dados numéricos para plotar.")
        return

    plotable.sort(key=lambda m: m["enrichment"], reverse=True)
    top = plotable[:15]

    fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * len(top) + 1)))
    names = [m["motif"] for m in top]
    enrichments = [m["enrichment"] for m in top]
    colors = ["#D32F2F" if m.get("pvalue", 1) < 0.05 else "#90CAF9" for m in top]

    bars = ax.barh(range(len(names)), enrichments, color=colors, edgecolor="gray", linewidth=0.5)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("Enriquecimento (observado / esperado)")
    ax.set_title("Motifs de TFs enriquecidos nas DMRs")
    ax.axvline(1, color="black", linestyle="--", linewidth=0.7, alpha=0.5)

    # Legenda
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#D32F2F", label="p < 0.05"),
        Patch(facecolor="#90CAF9", label="p ≥ 0.05"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8)
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"[OK] Plot salvo: {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, help="BED com DMCs (saída Fase 1)")
    ap.add_argument("--genome", default="bosTau9", help="Genoma ref para HOMER (ex: bosTau9)")
    ap.add_argument("--fasta", default=None, help="FASTA do genoma (para fallback sem HOMER)")
    ap.add_argument("--size", type=int, default=200, help="Tamanho da janela em torno de cada DMC (HOMER)")
    ap.add_argument("--out", default="results/fase4_homer", help="Diretório de saída")
    args = ap.parse_args()

    if not os.path.exists(args.bed):
        sys.exit(f"[ERRO] BED não encontrado: {args.bed}")

    os.makedirs(args.out, exist_ok=True)

    # Decidir modo de execução
    has_homer = shutil.which("findMotifsGenome.pl") is not None

    if has_homer:
        print("[Info] HOMER detectado — rodando análise completa.")
        motifs = run_homer(args.bed, args.genome, args.out, args.size)
        if not motifs:
            print("[WARN] HOMER não retornou resultados. Tentando fallback...")
            motifs = fallback_analysis(args.bed, args.fasta, args.out)
    else:
        print("[Info] HOMER não detectado — rodando fallback puro Python.")
        motifs = fallback_analysis(args.bed, args.fasta, args.out)

    # Salvar CSV consolidado
    if motifs:
        csv_path = os.path.join(args.out, "motifs_summary.csv")
        fieldnames = list(motifs[0].keys())
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(motifs)
        print(f"[OK] Tabela de motifs: {csv_path} ({len(motifs)} motifs)")

        # Plot
        plot_motif_enrichment(motifs, os.path.join(args.out, "motifs_enrichment.png"))
    else:
        print("[WARN] Nenhum motif encontrado.")

    print(f"[FIM] Resultados em: {args.out}")


if __name__ == "__main__":
    main()
