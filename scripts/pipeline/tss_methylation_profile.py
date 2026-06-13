#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fase 3 — Perfis de Metilação em torno do TSS (substituto deepTools)
===================================================================
Calcula, para cada gene em `all_target_genes_coordinates.json`, o sinal médio
de metilação em janelas de ±N kb centradas no Transcription Start Site (TSS),
respeitando a fita (+ ou −). Produz:

    * matriz.npy / matriz.csv  — matriz genes × bins (base para heatmap)
    * profile_<grupo>.png      — perfil médio por grupo (ex.: Angus vs Nelore)
    * heatmap_<grupo>.png      — heatmap de todas as regiões

Entradas:
    --methylation : arquivo(s) por amostra, formato bismark .cov.gz:
                    chr  start  end  methylation_pct  numC  numT
                    (colunas: 0=chrom, 1=pos, 3=meth%, 4=numC, 5=numT)
    --breeds      : lista CSV "Angus,Angus,Nelore,Nelore" alinhada aos arquivos
    --json        : JSON com coordenadas dos genes (padrão ARS-UCD1.3)
    --window      : tamanho da janela em bp (padrão 10000; alinhado à janela de análise da Fase 1 ±10 kb)
    --bins        : número de bins (padrão 40 = resolução de 100 bp/bin)

Uso:
    python tss_methylation_profile.py \
      --methylation sample1.tsv sample2.tsv sample3.tsv sample4.tsv \
      --breeds Angus,Angus,Nelore,Nelore \
      --json all_target_genes_coordinates.json \
      --out results/fase3_deeptools

Dependências: numpy, matplotlib (instalar com --break-system-packages se preciso)
Autor: Luiz Dione Barbosa De-Melo — GenOmics Epigenoma 2026
"""
from __future__ import annotations
import argparse
import csv
import gzip
import json
import os
import sys
from collections import defaultdict
from typing import Dict, List, Tuple

try:
    import numpy as np
except ImportError:
    sys.exit("[ERRO] numpy é requerido: pip install numpy --break-system-packages")


def _open(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")


def norm_chrom(c) -> str:
    return str(c).strip().lower().removeprefix("chr")


# --------------------------------------------------------------------- IO
def load_targets(json_path: str) -> List[dict]:
    """Retorna lista de genes com TSS absoluto e metadados."""
    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    genes = []
    for name, info in data.get("genes", {}).items():
        raw_strand = info.get("strand", 1)
        # Aceita: 1, -1, "+", "-", "1", "-1"
        if isinstance(raw_strand, str):
            raw_strand = raw_strand.strip()
            if raw_strand in ("+", "1"):
                strand = 1
            elif raw_strand in ("-", "-1"):
                strand = -1
            else:
                strand = 1   # fallback
        else:
            strand = int(raw_strand) if raw_strand else 1
        start, end = info.get("start"), info.get("end")
        if start is None or end is None:
            continue
        tss = start if strand == 1 else end
        genes.append(
            {
                "gene": name,
                "chrom": norm_chrom(info.get("chromosome")),
                "tss": int(tss),
                "strand": strand,
                "module": info.get("module", ""),
            }
        )
    return genes


def read_sample(path: str) -> Dict[str, List[Tuple[int, float]]]:
    """Lê 1 arquivo methylation: retorna dict chrom -> [(pos, freqC%), ...]"""
    out: Dict[str, List[Tuple[int, float]]] = defaultdict(list)
    with _open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        first = next(reader, None)
        if first is None:
            return out
        # tenta detectar se é header
        try:
            _ = int(first[1])
            rows = [first] + [r for r in reader]
        except (ValueError, IndexError):
            rows = list(reader)
        for r in rows:
            if len(r) < 4:
                continue
            try:
                chrom = norm_chrom(r[0])
                pos = int(r[1])
                # Bismark .cov.gz: chr start end meth_pct numC numT
                # FIX: usar col[3] como % metilação; col[4]+col[5] como cobertura
                if len(r) >= 6:
                    cov = int(float(r[4])) + int(float(r[5]))
                    freqC = float(r[3])        # coluna 3 = methylation %
                else:
                    # fallback methylKit: chr pos strand coverage numC numT
                    cov = int(float(r[3]))
                    freqC = (int(float(r[4])) / cov * 100) if cov > 0 else 0.0
                if cov <= 0:
                    continue
            except (ValueError, IndexError):
                continue
            out[chrom].append((pos, freqC))
    # ordena por posição para busca binária
    for c in out:
        out[c].sort()
    return out


# --------------------------------------------------------------------- core
def compute_matrix(
    sample_data: Dict[str, List[Tuple[int, float]]],
    genes: List[dict],
    window: int,
    bins: int,
) -> np.ndarray:
    """Retorna matriz (n_genes, bins) com metilação média por bin."""
    bin_size = (2 * window) / bins
    mat = np.full((len(genes), bins), np.nan, dtype=np.float32)
    for gi, g in enumerate(genes):
        cpgs = sample_data.get(g["chrom"], [])
        if not cpgs:
            continue
        start = g["tss"] - window
        end = g["tss"] + window
        # busca linear simples (cpgs é ordenada)
        buckets = [[] for _ in range(bins)]
        for pos, val in cpgs:
            if pos < start:
                continue
            if pos > end:
                break
            rel = pos - start
            if g["strand"] == -1:
                rel = (2 * window) - rel  # espelha para TSS orientado
            bi = min(int(rel / bin_size), bins - 1)
            buckets[bi].append(val)
        for bi, vals in enumerate(buckets):
            if vals:
                mat[gi, bi] = float(np.mean(vals))
    return mat


# --------------------------------------------------------------------- plotting
def plot_profile(matrices_by_group: Dict[str, np.ndarray], out: str, window: int):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib não disponível — pulando plots.")
        return
    plt.figure(figsize=(8, 5))
    x = np.linspace(-window, window, next(iter(matrices_by_group.values())).shape[1])
    for grp, m in matrices_by_group.items():
        mean = np.nanmean(m, axis=0)
        sem = np.nanstd(m, axis=0) / np.sqrt(np.sum(~np.isnan(m), axis=0))
        plt.plot(x, mean, label=grp, linewidth=2)
        plt.fill_between(x, mean - sem, mean + sem, alpha=0.2)
    plt.axvline(0, color="k", linestyle="--", linewidth=0.7)
    plt.xlabel("Distância do TSS (bp)")
    plt.ylabel("Metilação média (%)")
    plt.title("Perfil de metilação em torno do TSS")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()


def plot_heatmap(mat: np.ndarray, genes: List[dict], out: str, window: int, title: str,
                 max_genes: int = 80):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    # FIX: filtrar genes com dados e limitar a max_genes (evita figura de milhares de pol)
    has_data = ~np.all(np.isnan(mat), axis=1)
    mat_filt = mat[has_data]
    genes_filt = [g for g, h in zip(genes, has_data) if h]

    if len(genes_filt) == 0:
        print(f"[WARN] Nenhum gene com dados — heatmap {out} não gerado.")
        return

    # Ordenar por média de metilação decrescente e limitar
    if len(genes_filt) > max_genes:
        means = np.nanmean(mat_filt, axis=1)
        order = np.argsort(-means)[:max_genes]
        mat_filt = mat_filt[order]
        genes_filt = [genes_filt[i] for i in order]
        print(f"[Info] Heatmap: exibindo top {max_genes}/{sum(has_data)} genes com dados")

    n = len(genes_filt)
    height = min(max(4, n * 0.22 + 2), 28)   # máx 28 pol de altura
    fig, ax = plt.subplots(figsize=(8, height))
    im = ax.imshow(mat_filt, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=100,
                   extent=(-window, window, 0, n))
    ax.set_yticks(np.arange(n) + 0.5)
    ax.set_yticklabels([g["gene"] for g in genes_filt][::-1], fontsize=7)
    ax.set_xlabel("Distância do TSS (bp)")
    ax.set_title(title)
    plt.colorbar(im, ax=ax, label="Metilação (%)")
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()


# --------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--methylation", nargs="+", required=True)
    ap.add_argument("--breeds", required=True, help="Lista CSV: Angus,Angus,Nelore,Nelore")
    ap.add_argument("--json", default="all_target_genes_coordinates.json")
    ap.add_argument("--window", type=int, default=10000)
    ap.add_argument("--bins", type=int, default=40)
    ap.add_argument("--out", default="results/fase3_deeptools")
    args = ap.parse_args()

    breeds = [b.strip() for b in args.breeds.split(",")]
    if len(breeds) != len(args.methylation):
        sys.exit(f"[ERRO] Nº de breeds ({len(breeds)}) ≠ nº de arquivos ({len(args.methylation)})")

    os.makedirs(args.out, exist_ok=True)
    genes = load_targets(args.json)
    print(f"[Info] Genes-alvo carregados: {len(genes)}")

    # agrupa amostras por grupo e calcula matriz média
    matrices_by_group: Dict[str, List[np.ndarray]] = defaultdict(list)
    for path, breed in zip(args.methylation, breeds):
        print(f"[Info] Lendo {path} → grupo {breed}")
        sd = read_sample(path)
        m = compute_matrix(sd, genes, args.window, args.bins)
        matrices_by_group[breed].append(m)

    avg_by_group = {grp: np.nanmean(np.stack(ms), axis=0) for grp, ms in matrices_by_group.items()}

    # outputs
    for grp, m in avg_by_group.items():
        np.savetxt(os.path.join(args.out, f"matrix_{grp}.csv"), m, delimiter=",",
                   header=",".join(f"bin{i}" for i in range(m.shape[1])), comments="")
        plot_heatmap(m, genes, os.path.join(args.out, f"heatmap_{grp}.png"),
                     args.window, f"Metilação TSS ±{args.window}bp — {grp}")

    plot_profile(avg_by_group, os.path.join(args.out, "profile_comparison.png"), args.window)

    # gene table — inclui média de metilação por grupo (nanmean sobre todos os bins do perfil)
    group_names = sorted(avg_by_group.keys())
    with open(os.path.join(args.out, "genes_ordered.csv"), "w", encoding="utf-8") as fh:
        header_cols = ["gene", "chromosome", "tss", "strand", "module"]
        for grp in group_names:
            header_cols.append(f"{grp}_mean_pct")
        if len(group_names) == 2:
            header_cols.append("delta_pct")   # grupo1 - grupo0 (taurus - indicus)
        fh.write(",".join(header_cols) + "\n")
        for i, g in enumerate(genes):
            strand_str = "+" if g['strand'] in (1, "1", "+") else "-" if g['strand'] in (-1, "-1", "-") else str(g['strand'])
            row = [g['gene'], g['chrom'], str(g['tss']), strand_str, g['module']]
            grp_means = []
            for grp in group_names:
                col_vals = avg_by_group[grp][i, :]  # vetor de bins para este gene
                mean_val = float(np.nanmean(col_vals)) if not np.all(np.isnan(col_vals)) else 0.0
                row.append(f"{mean_val:.4f}")
                grp_means.append(mean_val)
            if len(group_names) == 2:
                row.append(f"{grp_means[1] - grp_means[0]:.4f}")
            fh.write(",".join(row) + "\n")

    print(f"[OK] Resultados salvos em: {args.out}")


if __name__ == "__main__":
    main()
