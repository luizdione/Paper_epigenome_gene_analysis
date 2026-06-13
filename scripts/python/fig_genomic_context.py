#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig_genomic_context.py — figura de 3 painéis do contexto genômico e validação
em amostras pareadas:
  (A) cis-SNP (DMC-level)      <- _cissnp/cissnp_table.tsv
  (B) tri-camada (gene-level)  <- _triple_layer/triple_layer_table.tsv
  (C) atlas pareado            <- _matched_atlas/matched_correlation.tsv
Lê apenas resultados reais. Saída: figures/Fig06_genomic_context.png
"""
from __future__ import annotations
import csv
from collections import Counter
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[2]
RUNS = BASE / "results" / "SUB02_runs"
FIG = BASE / "PAPER_UNIFICADO" / "figuras" / "Fig06_genomic_context.png"
HIGH = {"PRLR","GHSR","HSPA2","HSPB1","ISG15","MX1","ITGB3","KRT18"}


def read_tsv(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    cis = read_tsv(RUNS / "_cissnp" / "cissnp_table.tsv")
    tri = read_tsv(RUNS / "_triple_layer" / "triple_layer_table.tsv")
    mp = RUNS / "_matched_atlas" / "matched_correlation.tsv"
    matched = read_tsv(mp) if mp.exists() else []

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.3))

    # (A) cis-SNP
    cc = Counter(r["class"] for r in cis)
    order = ["no-SNP", "near-SNP", "CpG-SNP"]
    vals = [cc.get(k, 0) for k in order]
    cols = ["#1e8449", "#e1a93a", "#c0392b"]
    n = sum(vals)
    ax[0].bar(order, vals, color=cols, edgecolor="white")
    for i, v in enumerate(vals):
        ax[0].text(i, v + n*0.01, f"{v}\n({100*v/n:.0f}%)", ha="center", fontsize=9)
    ax[0].set_title(f"A. cis-genetic context of DMCs (n={n})\nno-SNP = genuinely epigenetic",
                    fontsize=10)
    ax[0].set_ylabel("subspecies DMCs (positions)")
    ax[0].spines[["top", "right"]].set_visible(False)

    # (B) tri-camada
    tc = Counter(r["layer_class"] for r in tri)
    labels = ["triple\n(genomic+epi+expr)", "epi-dominant\n(no strong Fst)"]
    tvals = [tc.get("triple", 0), tc.get("epi-dominant", 0)]
    ax[1].bar(labels, tvals, color=["#7d3c98", "#2e86c1"], edgecolor="white")
    for i, v in enumerate(tvals):
        ax[1].text(i, v + 0.3, str(v), ha="center", fontsize=10, fontweight="bold")
    ht_epi = sum(1 for r in tri if r["high_tier"] == "yes" and r["layer_class"] == "epi-dominant")
    ht_tot = sum(1 for r in tri if r["high_tier"] == "yes")
    ax[1].set_title(f"B. Genomic differentiation of candidates\nall {ht_epi}/{ht_tot} high-tier are epi-dominant",
                    fontsize=10)
    ax[1].set_ylabel("bridge candidate genes")
    ax[1].spines[["top", "right"]].set_visible(False)

    # (C) atlas pareado
    if matched:
        rs = [float(r["spearman_r"]) for r in matched]
        neg = sum(1 for v in rs if v < 0)
        ax[2].hist(rs, bins=12, color="#9bbcd8", edgecolor="white")
        for r in matched:
            if r["high_tier"] == "yes":
                ax[2].axvline(float(r["spearman_r"]), color="#1e8449", lw=1.2)
        ax[2].axvline(0, color="#444", lw=1)
        ax[2].set_title(f"C. Matched-sample test (Holstein atlas)\n"
                        f"{neg}/{len(rs)} genes negative (repressive)", fontsize=10)
        ax[2].set_xlabel("Spearman r (meth vs expr, across tissues)")
        ax[2].set_ylabel("panel genes")
    else:
        ax[2].text(0.5, 0.5, "matched_correlation.tsv\nainda não gerado", ha="center",
                   va="center", transform=ax[2].transAxes)
        ax[2].set_title("C. Matched-sample test (pendente)", fontsize=10)
    ax[2].spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=200)
    print("figura salva:", FIG, "| matched rows:", len(matched))


if __name__ == "__main__":
    main()
