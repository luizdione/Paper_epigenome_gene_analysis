#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_fig06_atlas_108.py — Fig. 6 (matched-sample validation, Zhou et al. 2020 atlas).

Re-plota a figura a partir do resultado já calculado
  results/SUB02_runs_108/_matched_atlas/matched_correlation.tsv  (88 genes; colunas:
  gene, module, n_tissues, spearman_r, p_value, direction, high_tier)
sem re-executar a análise. Diferenças em relação à versão submetida: 300 dpi, símbolos
de gene em itálico e as linhas verticais marcam os DEZ candidatos high-tier do painel 108
(Tabela 2), não a lista antiga de 8 genes gravada na coluna `high_tier`.
Saída: <FIG_OUT>/Fig06_matched_atlas.{png,pdf} + Fig06.png (via fig_paths).
Tamanhos de tela/letra vêm do perfil de estilo em fig_paths.py (FIG_PROFILE=default|print).
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fig_paths import SUB02, save_figure, P, wrap, set_title, layout  # noqa: E402

HIGH_TIER = ["GHSR", "HSPA2", "HSPB1", "PRLR", "ITGB3",
             "ERN1", "GRHL3", "PADI4", "LGALS3", "ISG15"]
GREEN = "#1e8449"; GREEN_D = "#1e6e3a"; GREY = "#9AA0A6"; INK = "#1F1F1F"; RED = "#C0392B"


def main():
    S = P("fig6")
    path = os.path.join(SUB02, "_matched_atlas", "matched_correlation.tsv")
    with open(path, encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t") if r.get("spearman_r") not in ("", "NA", None)]
    r = [float(x["spearman_r"]) for x in rows]
    n = len(r)
    neg = sum(1 for v in r if v < 0)
    by = {x["gene"]: x for x in rows}
    # asserções contra os números do relatório (lidos do mesmo TSV)
    assert n == 88, n
    assert neg == 61, neg
    assert abs(float(by["ISG15"]["spearman_r"]) + 0.673) < 1e-6, by["ISG15"]["spearman_r"]

    fig, ax = plt.subplots(figsize=S["figsize"], dpi=300)
    ax.axvspan(-1, 0, color="#E8F5E9", zorder=0)
    ax.hist(r, bins=16, range=(-1, 1), color=GREY, edgecolor="white", zorder=2)
    ax.axvline(0, color="#777777", lw=0.9, zorder=3)
    ymax = ax.get_ylim()[1] * 1.18
    ax.set_ylim(0, ymax)
    # candidatos high-tier: linha vertical + símbolo em itálico (deslocamento alternado p/ legibilidade)
    ht = sorted((g for g in HIGH_TIER if g in by), key=lambda g: float(by[g]["spearman_r"]))
    for i, g in enumerate(ht):
        v = float(by[g]["spearman_r"])
        ax.axvline(v, color=GREEN, lw=1.3, zorder=4)
        ax.text(v, ymax * (0.97 - S["stagger"] * (i % 3)), g, rotation=90, fontsize=S["label"],
                fontstyle="italic", fontweight="bold", va="top", ha="right", color=GREEN_D, zorder=5)
    ax.text(-0.98, ymax * 0.60, wrap("%d of %d genes\nnegative (repressive)" % (neg, n), S["annot_wrap"]),
            fontsize=S["annot"], color=INK, va="top", ha="left")
    ax.set_xlabel(wrap("Spearman r, promoter methylation vs expression (10 tissues, one Holstein cow)",
                       S["xlabel_wrap"]), fontsize=S["axis"])
    ax.set_ylabel("Panel genes", fontsize=S["axis"])
    top = set_title(fig, ax, "Matched-sample validation in a taurine WGBS+RNA-seq atlas (Zhou et al. 2020)\n"
                             "panel-level repressive direction; no single correlation survives Benjamini–Hochberg", S)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    layout(fig, top)
    save_figure(fig, "Fig06_matched_atlas")
    plt.close(fig)


if __name__ == "__main__":
    main()
