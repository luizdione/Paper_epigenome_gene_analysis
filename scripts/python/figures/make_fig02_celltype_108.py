#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_fig02_celltype_108.py — Figura 2 (painel 108 genes / +-15 kb).

DMCs por contraste de sistematizacao nos dois limiares (|db|>15% e |db|>25%),
destacando o gradiente subespecie-sobre-estacao e o gradiente linfoide>mieloide.
Numeros lidos de results/SUB02_runs_108/AGGREGATE_108.json (fonte-verdade).
Saida: <FIG_OUT>/Fig02_celltype_dmc.{png,pdf} + Fig02.png.
Tamanhos de tela/letra vem do perfil de estilo em fig_paths.py (FIG_PROFILE=default|print);
no perfil "print" a legenda sai para baixo do grafico (nao cabe mais dentro dele).
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fig_paths import SUB02, save_figure, P, wrap, set_title, layout  # noqa: E402
AGG = os.path.join(SUB02, "AGGREGATE_108.json")
S = P("fig2")

# rotulo legivel + eixo biologico por contraste
LABELS = {
    "C1_zenodo_NEL_vs_HF":          ("C1  Nellore vs Holstein (immune pool)", "constitutive"),
    "C7_bonus_zenodo_NEL_vs_NDama": ("C7  Nellore vs N'Dama (whole blood)",   "constitutive"),
    "C2_GSE_NelFev_vs_AngFev_HS":   ("C2  Nellore vs Angus (February)",       "season"),
    "C3_GSE_NelJun_vs_AngJun":      ("C3  Nellore vs Angus (June)",           "season"),
    "C4_GSE_NelJun_vs_NelFev_HS":   ("C4  Nellore June vs February",          "plastic"),
    "C5_GSE_AngJun_vs_AngFev_HS":   ("C5  Angus June vs February",           "plastic"),
    "C8_zenodo_NEL_vs_HF_CD4":      ("C8  CD4 T",   "lymphoid"),
    "C9_zenodo_NEL_vs_HF_CD8":      ("C9  CD8 T",   "lymphoid"),
    "C10_zenodo_NEL_vs_HF_Bcell":   ("C10 B cell",  "lymphoid"),
    "C11_zenodo_NEL_vs_HF_gdT":     ("C11 gamma-delta T", "lymphoid"),
    "C12_zenodo_NEL_vs_HF_gran":    ("C12 granulocyte", "myeloid"),
    "C13_zenodo_NEL_vs_HF_mono":    ("C13 monocyte",    "myeloid"),
    "C14_zenodo_NEL_vs_HF_NK":      ("C14 NK cell",     "lymphoid"),
}
# ordem de exibicao (topo -> base): constitutivo, estacao, plastico, tipos celulares
ORDER = ["C1_zenodo_NEL_vs_HF", "C7_bonus_zenodo_NEL_vs_NDama",
         "C2_GSE_NelFev_vs_AngFev_HS", "C3_GSE_NelJun_vs_AngJun",
         "C4_GSE_NelJun_vs_NelFev_HS", "C5_GSE_AngJun_vs_AngFev_HS",
         "C8_zenodo_NEL_vs_HF_CD4", "C9_zenodo_NEL_vs_HF_CD8",
         "C10_zenodo_NEL_vs_HF_Bcell", "C11_zenodo_NEL_vs_HF_gdT",
         "C12_zenodo_NEL_vs_HF_gran", "C13_zenodo_NEL_vs_HF_mono",
         "C14_zenodo_NEL_vs_HF_NK"]

COL = {"constitutive": "#2C6FBB", "season": "#7B5EA7", "plastic": "#B0B0B0",
       "lymphoid": "#2E9E54", "myeloid": "#BFE6CC"}


def main():
    agg = json.load(open(AGG, encoding="utf-8"))
    pc = agg["per_contrast"]

    labels, d15, d25, cols = [], [], [], []
    for k in ORDER:
        lab, axis = LABELS[k]
        labels.append(lab)
        d15.append(pc[k]["dmc_diff15"])
        d25.append(pc[k]["dmc_diff25"])
        cols.append(COL[axis])

    y = list(range(len(labels)))[::-1]  # topo = primeiro
    fig, ax = plt.subplots(figsize=S["figsize"], dpi=300)
    # barra clara = |db|>15% ; barra escura sobreposta = |db|>25%
    ax.barh(y, d15, color=cols, edgecolor="none", height=0.70, alpha=0.55, zorder=2)
    ax.barh(y, d25, color=cols, edgecolor="none", height=0.42, alpha=1.0, zorder=3)
    for yi, v15, v25 in zip(y, d15, d25):
        ax.text(v15 + 8, yi, "%d / %d" % (v15, v25), va="center", ha="left",
                fontsize=S["annot"], color="#1F1F1F")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=S["tick"])
    ax.set_xlabel(wrap("Differentially methylated cytosines (|Δβ| > 15%  /  > 25%)", S["xlabel_wrap"]),
                  fontsize=S["axis"])
    ax.set_xlim(0, max(d15) * S["xlim_factor"])
    top = set_title(fig, ax, "DMCs per systematisation contrast (108-gene panel, ±15 kb)\n"
                             "constitutive subspecies divergence dominates; lymphoid > myeloid", S)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    legend = [Patch(fc=COL["constitutive"], label="subspecies (blood / immune pool)"),
              Patch(fc=COL["season"], label="subspecies by season"),
              Patch(fc=COL["plastic"], label="within-breed seasonal (plastic)"),
              Patch(fc=COL["lymphoid"], label="lymphoid cell type"),
              Patch(fc=COL["myeloid"], label="myeloid cell type")]
    if S["legend_outside"]:
        fig.legend(handles=legend, loc="lower center", ncol=2, fontsize=S["legend"],
                   frameon=True, bbox_to_anchor=(0.5, 0.0), borderaxespad=0.4)
        layout(fig, top, 0.13)
    else:
        ax.legend(handles=legend, loc="lower right", fontsize=S["legend"], frameon=True)
        layout(fig, top)
    save_figure(fig, "Fig02_celltype_dmc")   # Fig02_celltype_dmc.{png,pdf} + Fig02.png em FIG_OUT
    plt.close(fig)


if __name__ == "__main__":
    main()
