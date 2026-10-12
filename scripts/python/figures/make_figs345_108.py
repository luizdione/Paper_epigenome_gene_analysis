#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_figs345_108.py — Figuras 3, 4, 5 do painel 108 (English labels, PNG+PDF).

Lê apenas artefatos versionados de results/SUB02_runs_108/:
  Fig 3 (multi-omic bridge)      <- _bridge_multiomic/bridge_gene_table.tsv
  Fig 4 (module enrichment)      <- _bridge_multiomic/bridge_module_enrichment.tsv (N=108, n=19)
  Fig 5 (fetal-liver validation) <- C6_.../c6_validation_table.tsv (+ RRBS Δβ do bridge)
Saída: <FIG_OUT>/Fig0{3,4,5}_*.{png,pdf} + Fig0{3,4,5}.png. Idempotente.
Tamanhos de tela/letra vêm do perfil de estilo em fig_paths.py (FIG_PROFILE=default|print).
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fig_paths import SUB02 as R, FIG_OUT as FIG, save_figure, P, wrap, set_title, layout  # noqa: E402

GREEN = "#2E9E54"; GREY = "#9AA0A6"; BLUE = "#2C6FBB"; PURPLE = "#7B5EA7"
RED = "#C0392B"; GOLD = "#C9971B"; INK = "#1F1F1F"

HIGH_TIER = {"GHSR", "HSPA2", "HSPB1", "PRLR", "ITGB3",
             "ERN1", "GRHL3", "PADI4", "LGALS3", "ISG15"}

MOD_LABEL = {
    "thermotolerance": "Thermotolerance (HSP)",
    "UPR_heat_stress": "Unfolded-protein response",
    "isg_panel": "Interferon-stimulated (ISG)",
    "tick_resistance": "Tick resistance",
    "milk_production": "Milk production",
    "thermo_reproduction": "Thermo-reproduction",
    "thermotolerance_extra": "Thermotolerance (extra)",
    "uterine_receptivity": "Uterine receptivity",
    "weight_gain": "Somatotropic (weight gain)",
    "embryo_development": "Embryo development",
    "oocyte_quality": "Oocyte quality",
    "galectin": "Galectin",
    "housekeeping": "Housekeeping",
    "implantation": "Implantation",
    "trophoblast_adhesion": "Trophoblast adhesion",
    "weight_gain_reproduction": "Weight gain / reproduction",
    "qc_epithelial_marker": "QC epithelial marker",
    "qc_neg_control": "QC negative control",
}
THERMAL = {"thermotolerance", "UPR_heat_stress", "thermotolerance_extra"}


def rd(path, delim="\t"):
    with open(path, encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter=delim))


def fnum(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def T(v):
    return str(v).strip().lower() in ("1", "true", "yes", "y", "sim")


def save(fig, stem):
    save_figure(fig, stem)   # grava Fig0N_<desc>.{png,pdf} e Fig0N.png em FIG_OUT
    plt.close(fig)


# --------------------------------------------------------------- Fig 3: bridge
def fig3():
    S = P("fig3")
    rows = rd(os.path.join(R, "_bridge_multiomic", "bridge_gene_table.tsv"))
    both = [r for r in rows if T(r.get("in_methylation_SUB02")) and T(r.get("in_expression_sub03"))]
    fig, ax = plt.subplots(figsize=S["figsize"], dpi=300)
    ax.axhline(0, color="#CCCCCC", lw=0.8, zorder=1)
    ax.axvline(0, color="#CCCCCC", lw=0.8, zorder=1)
    pts = {}
    for r in both:
        x = fnum(r.get("reg_methdiff")) or fnum(r.get("mean_methdiff")) or 0.0
        y = fnum(r.get("expr_best_log2fc")) or 0.0
        score = fnum(r.get("priority_score"), 1) or 1
        conc = r.get("concordance") == "concordant_repressive"
        prom = T(r.get("has_promoter_dmc"))
        col = GREEN if conc else GREY
        ax.scatter(x, y, s=28 + 34 * score, c=col, alpha=0.85,
                   edgecolors=INK if prom else "none", linewidths=1.0, zorder=3)
        pts[r["gene"]] = (x, y)
    # Manual label placement: the central cluster (ISG15/HSPA2/LGALS3/ERN1/HSPB1/PRLR)
    # collides under a single fixed offset. (dx, dy) in points (scaled by S["off"] in the
    # print profile, where the labels are bigger); arrow=True adds a thin leader to the
    # marker; a white halo keeps every label legible over any circle.
    LAB = {
        "GHSR":   (-3,  13, "center", "bottom", False),
        "GRHL3":  (7,    4, "left",   "bottom", False),
        "ITGB3":  (8,    1, "left",   "center", False),
        "PADI4":  (8,    1, "left",   "center", False),
        "ISG15":  (-11,  8, "right",  "center", True),
        "HSPA2":  (-11, -9, "right",  "center", True),
        "LGALS3": (-5, -15, "center", "top",    True),
        "ERN1":   (1,   15, "center", "bottom", True),
        "HSPB1":  (7,  -14, "center", "top",    True),
        "PRLR":   (12,   3, "left",   "center", True),
    }
    k = S["off"]
    for g, (dx, dy, ha, va, arrow) in LAB.items():
        if g not in pts:
            continue
        gx, gy = pts[g]
        ann = ax.annotate(
            g, (gx, gy), fontsize=S["label_fs"], fontweight="bold", fontstyle="italic",
            xytext=(dx * k, dy * k), textcoords="offset points",
            ha=ha, va=va, zorder=5,
            arrowprops=dict(arrowstyle="-", lw=0.6, color="#666666",
                            shrinkA=0.0, shrinkB=3.0) if arrow else None,
        )
        ann.set_path_effects([pe.withStroke(linewidth=S["halo"], foreground="white")])
    ax.set_xlabel("Regulatory methylation Δβ (%)  ·  hypo ← taurine | indicine → hyper", fontsize=S["axis"])
    ax.set_ylabel("Heat-stress expression (best log2 fold-change)", fontsize=S["axis"])
    top = set_title(fig, ax, "Multi-omic bridge: methylation Δβ vs heat-stress expression\n"
                             "47 genes marked at both layers · 19 concordantly repressive (108-gene panel)", S)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=GREEN, markersize=S["ms"],
                  label="concordant (repressive)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=GREY, markersize=S["ms"],
                  label="discordant"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="w", markeredgecolor=INK,
                  markersize=S["ms"], label="carries a promoter DMC"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=GREY, markersize=S["ms_big"],
                  label="size ∝ priority score")]
    ax.legend(handles=leg, loc=S["legend_loc"], fontsize=S["legend"], frameon=True)
    layout(fig, top)
    save(fig, "Fig03_bridge_concordance")


# --------------------------------------------------------------- Fig 4: modules
def fig4():
    S = P("fig4")
    rows = rd(os.path.join(R, "_bridge_multiomic", "bridge_module_enrichment.tsv"))
    rows = [r for r in rows if int(r["concordant"]) > 0]
    rows.sort(key=lambda r: float(r["hypergeom_p"]))
    labels, k, K, ps, cols = [], [], [], [], []
    for r in rows:
        mod = r["module"]
        labels.append(MOD_LABEL.get(mod, mod))
        k.append(int(r["concordant"])); K.append(int(r["genes_in_panel"]))
        ps.append(float(r["hypergeom_p"]))
        cols.append(GREEN if mod in THERMAL else RED if mod == "isg_panel"
                    else PURPLE if mod == "tick_resistance" else GREY)
    y = list(range(len(labels)))[::-1]
    fig, ax = plt.subplots(figsize=(S["width"], S["row"] * len(labels) + S["base"]), dpi=300)
    ax.barh(y, K, color="#E3E7EA", zorder=2, label="genes in panel")
    ax.barh(y, k, color=cols, zorder=3)
    for yi, kk, KK, p in zip(y, k, K, ps):
        ax.text(KK + 0.25, yi, "%d/%d  (p = %.2f)" % (kk, KK, p), va="center", fontsize=S["annot"], color=INK)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=S["tick"])
    ax.set_xlabel("Panel genes  (dark = concordant at both omic layers)", fontsize=S["axis"])
    ax.set_xlim(0, max(K) * S["xlim_factor"])
    top = set_title(fig, ax, "Within-panel module representation among concordant genes\n"
                             "(hypergeometric, universe = 108-gene panel, n = 19; no module survives BH)", S)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    leg = [Patch(fc=GREEN, label="thermal / heat-shock / UPR"),
           Patch(fc=RED, label="interferon (ISG)"),
           Patch(fc=PURPLE, label="tick resistance"),
           Patch(fc=GREY, label="other module")]
    lkw = {"bbox_to_anchor": S["legend_anchor"]} if S["legend_anchor"] else {}
    ax.legend(handles=leg, loc=S["legend_loc"], fontsize=S["legend"], frameon=True, **lkw)
    layout(fig, top)
    save(fig, "Fig04_module_enrichment")


# --------------------------------------------------------------- Fig 5: fetal liver
def fig5():
    S = P("fig5")
    c6 = rd(os.path.join(R, "C6_GIGAS_BIBI_vs_BTBT_fetalLiver", "c6_validation_table.tsv"))
    bt = {r["gene"]: r for r in rd(os.path.join(R, "_bridge_multiomic", "bridge_gene_table.tsv"))}
    # ordena: high-tier primeiro, depois com DMR, depois demais
    def key(r):
        return (0 if r["gene"] in HIGH_TIER else 1,
                0 if r["C6_fetalliver_WGBS"] == "yes" else 1, r["gene"])
    c6.sort(key=key, reverse=True)  # invert p/ barh (topo = high-tier)
    genes = [r["gene"] for r in c6]
    rrbs = [fnum((bt.get(g, {}) or {}).get("reg_methdiff")) or
            fnum((bt.get(g, {}) or {}).get("mean_methdiff")) or 0.0 for g in genes]
    c6d = [fnum(r["C6_reg_methdiff"]) if r["C6_fetalliver_WGBS"] == "yes" else None for r in c6]
    y = list(range(len(genes)))
    fig, ax = plt.subplots(figsize=(S["width"], S["row"] * len(genes) + S["base"]), dpi=300)
    h = 0.38
    ax.barh([yi + h / 2 for yi in y], rrbs, height=h, color=BLUE, label="RRBS blood/immune (SUB02)", zorder=3)
    ax.barh([yi - h / 2 for yi in y], [v if v is not None else 0 for v in c6d],
            height=h, color=PURPLE, label="WGBS fetal liver (C6)", zorder=3)
    ax.axvline(0, color="#888888", lw=0.8, zorder=2)
    labels = [("★ " + g) if g in HIGH_TIER else g for g in genes]
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=S["tick"], fontstyle="italic")
    ax.set_ylim(-0.8, len(genes) - 0.2)
    ax.set_xlabel("Regulatory Δβ (%)  ·  hypo ← taurine | indicine → hyper", fontsize=S["axis"])
    top = set_title(fig, ax, "Orthogonal validation in fetal-liver WGBS (MacPhillamy et al. 2024)\n"
                             "25 of 47 bridge candidates carry a fetal-liver DMR · 15 in the same direction", S)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if S["legend_outside"]:
        # perfil print: 47 linhas apertadas -> a legenda maior cobriria barras; vai para baixo do gráfico
        h_, l_ = ax.get_legend_handles_labels()
        fig.legend(h_, l_, loc="lower center", ncol=2, fontsize=S["legend"], frameon=True,
                   bbox_to_anchor=(0.5, 0.0), borderaxespad=0.4)
        layout(fig, top, 0.05)
    else:
        ax.legend(loc="lower right", fontsize=S["legend"], frameon=True)
        layout(fig, top)
    save(fig, "Fig05_C6_validation")


if __name__ == "__main__":
    fig3(); fig4(); fig5()
