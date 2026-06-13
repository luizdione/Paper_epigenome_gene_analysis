#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Monta as figuras do artigo UNIFICADO em PAPER_UNIFICADO/figuras/:
  Fig01_discovery_dmc.png      <- SUB01 run04 (distribuição funcional das DMCs)
  Fig02_celltype_dmc.png       <- NOVO: DMCs por tipo celular (lido dos resultados)
  Fig03_bridge_concordance.png <- ponte re-rodada (_bridge_multiomic)
  Fig04_module_enrichment.png  <- enriquecimento de módulos (_bridge_multiomic)
  Fig05_C6_validation.png      <- validação WGBS re-rodada (PAPER_SUB02/figuras)
Nenhum número inserido à mão: as contagens do Fig02 são lidas dos CSVs reais.
"""
from pathlib import Path
import shutil
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[2]
RUNS = BASE / "results" / "SUB02_runs"
BRIDGE = RUNS / "_bridge_multiomic"
FIGOUT = BASE / "PAPER_UNIFICADO" / "figuras"
FIGOUT.mkdir(parents=True, exist_ok=True)

CELLTYPES = [
    ("C8_zenodo_NEL_vs_HF_CD4", "CD4 T"),
    ("C9_zenodo_NEL_vs_HF_CD8", "CD8 T"),
    ("C10_zenodo_NEL_vs_HF_Bcell", "B cell"),
    ("C11_zenodo_NEL_vs_HF_gdT", "γδ-T"),
    ("C12_zenodo_NEL_vs_HF_gran", "Granulocyte"),
    ("C13_zenodo_NEL_vs_HF_mono", "Monocyte"),
    ("C14_zenodo_NEL_vs_HF_NK", "NK"),
]


def count_dmcs(comp, diff):
    f = RUNS / comp / f"fase1_methylkit_diff{diff}" / "DMCs_significant.csv"
    if not f.exists():
        return 0
    with f.open(encoding="utf-8", errors="ignore") as fh:
        return max(sum(1 for _ in fh) - 1, 0)


def fig_celltype():
    labels, d15, d25 = [], [], []
    for comp, lab in CELLTYPES:
        labels.append(lab)
        d15.append(count_dmcs(comp, 15))
        d25.append(count_dmcs(comp, 25))
    # ordena por linfoide->mieloide ja na ordem dada; cor por linhagem
    lymphoid = {"CD4 T", "CD8 T", "B cell", "γδ-T", "NK"}
    colors = ["#7d3c98" if l in lymphoid else "#b9770e" for l in labels]
    x = range(len(labels))
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    ax.bar([i - 0.2 for i in x], d15, width=0.4, color=colors, alpha=0.95,
           label="|Δβ| > 15%")
    ax.bar([i + 0.2 for i in x], d25, width=0.4, color=colors, alpha=0.5,
           label="|Δβ| > 25%")
    for i, v in enumerate(d15):
        ax.text(i - 0.2, v + 6, str(v), ha="center", fontsize=8.5, fontweight="bold")
    ax.set_xticks(list(x)); ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Significant DMCs (Nelore vs Holstein)")
    ax.set_title("Cell-type-resolved subspecies methylation divergence\n"
                 "purple = lymphoid · orange = myeloid", fontsize=10.5)
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = FIGOUT / "Fig02_celltype_dmc.png"
    fig.savefig(out, dpi=200); plt.close(fig)
    print("[Fig02] DMCs/tipo celular (diff15):", d15, "->", out.name)


def copy_fig(src, dst):
    if Path(src).exists():
        shutil.copy(src, FIGOUT / dst)
        print(f"[copia] {dst}")
    else:
        print(f"[aviso] ausente: {src}")


def main():
    # Fig01 — descoberta (SUB01)
    sub01 = BASE / "results" / "SUB01_GSE161113_run04" / "fase2_intersect"
    for cand in ["dmc_functional_distribution.png", "dmc_functional_by_gene.png"]:
        if (sub01 / cand).exists():
            copy_fig(sub01 / cand, "Fig01_discovery_dmc.png"); break
    fig_celltype()
    copy_fig(BRIDGE / "bridge_concordance.png", "Fig03_bridge_concordance.png")
    copy_fig(BRIDGE / "bridge_modules.png", "Fig04_module_enrichment.png")
    copy_fig(BASE / "PAPER_SUB02" / "figuras" / "Fig04_C6_validation.png",
             "Fig05_C6_validation.png")
    print("Figuras em:", FIGOUT)


if __name__ == "__main__":
    main()
