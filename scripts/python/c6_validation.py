#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
c6_validation.py — C6 como VALIDAÇÃO ORTOGONAL (não na ponte primária):
cruza os genes com DMR de fígado fetal WGBS (Angus−Brahman; MacPhillamy et al.
2024, GigaScience giae061) com os candidatos da ponte multi-ômica (SUB02 RRBS +
sub03 RNA-seq). Confere se a DIREÇÃO da metilação raça é consistente entre
estudos/tecidos/ensaios (RRBS sangue/imune vs WGBS fígado fetal) — ou seja, se o
achado é robusto a tecido e a tipo de assay.

Saída: results/SUB02_runs/C6_.../c6_validation_table.tsv
       PAPER_SUB02/figuras/Fig04_C6_validation.png
       PAPER_SUB02/tabelas/TableS5_C6_validation.tsv
"""
from __future__ import annotations
import csv
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
C6 = BASE / "results" / "SUB02_runs" / "C6_GIGAS_BIBI_vs_BTBT_fetalLiver"
C6_ANN = C6 / "fase2_annotation_diff15" / "annotation_table.csv"
BRIDGE = BASE / "results" / "SUB02_runs" / "_bridge_multiomic" / "bridge_gene_table.tsv"
FIG = BASE / "PAPER_SUB02" / "figuras" / "Fig04_C6_validation.png"
TAB = BASE / "PAPER_SUB02" / "tabelas" / "TableS5_C6_validation.tsv"


def c6_by_gene():
    agg = defaultdict(lambda: {"md": [], "reg": [], "prom": 0, "n": 0})
    with C6_ANN.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            g = r["gene"].strip().strip('"')
            md = float(r["meth.diff"])
            feat = r["feature"].strip().strip('"')
            a = agg[g]
            a["n"] += 1
            a["md"].append(md)
            if feat in ("promoter", "upstream_distal"):
                a["reg"].append(md)
            if feat == "promoter":
                a["prom"] += 1
    out = {}
    for g, a in agg.items():
        mean = sum(a["md"]) / len(a["md"])
        reg = (sum(a["reg"]) / len(a["reg"])) if a["reg"] else mean
        out[g] = {"n_dmr": a["n"], "mean_md": round(mean, 1),
                  "reg_md": round(reg, 1),
                  "dir": "hyper" if reg > 0 else "hypo",
                  "has_prom": "yes" if a["prom"] else "no"}
    return out


def main():
    c6 = c6_by_gene()
    rows = list(csv.DictReader(BRIDGE.open(encoding="utf-8"), delimiter="\t"))
    both = [r for r in rows if r["in_methylation_SUB02"] == "yes"
            and r["in_expression_sub03"] == "yes"]

    out_rows = []
    corrob = 0
    dir_agree = 0
    for r in both:
        g = r["gene"]
        c = c6.get(g)
        in_c6 = "yes" if c else "no"
        agree = ""
        if c:
            corrob += 1
            agree = "yes" if c["dir"] == r["meth_direction"] else "no"
            if agree == "yes":
                dir_agree += 1
        out_rows.append({
            "gene": g,
            "priority_tier": r["priority_tier"],
            "sub02_meth_dir": r["meth_direction"],
            "expr_dir": r["expr_direction"],
            "concordance": r["concordance"],
            "C6_fetalliver_WGBS": in_c6,
            "C6_n_DMR": c["n_dmr"] if c else "",
            "C6_reg_methdiff": c["reg_md"] if c else "",
            "C6_dir": c["dir"] if c else "",
            "C6_has_promoter": c["has_prom"] if c else "",
            "C6_direction_agrees": agree,
        })
    # ordena: prioridade alta primeiro, depois corroborados
    order = {"alta": 0, "média": 1, "baixa": 2, "—": 3}
    out_rows.sort(key=lambda x: (order.get(x["priority_tier"], 9),
                                 x["C6_fetalliver_WGBS"] != "yes", x["gene"]))
    TAB.parent.mkdir(parents=True, exist_ok=True)
    with TAB.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(out_rows)
    # cópia auditável junto ao run
    (C6 / "c6_validation_table.tsv").write_text(TAB.read_text(encoding="utf-8"),
                                                encoding="utf-8")

    # figura: candidatos da ponte (tier alta/média) com corroboração C6
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        focus = [r for r in out_rows if r["priority_tier"] in ("alta", "média")]
        focus = focus[::-1]
        labels = [r["gene"] + (" ★" if r["priority_tier"] == "alta" else "") for r in focus]
        y = np.arange(len(focus))
        sub02 = [float(rows_lookup(both, r["gene"], "reg_methdiff") or
                       rows_lookup(both, r["gene"], "mean_methdiff") or 0) for r in focus]
        c6v = [float(r["C6_reg_methdiff"]) if r["C6_reg_methdiff"] != "" else 0 for r in focus]
        fig, ax = plt.subplots(figsize=(8.6, 0.42 * len(focus) + 1.6))
        ax.barh(y + 0.2, sub02, height=0.38, color="#2c7fb8", label="SUB02 RRBS (sangue/imune)")
        ax.barh(y - 0.2, c6v, height=0.38, color="#8e44ad", label="C6 WGBS (fígado fetal)")
        ax.axvline(0, color="#555", lw=0.8)
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.5)
        ax.set_xlabel("Δβ regulatório  (hipo ← taurino | indicino → hiper)")
        ax.set_title("C6 — validação ortogonal: metilação raça em fígado fetal (WGBS)\n"
                     "vs candidatos da ponte (RRBS). Mesma direção = robusto a tecido/ensaio.",
                     fontsize=9.5)
        ax.legend(frameon=False, fontsize=8.5, loc="lower right")
        fig.tight_layout(); fig.savefig(FIG, dpi=200); plt.close(fig)
    except Exception as exc:
        print("[warn] figura:", exc)

    print(f"[C6] genes com DMR fígado fetal: {len(c6)}")
    print(f"[C6] candidatos da ponte (both-layer): {len(both)}; "
          f"corroborados por C6: {corrob}; com mesma direção: {dir_agree}")
    print(f"[C6] tabela: {TAB}")
    print(f"[C6] figura: {FIG}")
    return 0


def rows_lookup(rows, gene, col):
    for r in rows:
        if r["gene"] == gene:
            v = r.get(col, "")
            return v if v != "" else None
    return None


if __name__ == "__main__":
    raise SystemExit(main())
