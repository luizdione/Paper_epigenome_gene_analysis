#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bridge_deg_dmr.py — Ponte multi-ômica: cruza os genes diferencialmente
METILADOS do epigenoma (SUB02, comparações raça×calor) com os genes
diferencialmente EXPRESSOS do genomas (sub03 RNA-seq heat stress) e com o
painel-alvo do sub04. Integração no nível do GENE (símbolo).

Saída: results/SUB02_runs/_bridge_multiomic/
    bridge_gene_table.tsv     uma linha por gene com sinais metilação×expressão
    bridge_summary.txt        contagens, concordância, genes-chave
    bridge_concordance.png     hipermet/hipomet × up/down regulado

Sem números inventados: tudo lido dos CSVs/TSVs reais. Genes sem dado de um
lado ficam com NA nesse lado.
"""
from __future__ import annotations

import csv
import glob
import os
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]               # raiz epigenoma
SUB02 = BASE / "results" / "SUB02_runs"
# saídas do genomas (montadas em E:/data; fallback p/ caminho do repo)
GENOMAS_DEG_DIRS = [
    Path("E:/data/results/sub03_deseq2"),
    BASE.parent / "genomas" / "results" / "sub03_rnaseq_heat_stress",
]
SUB04_PANEL = [
    Path("E:/data/results/sub04_reverse_aims/panel73_deg_status.tsv"),
]
OUT = SUB02 / "_bridge_multiomic"
OUT.mkdir(parents=True, exist_ok=True)

DEG_PADJ = 0.05


# Peso regulatório por feature (relevância para leitura repressiva da metilação).
# Promotor é o local canônico de repressão por metilação; corpo gênico associa-se
# frequentemente a metilação POSITIVAMENTE correlacionada à expressão.
FEATURE_WEIGHT = {"promoter": 1.0, "upstream_distal": 0.6,
                  "downstream_distal": 0.3, "gene_body": 0.2}
REGULATORY_FEATURES = ("promoter", "upstream_distal")  # usados p/ direção regulatória


# ----------------------------------------------------------------- SUB02 (metilação)
def load_sub02_methylation():
    """gene -> dict com sinais por feature/comparação, além do agregado."""
    genes = defaultdict(lambda: {"n_dmcs": 0, "methdiffs": [],
                                 "reg_methdiffs": [], "comparisons": defaultdict(int),
                                 "features": defaultdict(int)})
    for ann in sorted(SUB02.glob("C*/fase2_annotation_diff15/annotation_table.csv")):
        comp = ann.parents[1].name
        if comp.startswith("C6_"):
            continue  # C6 = WGBS/fígado fetal (estudo externo) — validação suplementar, não na ponte primária
        with ann.open(encoding="utf-8", errors="ignore") as fh:
            for row in csv.DictReader(fh):
                g = (row.get("gene") or "").strip().strip('"')
                if not g or g in ("NA", "."):
                    continue
                try:
                    md = float(row.get("meth.diff"))
                except (TypeError, ValueError):
                    continue
                feat = (row.get("feature") or "").strip().strip('"')
                rec = genes[g]
                rec["n_dmcs"] += 1
                rec["methdiffs"].append(md)
                rec["comparisons"][comp] += 1
                if feat:
                    rec["features"][feat] += 1
                if feat in REGULATORY_FEATURES:
                    rec["reg_methdiffs"].append(md)
    out = {}
    for g, r in genes.items():
        mean_md = sum(r["methdiffs"]) / len(r["methdiffs"]) if r["methdiffs"] else 0.0
        reg_md = (sum(r["reg_methdiffs"]) / len(r["reg_methdiffs"])
                  if r["reg_methdiffs"] else None)
        # direção regulatória (promotor/upstream) tem prioridade sobre o agregado
        dir_md = reg_md if reg_md is not None else mean_md
        reg_w = round(sum(FEATURE_WEIGHT.get(f, 0.2) * n
                          for f, n in r["features"].items()), 2)
        out[g] = {
            "n_dmcs": r["n_dmcs"],
            "mean_methdiff": round(mean_md, 2),
            "reg_methdiff": round(reg_md, 2) if reg_md is not None else "",
            "meth_direction": "hyper" if dir_md > 0 else "hypo",
            "has_promoter_dmc": "yes" if r["features"].get("promoter") else "no",
            "regulatory_weight": reg_w,
            "comparisons": ",".join(f"{c}:{n}" for c, n in sorted(r["comparisons"].items())),
            "n_comparisons": len(r["comparisons"]),
            "features": ",".join(f"{f}:{n}" for f, n in sorted(r["features"].items())),
        }
    return out


# ----------------------------------------------------------------- sub03 (expressão)
def _deg_dir():
    for d in GENOMAS_DEG_DIRS:
        if d.exists() and list(d.glob("*_DEGs.tsv")):
            return d
    return None


def load_sub03_expression():
    """gene -> dict(tissues_sig, best_log2fc, expr_direction, n_tissues)."""
    d = _deg_dir()
    genes = defaultdict(lambda: {"tissues": {}, "log2fcs": []})
    if not d:
        return {}, []
    tissues_seen = []
    for f in sorted(d.glob("*_DEGs.tsv")):
        tissue = f.name.replace("_DEGs.tsv", "")
        tissues_seen.append(tissue)
        with f.open(encoding="utf-8", errors="ignore") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                g = (row.get("gene_name") or "").strip()
                if not g:
                    continue
                try:
                    lfc = float(row.get("log2FoldChange"))
                    padj = float(row.get("padj"))
                except (TypeError, ValueError):
                    continue
                if padj < DEG_PADJ:
                    genes[g]["tissues"][tissue] = round(lfc, 2)
                    genes[g]["log2fcs"].append(lfc)
    out = {}
    for g, r in genes.items():
        if not r["log2fcs"]:
            continue
        best = max(r["log2fcs"], key=abs)
        out[g] = {
            "expr_tissues_sig": ",".join(sorted(r["tissues"])),
            "n_tissues_sig": len(r["tissues"]),
            "best_log2fc": round(best, 2),
            "expr_direction": "up" if best > 0 else "down",
        }
    return out, tissues_seen


# ----------------------------------------------------------------- sub04 (painel)
def load_panel():
    for f in SUB04_PANEL:
        if f.exists():
            panel = {}
            with f.open(encoding="utf-8", errors="ignore") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    g = (row.get("gene") or "").strip()
                    if g:
                        panel[g] = {"n_tissues_sig": row.get("n_tissues_sig", ""),
                                    "is_AIM015": row.get("is_AIM015", "")}
            return panel
    return {}


# ----------------------------------------------------------------- integração
def concordance(meth_dir, expr_dir):
    """Concordância canônica (metilação repressiva): hyper+down / hypo+up."""
    if meth_dir is None or expr_dir is None:
        return "NA"
    if (meth_dir == "hyper" and expr_dir == "down") or \
       (meth_dir == "hypo" and expr_dir == "up"):
        return "concordant_repressive"
    return "discordant"


def main():
    meth = load_sub02_methylation()
    expr, tissues = load_sub03_expression()
    panel = load_panel()

    def priority_score(m, e, conc_label):
        """Escore de prioridade (0-6) p/ candidatos de validação."""
        if not (m and e):
            return 0
        s = 0
        s += 2 if conc_label == "concordant_repressive" else 0
        s += 1 if m["has_promoter_dmc"] == "yes" else 0          # promotor = repressão canônica
        s += min(m["n_comparisons"], 2)                          # robustez entre comparações
        s += 1 if e["n_tissues_sig"] >= 2 else 0                 # robustez entre tecidos
        return s

    all_genes = sorted(set(meth) | set(expr))
    rows = []
    for g in all_genes:
        m = meth.get(g)
        e = expr.get(g)
        p = panel.get(g, {})
        conc_label = concordance(m["meth_direction"] if m else None,
                                 e["expr_direction"] if e else None)
        score = priority_score(m, e, conc_label)
        tier = ("—" if not (m and e) else
                "alta" if score >= 5 else "média" if score >= 3 else "baixa")
        rows.append({
            "gene": g,
            "in_methylation_SUB02": "yes" if m else "no",
            "n_DMCs": m["n_dmcs"] if m else "",
            "mean_methdiff": m["mean_methdiff"] if m else "",
            "reg_methdiff": m["reg_methdiff"] if m else "",
            "meth_direction": m["meth_direction"] if m else "",
            "has_promoter_dmc": m["has_promoter_dmc"] if m else "",
            "regulatory_weight": m["regulatory_weight"] if m else "",
            "meth_features": m["features"] if m else "",
            "meth_comparisons": m["comparisons"] if m else "",
            "n_comparisons": m["n_comparisons"] if m else "",
            "in_expression_sub03": "yes" if e else "no",
            "expr_best_log2fc": e["best_log2fc"] if e else "",
            "expr_direction": e["expr_direction"] if e else "",
            "expr_tissues_sig": e["expr_tissues_sig"] if e else "",
            "n_tissues_sig": e["n_tissues_sig"] if e else "",
            "concordance": conc_label,
            "priority_score": score,
            "priority_tier": tier,
            "panel73_n_tissues_sig": p.get("n_tissues_sig", ""),
            "panel73_is_AIM": p.get("is_AIM015", ""),
        })

    # ordena: cruzamentos primeiro, depois por prioridade desc
    rank = {"yes": 0, "no": 1}
    rows.sort(key=lambda r: (rank[r["in_methylation_SUB02"]] + rank[r["in_expression_sub03"]],
                             -r["priority_score"], r["gene"]))
    cols = list(rows[0].keys()) if rows else []
    with (OUT / "bridge_gene_table.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    both = [r for r in rows if r["in_methylation_SUB02"] == "yes"
            and r["in_expression_sub03"] == "yes"]
    conc = [r for r in both if r["concordance"] == "concordant_repressive"]
    prom_conc = [r for r in conc if r["has_promoter_dmc"] == "yes"]

    # sumário
    lines = [
        "PONTE MULTI-ÔMICA  —  epigenoma SUB02 (DMR) × genomas sub03 (DEG) + sub04",
        "=" * 72,
        f"Comparações SUB02 lidas: {sorted({c.split(':')[0] for r in both for c in (r['meth_comparisons'] or '').split(',') if c})}",
        f"Tecidos RNA-seq (sub03): {', '.join(tissues) or '—'}",
        f"Genes com DMC (SUB02, |Δβ|>15%): {len(meth)}",
        f"Genes DE (sub03, padj<{DEG_PADJ}): {len(expr)}",
        f"Genes em AMBAS as camadas: {len(both)}",
        f"  concordância repressiva (direção regulatória hiper↓/hipo↑): {len(conc)}",
        f"  dos quais com DMC em PROMOTOR (sinal mais interpretável): {len(prom_conc)}",
        "",
        "Refino: a direção da metilação usa preferencialmente promotor/upstream;",
        "priority_score (0-6) = concordância(2) + promotor(1) + nº comparações(≤2) + ≥2 tecidos(1).",
        "",
        "Candidatos prioritários (tier alta/média), ordenados:",
    ]
    for r in [x for x in both if x["priority_tier"] in ("alta", "média")]:
        lines.append(
            f"  [{r['priority_tier']:<5} {r['priority_score']}] {r['gene']:<10} "
            f"meth={r['meth_direction']:<5}(reg Δβ={r['reg_methdiff']}, "
            f"prom={r['has_promoter_dmc']}, {r['n_comparisons']} comp) | "
            f"expr={r['expr_direction']:<4}(log2FC={r['expr_best_log2fc']}, "
            f"{r['n_tissues_sig']} tec) | {r['concordance']}")
    (OUT / "bridge_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # figura de concordância
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        if both:
            fig, ax = plt.subplots(figsize=(8.4, 6.4))
            for r in both:
                # eixo x: direção regulatória (promotor/upstream) quando houver
                x = r["reg_methdiff"] if r["reg_methdiff"] != "" else r["mean_methdiff"]
                y = r["expr_best_log2fc"]
                is_conc = r["concordance"] == "concordant_repressive"
                col = "#1e8449" if is_conc else "#c0392b"
                size = 50 + 35 * r["priority_score"]          # tamanho ∝ prioridade
                edge = "#111" if r["has_promoter_dmc"] == "yes" else "white"
                ax.scatter(x, y, s=size, color=col, edgecolor=edge,
                           linewidth=1.3, zorder=3, alpha=0.9)
                lbl = r["gene"] + ("*" if r["has_promoter_dmc"] == "yes" else "")
                ax.annotate(lbl, (x, y), fontsize=8, xytext=(4, 4),
                            textcoords="offset points")
            ax.axhline(0, color="#888", lw=0.8)
            ax.axvline(0, color="#888", lw=0.8)
            ax.set_xlabel("Metilação SUB02 — Δβ regulatório (promotor/upstream) %  (hipo ← | → hiper)")
            ax.set_ylabel("Expressão sub03 — log2FC (down ← | → up)")
            ax.set_title("Ponte multi-ômica: metilação (raça×calor) × expressão (heat stress)\n"
                         "verde = concordância repressiva · tamanho ∝ prioridade · "
                         "* = DMC em promotor", fontsize=9.5)
            fig.tight_layout()
            fig.savefig(OUT / "bridge_concordance.png", dpi=200)
            plt.close(fig)
    except Exception as exc:
        print(f"[warn] figura não gerada: {exc}")

    print(f"[bridge] genes metilados SUB02 : {len(meth)}")
    print(f"[bridge] genes DE sub03        : {len(expr)}")
    print(f"[bridge] interseção            : {len(both)}  (concordantes: {len(conc)})")
    print(f"[bridge] saída                 : {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
