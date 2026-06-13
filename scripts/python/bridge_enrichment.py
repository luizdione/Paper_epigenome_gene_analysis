#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bridge_enrichment.py — enriquece a ponte multi-ômica:
  (a) sobre-representação por MÓDULO FUNCIONAL dos genes da interseção
      metilação×expressão (teste hipergeométrico vs o painel-alvo);
  (b) breakdown POR COMPARAÇÃO (quantos genes da ponte cada contraste sustenta);
  (c) figura de módulos.

Lê apenas artefatos versionados:
  results/SUB02_runs/_bridge_multiomic/bridge_gene_table.tsv
  subprojeto/SUB02/SUB02_target_genes_coordinates.json  (gene -> module)
Saída: results/SUB02_runs/_bridge_multiomic/
  bridge_module_enrichment.tsv, bridge_by_comparison.tsv,
  bridge_modules.png, bridge_enrichment_summary.txt
"""
from __future__ import annotations
import csv, json
from collections import defaultdict
from pathlib import Path
from scipy.stats import hypergeom

BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "results" / "SUB02_runs" / "_bridge_multiomic"
TABLE = OUT / "bridge_gene_table.tsv"
JSON = BASE / "subprojeto" / "SUB02" / "SUB02_target_genes_coordinates.json"

MODULE_LABEL = {
    "thermotolerance": "Termotolerância (HSP)",
    "thermo_reproduction": "Termo-reprodução",
    "milk_production": "Produção de leite",
}


def gene_modules():
    d = json.loads(JSON.read_text(encoding="utf-8"))
    return {g: (info.get("module") or "outros")
            for g, info in d.get("genes", {}).items()}


def main():
    mods = gene_modules()
    panel = sorted(mods)                      # universo = painel-alvo
    rows = list(csv.DictReader(TABLE.open(encoding="utf-8"), delimiter="\t"))

    both = [r for r in rows if r["in_methylation_SUB02"] == "yes"
            and r["in_expression_sub03"] == "yes"]
    conc = [r for r in both if r["concordance"] == "concordant_repressive"]
    conc_genes = [r["gene"] for r in conc]
    meth_genes = [r["gene"] for r in rows if r["in_methylation_SUB02"] == "yes"]

    # ---- (a) enriquecimento por módulo (hipergeométrico) ----
    N = len(panel)
    n = len([g for g in conc_genes if g in mods])     # amostra (concordantes no painel)
    by_mod_panel = defaultdict(list)
    for g in panel:
        by_mod_panel[mods[g]].append(g)
    by_mod_conc = defaultdict(list)
    for g in conc_genes:
        if g in mods:
            by_mod_conc[mods[g]].append(g)

    enr = []
    for mod, gs in sorted(by_mod_panel.items(), key=lambda x: -len(by_mod_conc.get(x[0], []))):
        K = len(gs)                  # genes do módulo no painel
        k = len(by_mod_conc.get(mod, []))
        # P(X >= k) sob hipergeométrica
        p = float(hypergeom.sf(k - 1, N, K, n)) if k > 0 else 1.0
        frac = round(k / K, 3) if K else 0.0
        enr.append({"module": mod, "label": MODULE_LABEL.get(mod, mod),
                    "genes_in_panel": K, "concordant": k,
                    "frac_module_concordant": frac,
                    "concordant_genes": ",".join(sorted(by_mod_conc.get(mod, []))),
                    "hypergeom_p": round(p, 4)})
    with (OUT / "bridge_module_enrichment.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(enr[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(enr)

    # ---- (b) breakdown por comparação ----
    comp_genes = defaultdict(set)
    comp_conc = defaultdict(set)
    for r in both:
        for tok in (r["meth_comparisons"] or "").split(","):
            if ":" in tok:
                comp = tok.split(":")[0]
                comp_genes[comp].add(r["gene"])
                if r["concordance"] == "concordant_repressive":
                    comp_conc[comp].add(r["gene"])
    comp_rows = []
    for comp in sorted(comp_genes, key=lambda c: -len(comp_genes[c])):
        comp_rows.append({"comparison": comp,
                          "bridge_genes": len(comp_genes[comp]),
                          "concordant_genes": len(comp_conc[comp]),
                          "genes": ",".join(sorted(comp_genes[comp]))})
    with (OUT / "bridge_by_comparison.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(comp_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(comp_rows)

    # ---- (c) figura de módulos ----
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        labels = [e["label"] for e in enr]
        kk = [e["concordant"] for e in enr]
        KK = [e["genes_in_panel"] for e in enr]
        import numpy as np
        y = np.arange(len(labels))
        fig, ax = plt.subplots(figsize=(7.6, 0.7 * len(labels) + 1.6))
        ax.barh(y, KK, color="#d6dbdf", label="genes no painel")
        ax.barh(y, kk, color="#1e8449", label="concordantes na ponte")
        for yi, k, K, e in zip(y, kk, KK, enr):
            ax.text(K + 0.1, yi, f"{k}/{K}  (p={e['hypergeom_p']})", va="center", fontsize=8)
        ax.set_yticks(y); ax.set_yticklabels(labels); ax.invert_yaxis()
        ax.set_xlabel("nº de genes-alvo")
        ax.set_title("Ponte: enriquecimento por módulo funcional\n"
                     "(genes concordantes metilação×expressão)", fontsize=10)
        ax.legend(frameon=False, fontsize=8.5, loc="lower right")
        fig.tight_layout(); fig.savefig(OUT / "bridge_modules.png", dpi=200); plt.close(fig)
    except Exception as exc:
        print("[warn] figura:", exc)

    # ---- sumário ----
    lines = ["ENRIQUECIMENTO DA PONTE — módulos funcionais + por comparação",
             "=" * 64,
             f"Genes concordantes (universo = painel de {N} genes): {n}", ""]
    for e in enr:
        lines.append(f"  {e['label']:<24} {e['concordant']}/{e['genes_in_panel']} "
                     f"(p={e['hypergeom_p']})  {e['concordant_genes']}")
    lines += ["", "Por comparação (genes da ponte / concordantes):"]
    for r in comp_rows:
        lines.append(f"  {r['comparison']:<34} {r['bridge_genes']:>2} / {r['concordant_genes']:>2}")
    (OUT / "bridge_enrichment_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
