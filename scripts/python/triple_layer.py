#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
triple_layer.py — SINERGIA com a plataforma genomas.

Cruza os candidatos da ponte (metilação SUB02 × expressão sub03) com a camada
GENÔMICA do genomas sub04 (Fst Gir×Holstein / AIMs — marcadores de ancestralidade).
Classifica cada candidato em:
  - "triple"  -> também sob diferenciação genômica (AIM Fst≥0.30): convergência
                 genômica + epigenômica + transcriptômica;
  - "epi-dominant" -> metilado e expresso, mas SEM diferenciação genômica forte
                 (não-AIM): a divergência é regulatória, não um rastro passivo da
                 seleção de sequência — argumento direto contra o caveat passivo.

Sem números inventados: lê bridge_gene_table.tsv e degs_x_aims_genes.tsv reais.
Saída: results/SUB02_runs/_triple_layer/
"""
from __future__ import annotations
import csv
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
BRIDGE = BASE / "results" / "SUB02_runs" / "_bridge_multiomic" / "bridge_gene_table.tsv"
AIMS = Path("E:/data/results/sub04_reverse_aims/degs_x_aims_genes.tsv")
OUT = BASE / "results" / "SUB02_runs" / "_triple_layer"; OUT.mkdir(parents=True, exist_ok=True)
HIGH_TIER = {"PRLR","GHSR","HSPA2","HSPB1","ISG15","MX1","ITGB3","KRT18"}


def load_aims():
    d = {}
    with AIMS.open(encoding="utf-8", errors="ignore") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            d[row["gene"]] = row
    return d


def main():
    aims = load_aims()
    rows = []
    with BRIDGE.open(encoding="utf-8", errors="ignore") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["in_methylation_SUB02"] != "yes" or r["in_expression_sub03"] != "yes":
                continue
            g = r["gene"]; a = aims.get(g, {})
            in_aim = (a.get("in_AIM_fst030", "").strip() in ("True", "true", "1"))
            elite = (a.get("in_AIM_fst045_elite", "").strip() in ("True", "true", "1"))
            fst = a.get("max_Fst_GirHol", "")
            layer = "triple" if in_aim else "epi-dominant"
            rows.append({
                "gene": g, "priority_tier": r.get("priority_tier", ""),
                "meth_direction": r.get("meth_direction", ""),
                "concordance": r.get("concordance", ""),
                "expr_best_log2fc": r.get("expr_best_log2fc", ""),
                "in_AIM_fst030": "yes" if in_aim else "no",
                "AIM_elite_fst045": "yes" if elite else "no",
                "max_Fst_GirHol": fst,
                "layer_class": layer,
                "high_tier": "yes" if g in HIGH_TIER else "no",
            })
    rows.sort(key=lambda x: (x["high_tier"] != "yes", x["layer_class"], x["gene"]))
    with (OUT / "triple_layer_table.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    n = len(rows)
    triple = [r for r in rows if r["layer_class"] == "triple"]
    epi = [r for r in rows if r["layer_class"] == "epi-dominant"]
    ht = [r for r in rows if r["high_tier"] == "yes"]
    ht_epi = [r for r in ht if r["layer_class"] == "epi-dominant"]
    lines = [
        "SINERGIA genômica (genomas sub04 Fst/AIM) × ponte epigenoma×transcriptoma",
        "=" * 74,
        f"Candidatos em ambas as camadas ômicas (ponte): {n}",
        f"  TRIPLE (também AIM, Fst≥0.30 — genômico+epi+transcriptômico): {len(triple)}",
        f"    -> {', '.join(r['gene'] for r in triple)}",
        f"  EPI-DOMINANT (metilado+expresso, NÃO-AIM): {len(epi)}",
        f"    -> divergência regulatória sem diferenciação genômica forte (anti-caveat passivo)",
        "",
        f"High-tier: {len(ht)} candidatos; epi-dominant: {len(ht_epi)} "
        f"({', '.join(r['gene'] for r in ht_epi)})",
        "",
        "Detalhe high-tier (gene | tier | layer | Fst Gir×Hol):",
    ]
    for r in ht:
        lines.append(f"  {r['gene']:7s} | {r['priority_tier']:5s} | {r['layer_class']:12s} "
                     f"| Fst={r['max_Fst_GirHol'] or '—'}")
    (OUT / "triple_layer_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n[triple] saídas em {OUT}")


if __name__ == "__main__":
    main()
