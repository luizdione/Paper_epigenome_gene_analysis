#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
export_online_resources.py — exporta as tabelas suplementares (Online Resources 1..9) a partir
das saídas versionadas em results/SUB02_runs_108/ e results/ancestry/. Idempotente.

Renomeia colunas internas do pipeline para nomes publicáveis e alinha a flag `high_tier`
(OR06, OR08, OR09) aos dez candidatos high-tier finais (HIGH). OR10 (painel curado) e OR11
(enriquecimento de motivos) são tabelas curadas mantidas diretamente em tables/.
Uso: python scripts/python/export_online_resources.py   (na raiz do repositório)
"""
import csv
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
R = os.path.join(ROOT, "results", "SUB02_runs_108")
AIMS = os.path.join(ROOT, "results", "ancestry", "degs_x_aims_genes.tsv")
OUT = os.path.join(ROOT, "tables")

HIGH = {"GHSR", "HSPA2", "HSPB1", "PRLR", "ITGB3", "ERN1", "GRHL3", "PADI4", "LGALS3", "ISG15"}
RENAME = {"in_methylation_SUB02": "in_methylation", "in_expression_sub03": "in_expression",
          "sub02_meth_dir": "rrbs_meth_dir", "max_Fst_GirHol": "max_Fst_GyrHol"}
HIGH_TIER_TABLES = {6, 8, 9}
NAMES = {1: "OR01_bridge_overlap47.tsv", 2: "OR02_bridge_gene_table.tsv", 3: "OR03_concordant_by_contrast.tsv",
         4: "OR04_module_enrichment.tsv", 5: "OR05_fetalliver_C6.tsv", 6: "OR06_cissnp_1000bull.tsv",
         7: "OR07_Fst_AIM_context.tsv", 8: "OR08_matched_atlas_Zhou.tsv", 9: "OR09_motifbreakR.tsv"}


def rd(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def T(v):
    return str(v).strip().lower() in ("1", "true", "yes", "y")


def write(rows, cols, n):
    cols = [RENAME.get(c, c) for c in cols]
    out = []
    for r in rows:
        r = {RENAME.get(k, k): v for k, v in r.items()}
        if n in HIGH_TIER_TABLES and "high_tier" in r:
            r["high_tier"] = "yes" if r.get("gene") in HIGH else "no"
        out.append(r)
    if n in HIGH_TIER_TABLES:
        out.sort(key=lambda r: r.get("high_tier") != "yes")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, NAMES[n]), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(out)
    print("%s (%d rows)" % (NAMES[n], len(out)))


def main():
    br = rd(os.path.join(R, "_bridge_multiomic", "bridge_gene_table.tsv")); bcols = list(br[0].keys())
    write(br, bcols, 2)
    overlap = [r for r in br if T(r.get("in_methylation_SUB02")) and T(r.get("in_expression_sub03"))]
    overlap.sort(key=lambda r: -float(r.get("priority_score") or 0))
    write(overlap, bcols, 1)
    bc = rd(os.path.join(R, "_bridge_multiomic", "bridge_by_comparison.tsv")); write(bc, list(bc[0].keys()), 3)
    me = rd(os.path.join(R, "_bridge_multiomic", "bridge_module_enrichment.tsv")); write(me, list(me[0].keys()), 4)
    c6 = rd(os.path.join(R, "C6_GIGAS_BIBI_vs_BTBT_fetalLiver", "c6_validation_table.tsv")); write(c6, list(c6[0].keys()), 5)
    cs = rd(os.path.join(R, "_cissnp", "cissnp_table.tsv")); write(cs, list(cs[0].keys()), 6)
    at = rd(os.path.join(R, "_matched_atlas", "matched_correlation.tsv")); write(at, list(at[0].keys()), 8)
    mo = rd(os.path.join(R, "_motif_dmc", "motif_dmc_table.tsv")); write(mo, list(mo[0].keys()), 9)
    aim = {r["gene"]: r for r in rd(AIMS)}
    rows = []
    for r in overlap:
        a = aim.get(r["gene"], {}); in30 = T(a.get("in_AIM_fst030"))
        rows.append({"gene": r["gene"], "priority_tier": r.get("priority_tier", ""),
                     "concordance": r.get("concordance", ""), "in_AIM_fst030": a.get("in_AIM_fst030", "NA"),
                     "in_AIM_fst045_elite": a.get("in_AIM_fst045_elite", "NA"),
                     "max_Fst_GirHol": a.get("max_Fst_GirHol", "NA"),
                     "genomic_context": "triple_convergence" if in30 else "epi_dominant"})
    rows.sort(key=lambda r: (r["genomic_context"] != "triple_convergence", r["gene"]))
    write(rows, ["gene", "priority_tier", "concordance", "in_AIM_fst030", "in_AIM_fst045_elite",
                 "max_Fst_GirHol", "genomic_context"], 7)


if __name__ == "__main__":
    main()
