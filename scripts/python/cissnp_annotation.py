#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cissnp_annotation.py — camada CIS-GENÉTICA (atalho leve, sem download).

Para cada DMC do painel nos contrastes de SUBESPÉCIE, verifica no catálogo
1000 Bull Genomes (TAUIND, taurino+indicino, 1842 amostras; tabix por cromossomo)
se há um SNP segregante:
  - SOBRE o CpG (|POS - dmc| <= 1)  -> "CpG-SNP": cria/destrói o sítio = mecanismo
    passivo mais direto (a diferença de metilação pode ser consequência da sequência);
  - PRÓXIMO (|POS - dmc| <= 10)      -> "near-SNP";
  - senão                            -> "no-SNP" (compatível com metilação genuína).

Isso quantifica um LIMITE SUPERIOR da fração de DMCs potencialmente dirigida por
genótipo (caveat passivo; Gutiérrez-Arcelus et al. 2013), sem alegar causalidade.

Roda DENTRO do container (pysam + /data/1000bg_vcf montado).
Saídas: results/SUB02_runs/_cissnp/
"""
from __future__ import annotations
import csv, glob, os
from collections import defaultdict
from pathlib import Path
import pysam

PROJ = Path(".")
RUNS = PROJ / "results" / "SUB02_runs"
VCFDIR = Path("/data/1000bg_vcf")
OUT = RUNS / "_cissnp"; OUT.mkdir(parents=True, exist_ok=True)

BREED_CONTRASTS = [
    "C1_zenodo_NEL_vs_HF", "C2_GSE_NelFev_vs_AngFev_HS", "C3_GSE_NelJun_vs_AngJun",
    "C7_bonus_zenodo_NEL_vs_NDama", "C8_zenodo_NEL_vs_HF_CD4", "C9_zenodo_NEL_vs_HF_CD8",
    "C10_zenodo_NEL_vs_HF_Bcell", "C11_zenodo_NEL_vs_HF_gdT", "C12_zenodo_NEL_vs_HF_gran",
    "C13_zenodo_NEL_vs_HF_mono", "C14_zenodo_NEL_vs_HF_NK",
]
HIGH_TIER = {"GHSR", "HSPA2", "HSPB1", "PRLR", "ITGB3", "ERN1", "GRHL3", "PADI4", "LGALS3", "ISG15"}
CPG_WIN, NEAR_WIN = 1, 10


def load_dmcs():
    """(chrom,pos) -> {gene, feature, contrasts:set, methdiffs:[]}"""
    d = defaultdict(lambda: {"gene": "", "feature": "", "contrasts": set(), "md": []})
    for c in BREED_CONTRASTS:
        f = RUNS / c / "fase2_annotation_diff15" / "annotation_table.csv"
        if not f.exists():
            continue
        with f.open(encoding="utf-8", errors="ignore") as fh:
            for row in csv.DictReader(fh):
                chrom = (row.get("chr") or "").strip().strip('"')
                try:
                    pos = int(float(row.get("start")))
                    md = float(row.get("meth.diff"))
                except (TypeError, ValueError):
                    continue
                g = (row.get("gene") or "").strip().strip('"')
                feat = (row.get("feature") or "").strip().strip('"')
                rec = d[(chrom, pos)]
                rec["gene"] = rec["gene"] or g
                rec["feature"] = rec["feature"] or feat
                rec["contrasts"].add(c); rec["md"].append(md)
    return d


def vcf_for(chrom):
    p = VCFDIR / f"Chr{chrom}-Run8-TAUIND-public.vcf.gz"
    return str(p) if p.exists() else None


def main():
    dmcs = load_dmcs()
    print(f"[cissnp] DMCs únicos (posição) nos contrastes de subespécie: {len(dmcs)}")
    handles = {}
    rows = []
    for (chrom, pos), rec in dmcs.items():
        vp = vcf_for(chrom)
        cls = "no-SNP"; nearest = ""
        if vp:
            if chrom not in handles:
                try:
                    handles[chrom] = pysam.TabixFile(vp)
                except Exception:
                    handles[chrom] = None
            tb = handles[chrom]
            if tb is not None:
                best = None
                try:
                    for line in tb.fetch(chrom, max(pos - NEAR_WIN - 1, 0), pos + NEAR_WIN + 1):
                        vpos = int(line.split("\t", 2)[1])
                        dist = abs(vpos - pos)
                        if best is None or dist < best:
                            best = dist
                except (ValueError, OSError):
                    pass
                if best is not None:
                    nearest = best
                    cls = "CpG-SNP" if best <= CPG_WIN else ("near-SNP" if best <= NEAR_WIN else "no-SNP")
        rows.append({"chrom": chrom, "pos": pos, "gene": rec["gene"],
                     "feature": rec["feature"], "n_contrasts": len(rec["contrasts"]),
                     "mean_methdiff": round(sum(rec["md"]) / len(rec["md"]), 1),
                     "class": cls, "nearest_snp_bp": nearest,
                     "high_tier": "yes" if rec["gene"] in HIGH_TIER else "no"})

    rows.sort(key=lambda r: (r["high_tier"] != "yes", r["gene"], r["chrom"], r["pos"]))
    with (OUT / "cissnp_table.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    n = len(rows)
    cpg = sum(1 for r in rows if r["class"] == "CpG-SNP")
    near = sum(1 for r in rows if r["class"] == "near-SNP")
    none = sum(1 for r in rows if r["class"] == "no-SNP")
    # por gene high-tier
    ht = defaultdict(lambda: [0, 0])  # gene -> [cpg_or_near, total]
    for r in rows:
        if r["high_tier"] == "yes":
            ht[r["gene"]][1] += 1
            if r["class"] in ("CpG-SNP", "near-SNP"):
                ht[r["gene"]][0] += 1
    lines = [
        "CAMADA CIS-GENÉTICA — DMCs de subespécie × SNPs 1000 Bull (TAUIND, 1842 bovinos)",
        "=" * 76,
        f"DMCs (posições únicas) avaliados: {n}",
        f"  CpG-SNP (SNP sobre o CpG, |Δ|≤1 pb): {cpg} ({100*cpg/n:.1f}%) — passivo mais direto",
        f"  near-SNP (SNP a ≤10 pb): {near} ({100*near/n:.1f}%)",
        f"  no-SNP (sem SNP em ±10 pb): {none} ({100*none/n:.1f}%) — compatível c/ metilação genuína",
        "",
        f"Limite superior 'potencialmente dirigido por genótipo' (CpG-SNP+near): "
        f"{cpg+near}/{n} ({100*(cpg+near)/n:.1f}%)",
        f"Fração mais conservadora (só CpG-SNP): {cpg}/{n} ({100*cpg/n:.1f}%)",
        "",
        "High-tier (gene: DMCs com SNP em ±10pb / total):",
    ]
    for g in sorted(ht):
        a, t = ht[g]
        lines.append(f"  {g:7s} {a}/{t}")
    (OUT / "cissnp_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n[cissnp] saídas em {OUT}")


if __name__ == "__main__":
    main()
