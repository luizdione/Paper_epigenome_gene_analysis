#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
motif_dmc_cross.py — camada MECANÍSTICA (sinergia genomas sub01 motifbreakR).

Cruza os SNPs que ROMPEM os motifs de TF que enriquecemos nos flancos das DMCs
(SP1, NF-κB p65/RELA, NF-κB p50/NFKB1, KLF4, NFAT) — preditos pelo motifbreakR do
genomas sub01 — com as posições dos nossos DMCs de subespécie. Um DMC que coincide
(±25 pb) com um SNP que quebra um desses motifs é um sítio regulatório alterado
SIMULTANEAMENTE no nível de sequência (genético) e de metilação (epigenético):
a ligação mecanística direta entre as duas camadas.

Lê resultados reais; nada inventado. Saída: results/SUB02_runs/_motif_dmc/
"""
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import csv
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
MBR = Path("E:/data/results/sub01/regulatory_motifs/motifbreakr_results.tsv")
DMCS = BASE / "results" / "SUB02_runs" / "_cissnp" / "cissnp_table.tsv"
OUT = BASE / "results" / "SUB02_runs" / "_motif_dmc"; OUT.mkdir(parents=True, exist_ok=True)

# TFs enriquecidos nos nossos DMCs (nomes HOCOMOCO/motifbreakR)
TARGET_TF = {"SP1": "SP1", "KLF4": "KLF4", "RELA": "NF-κB p65", "NFKB1": "NF-κB p50",
             "NFKB2": "NF-κB", "REL": "NF-κB", "NFAC1": "NFAT", "NFATC1": "NFAT",
             "NFAC2": "NFAT", "NFATC2": "NFAT"}
TIGHT, LOOSE = 10, 25
HIGH_TIER = {"PRLR","GHSR","HSPA2","HSPB1","ISG15","MX1","ITGB3","KRT18"}


def norm_chr(c):
    c = str(c).strip()
    return c[3:] if c.lower().startswith("chr") else c


def load_motif_snps():
    """chrom -> sorted [(pos, TF_label)] para SNPs 'strong' nos TFs-alvo."""
    by_chrom = defaultdict(list)
    with MBR.open(encoding="utf-8", errors="ignore") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tf = (row.get("geneSymbol") or "").strip().upper()
            if tf not in TARGET_TF:
                continue
            if (row.get("effect") or "").strip().lower() != "strong":
                continue
            try:
                pos = int(row["pos"])
            except (TypeError, ValueError, KeyError):
                continue
            by_chrom[norm_chr(row["chr"])].append((pos, TARGET_TF[tf]))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom


def nearest(by_chrom, chrom, pos):
    arr = by_chrom.get(chrom)
    if not arr:
        return None, None
    i = bisect_left(arr, (pos,))
    best = None; bestlab = None
    for j in (i - 1, i, i + 1):
        if 0 <= j < len(arr):
            d = abs(arr[j][0] - pos)
            if best is None or d < best:
                best = d; bestlab = arr[j][1]
    return best, bestlab


def main():
    msnp = load_motif_snps()
    n_msnp = sum(len(v) for v in msnp.values())
    print(f"[motif] SNPs 'strong' que rompem SP1/NF-κB/KLF4/NFAT: {n_msnp}")

    rows = []
    with DMCS.open(encoding="utf-8", errors="ignore") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            chrom = norm_chr(r["chrom"]); pos = int(r["pos"])
            d, lab = nearest(msnp, chrom, pos)
            hit = "none"
            if d is not None:
                hit = "tight(<=10)" if d <= TIGHT else ("loose(<=25)" if d <= LOOSE else "none")
            if hit != "none":
                rows.append({"chrom": chrom, "pos": pos, "gene": r["gene"],
                             "feature": r["feature"], "motif_TF": lab,
                             "dist_bp": d, "proximity": hit,
                             "high_tier": r.get("high_tier", "no")})

    rows.sort(key=lambda x: (x["high_tier"] != "yes", x["dist_bp"]))
    with (OUT / "motif_dmc_table.tsv").open("w", newline="", encoding="utf-8") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(rows)

    n_dmc_total = sum(1 for _ in open(DMCS, encoding="utf-8")) - 1
    tight = [r for r in rows if r["proximity"].startswith("tight")]
    by_tf = defaultdict(int)
    for r in rows:
        by_tf[r["motif_TF"]] += 1
    ht_genes = sorted({r["gene"] for r in rows if r["high_tier"] == "yes"})
    lines = [
        "CAMADA MECANÍSTICA — DMCs de subespécie × SNPs que rompem motifs (motifbreakR, genomas sub01)",
        "=" * 80,
        f"DMCs de subespécie avaliados: {n_dmc_total}",
        f"DMCs coincidentes (±{LOOSE} pb) com SNP que rompe SP1/NF-κB/KLF4/NFAT: {len(rows)}",
        f"  dos quais a ≤{TIGHT} pb (mesmo motif): {len(tight)}",
        f"Por TF (DMCs com SNP-disruptor próximo): " + ", ".join(f"{k}={v}" for k, v in sorted(by_tf.items())),
        "",
        f"Genes high-tier com DMC em sítio de motif rompido por SNP: {', '.join(ht_genes) or '—'}",
        "",
        "Interpretação: nesses sítios, a divergência de metilação coincide com uma variante",
        "que altera a ligação do MESMO fator (SP1/NF-κB/KLF4) enriquecido nos flancos das DMCs —",
        "elo mecanístico direto entre as camadas genética e epigenética da divergência subespécie.",
    ]
    (OUT / "motif_dmc_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\n[motif] saídas em {OUT}")


if __name__ == "__main__":
    main()
