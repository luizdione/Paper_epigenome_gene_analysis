"""
DMC functional distribution
===========================
Classifica cada DMC (Fase 2 intersect) em 4 buckets relativos ao gene-alvo:

    upstream_distal  : dist_to_TSS < -2000
    proximal_promoter: -2000 <= dist_to_TSS <= +500
    gene_body        : +500 < dist_to_TSS < (gene_end - gene_start)
    downstream_distal: dist_to_TSS >= (gene_end - gene_start)

Gera:
    - dmc_functional_distribution.csv  (uma linha por DMC + bucket + gene_length)
    - dmc_functional_distribution.png  (barplot dos 4 buckets)
    - dmc_functional_by_gene.png       (stacked bar por gene)
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PROJ = Path(__file__).resolve().parents[2]
DMCS_CSV = PROJ / "results" / "SUB01_GSE161113_run04" / "fase2_intersect" / "dmrs_in_target_genes.csv"
GENES_JSON = PROJ / "genes_alvo" / "all_target_genes_coordinates.json"
OUT_DIR = PROJ / "results" / "SUB01_GSE161113_run04" / "fase2_intersect"


def classify(d_tss: int, gene_len: int) -> str:
    if d_tss < -2000:
        return "upstream_distal"
    if -2000 <= d_tss <= 500:
        return "proximal_promoter"
    if 500 < d_tss < gene_len:
        return "gene_body"
    return "downstream_distal"


def main() -> int:
    if not DMCS_CSV.exists():
        sys.exit(f"DMC CSV not found: {DMCS_CSV}")
    if not GENES_JSON.exists():
        sys.exit(f"Targets JSON not found: {GENES_JSON}")

    genes = json.loads(GENES_JSON.read_text(encoding="utf-8"))["genes"]
    gene_len = {g: (info["end"] - info["start"]) for g, info in genes.items()}

    rows: list[dict] = []
    with DMCS_CSV.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            gene = r["gene"]
            d = int(r["dist_to_TSS"])
            glen = gene_len.get(gene)
            if glen is None:
                continue
            bucket = classify(d, glen)
            rows.append({
                "chr": r["chr"], "pos": r["start"], "gene": gene,
                "module": r["module"], "meth_diff": r["meth_diff"],
                "qvalue": r["qvalue"], "dist_to_TSS": d,
                "gene_length": glen, "bucket": bucket,
            })

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    out_csv = OUT_DIR / "dmc_functional_distribution.csv"
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"[OK] {len(rows)} DMCs classified -> {out_csv}")

    order = ["upstream_distal", "proximal_promoter", "gene_body", "downstream_distal"]
    label_map = {
        "upstream_distal": "Upstream distal\n(< -2 kb TSS)",
        "proximal_promoter": "Proximal promoter\n(-2 kb to +0.5 kb)",
        "gene_body": "Gene body\n(+0.5 kb to TES)",
        "downstream_distal": "Downstream distal\n(> TES)",
    }
    color_map = {
        "upstream_distal": "#1f77b4",
        "proximal_promoter": "#d62728",
        "gene_body": "#2ca02c",
        "downstream_distal": "#ff7f0e",
    }

    counts = Counter(r["bucket"] for r in rows)
    vals = [counts.get(b, 0) for b in order]
    total = sum(vals)
    pcts = [v / total * 100 if total else 0 for v in vals]

    fig, ax = plt.subplots(figsize=(8, 5.5))
    bars = ax.bar(range(len(order)), vals,
                  color=[color_map[b] for b in order],
                  edgecolor="black", linewidth=0.5)
    for i, (v, p) in enumerate(zip(vals, pcts)):
        ax.text(i, v + max(vals) * 0.02, f"n={v}\n({p:.1f}%)",
                ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([label_map[b] for b in order], fontsize=10)
    ax.set_ylabel("Number of DMCs", fontsize=11)
    ax.set_title(f"Functional distribution of {total} DMCs\n"
                 f"SUB01 — Bos taurus (Angus) vs Bos indicus (Nelore) — RRBS GSE161113",
                 fontsize=11)
    ax.set_ylim(0, max(vals) * 1.18 if vals else 1)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    out_png = OUT_DIR / "dmc_functional_distribution.png"
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] -> {out_png}")

    by_gene: dict[str, Counter] = defaultdict(Counter)
    for r in rows:
        by_gene[r["gene"]][r["bucket"]] += 1
    genes_sorted = sorted(by_gene.keys(),
                          key=lambda g: sum(by_gene[g].values()),
                          reverse=True)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    bottoms = np.zeros(len(genes_sorted))
    x = np.arange(len(genes_sorted))
    for b in order:
        heights = np.array([by_gene[g].get(b, 0) for g in genes_sorted])
        ax.bar(x, heights, bottom=bottoms,
               label=label_map[b].replace("\n", " "),
               color=color_map[b], edgecolor="black", linewidth=0.4)
        bottoms = bottoms + heights
    ax.set_xticks(x)
    ax.set_xticklabels(genes_sorted, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Number of DMCs", fontsize=11)
    ax.set_title("Per-gene functional distribution of DMCs (SUB01)", fontsize=11)
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    out_png2 = OUT_DIR / "dmc_functional_by_gene.png"
    plt.savefig(out_png2, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] -> {out_png2}")

    print("\n=== SUMMARY ===")
    for b in order:
        print(f"  {b:<20s}: n={counts.get(b, 0)}  ({counts.get(b, 0)/total*100:.1f}%)")
    print(f"  Total: {total} DMCs across {len(by_gene)} genes\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
