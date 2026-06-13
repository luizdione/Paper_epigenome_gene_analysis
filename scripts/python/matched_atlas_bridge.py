#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
matched_atlas_bridge.py — VALIDAÇÃO EM AMOSTRAS PAREADAS da lógica da ponte.

Usa o atlas Zhou/Fang (mesmas vacas Holstein; WGBS GSE147087 + RNA-seq GSE137943),
restrito à vaca 3842 (tem o maior leque de tecidos com AMBOS os ensaios).
Para cada gene do painel, correlaciona — ENTRE TECIDOS DO MESMO ANIMAL — a
metilação do promotor (WGBS) com a expressão (FPKM). A relação esperada sob
metilação repressiva é NEGATIVA. Reporta a taxa de base (quantos genes do painel
mostram a correlação negativa esperada), análogo a Sajjanar 2024 (118/275).

Sem números inventados: tudo lido dos CpG_report e do FPKM reais.
Saídas em results/SUB02_runs/_matched_atlas/.
"""
from __future__ import annotations
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import gzip, json, csv, glob, re
from pathlib import Path
from bisect import bisect_right
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import math

BASE = Path(__file__).resolve().parents[2]
ATLAS = Path("E:/data/epigenomas/Zhou_atlas")
EXTR = ATLAS / "extracted"
FPKM = ATLAS / "GSE137943_fpkm.txt.gz"
PANEL = BASE / "subprojeto" / "SUB02" / "SUB02_target_genes_coordinates.json"
OUT = BASE / "results" / "SUB02_runs" / "_matched_atlas"
OUT.mkdir(parents=True, exist_ok=True)

ANIMAL = "3842"
PROMOTER_UP, PROMOTER_DOWN = 2000, 500     # janela do promotor relativa ao TSS
MIN_COV, MIN_CPG, MIN_TISSUES = 10, 3, 6    # critérios mínimos
HIGH_TIER = {"PRLR","GHSR","HSPA2","HSPB1","ISG15","MX1","ITGB3","KRT18"}

# tecido (token do nome do arquivo WGBS) -> código da coluna FPKM
T2F = {"Adipose":"ADIP","Heart":"HEAR","Ileum":"ILEU","Kidney":"KIDN",
       "Liver":"LIVE","Lung":"LUNG","Muscle":"LDOR","Ovary":"OVAR",
       "Rumen":"RUME","Uterus":"UTER"}


def load_panel():
    d = json.load(open(PANEL, encoding="utf-8"))["genes"]
    genes = {}
    by_chrom = defaultdict(list)   # chrom -> list of (pstart, pend, gene)
    for g, r in d.items():
        chrom = str(r["chromosome"]); strand = r.get("strand", 1)
        tss = r["start"] if strand == 1 else r["end"]
        if strand == 1:
            ps, pe = tss - PROMOTER_UP, tss + PROMOTER_DOWN
        else:
            ps, pe = tss - PROMOTER_DOWN, tss + PROMOTER_UP
        genes[g] = {"ens": r.get("ensembl_id", ""), "chrom": chrom,
                    "pstart": ps, "pend": pe, "module": r.get("module", "")}
        by_chrom[chrom].append((ps, pe, g))
    for c in by_chrom:
        by_chrom[c].sort()
    return genes, by_chrom


def meth_promoter(path, by_chrom):
    """Soma meth/unmeth por gene na janela do promotor a partir de 1 CpG_report."""
    starts = {c: [w[0] for w in v] for c, v in by_chrom.items()}
    acc = defaultdict(lambda: [0, 0, 0])  # gene -> [meth, unmeth, n_cpg_cov]
    panel_chroms = set(by_chrom)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            p = line.split("\t")
            if p[0] not in panel_chroms:
                continue
            m = int(p[3]); u = int(p[4])
            if m + u == 0:
                continue
            pos = int(p[1]); wins = by_chrom[p[0]]
            i = bisect_right(starts[p[0]], pos) - 1
            # checa janelas que podem conter pos (podem haver sobreposições leves)
            j = i
            while j >= 0 and wins[j][0] <= pos:
                if wins[j][1] >= pos:
                    a = acc[wins[j][2]]; a[0] += m; a[1] += u; a[2] += 1
                j -= 1
    out = {}
    for g, (m, u, n) in acc.items():
        if m + u >= MIN_COV and n >= MIN_CPG:
            out[g] = 100.0 * m / (m + u)
    return out


def load_expr(genes):
    """ens -> {tissue_code: fpkm} para a vaca ANIMAL.
    A matriz é estilo R: cabeçalho com N nomes de amostra (sem rótulo da coluna
    de gene), e cada linha de dado tem gene_id + N valores (offset de 1)."""
    with gzip.open(FPKM, "rt") as fh:
        # matriz estilo R: separada por ESPAÇO, nomes/rownames entre aspas
        header = [tok.strip('"') for tok in next(fh).split()]
        code_to_hidx = {}
        for j, h in enumerate(header):
            if h.startswith("E1_") and h.endswith(ANIMAL):
                code_to_hidx[h[3:-len(ANIMAL)]] = j
        want_ens = {g["ens"] for g in genes.values() if g["ens"]}
        expr = {}
        offset = None
        for line in fh:
            p = line.split()
            if offset is None:
                offset = 1 if len(p) == len(header) + 1 else 0
            ens = p[0].strip().strip('"')
            if ens not in want_ens:
                continue
            d = {}
            for code, j in code_to_hidx.items():
                idx = j + offset
                if idx < len(p):
                    v = p[idx].strip().strip('"')
                    if v not in ("", "NA"):
                        try:
                            d[code] = float(v)
                        except ValueError:
                            pass
            expr[ens] = d
    return expr, sorted(code_to_hidx)


def main():
    genes, by_chrom = load_panel()
    # arquivos WGBS da vaca 3842
    files = {}
    for f in glob.glob(str(EXTR / f"*{ANIMAL}*CpG_report.txt.gz")):
        name = Path(f).name
        for tok in T2F:
            if re.search(rf"_{tok}.*{ANIMAL}", name) or re.search(rf"_{tok}{ANIMAL}", name) \
               or (tok == "Muscle" and "Ldorsi" + ANIMAL in name):
                files[tok] = f; break
    print(f"[atlas] tecidos WGBS (vaca {ANIMAL}): {sorted(files)}")

    # metilação por tecido
    meth = {}  # tissue_code -> {gene: %}
    for tok, f in sorted(files.items()):
        code = T2F[tok]
        meth[code] = meth_promoter(f, by_chrom)
        print(f"  {tok:8s} ({code}): {len(meth[code])} genes c/ metilação de promotor")

    expr, expr_codes = load_expr(genes)
    print(f"[atlas] tecidos com FPKM (vaca {ANIMAL}): {expr_codes}")

    # tecidos com AMBOS
    common = [c for c in meth if c in expr_codes]
    print(f"[atlas] tecidos pareados (WGBS ∩ RNA-seq): {common} (n={len(common)})")

    rows = []
    for g, info in genes.items():
        ens = info["ens"]
        if not ens or ens not in expr:
            continue
        xs, ys, ts = [], [], []
        for c in common:
            mv = meth[c].get(g)
            ev = expr[ens].get(c)
            if mv is None or ev is None:
                continue
            xs.append(mv); ys.append(math.log2(ev + 1.0)); ts.append(c)
        if len(xs) < MIN_TISSUES:
            continue
        r, p = spearmanr(xs, ys)
        rows.append({"gene": g, "module": info["module"], "n_tissues": len(xs),
                     "spearman_r": round(r, 3), "p_value": round(p, 4),
                     "direction": "negative(repressive)" if r < 0 else "positive",
                     "high_tier": "yes" if g in HIGH_TIER else "no"})
    rows.sort(key=lambda x: (x["high_tier"] != "yes", x["spearman_r"]))

    with (OUT / "matched_correlation.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    tested = len(rows)
    neg = sum(1 for r in rows if r["spearman_r"] < 0)
    neg_sig = sum(1 for r in rows if r["spearman_r"] < 0 and r["p_value"] < 0.05)
    ht = [r for r in rows if r["high_tier"] == "yes"]
    ht_neg = sum(1 for r in ht if r["spearman_r"] < 0)

    summ = [
        "VALIDAÇÃO PAREADA (atlas Zhou/Fang — mesma vaca Holstein 3842, entre tecidos)",
        "=" * 74,
        f"Tecidos pareados WGBS∩RNA-seq: {len(common)} ({', '.join(common)})",
        f"Genes do painel testados (≥{MIN_TISSUES} tecidos): {tested}",
        f"  com correlação NEGATIVA (repressiva, esperada): {neg}/{tested} "
        f"({100*neg/tested:.0f}%)",
        f"  negativa e significativa (p<0,05): {neg_sig}/{tested}",
        f"Candidatos high-tier testados: {len(ht)}; com correlação negativa: {ht_neg}/{len(ht)}",
        "",
        "High-tier (gene: r, p, n_tecidos):",
    ]
    for r in ht:
        summ.append(f"  {r['gene']:7s} r={r['spearman_r']:+.2f} p={r['p_value']} "
                    f"(n={r['n_tissues']})")
    (OUT / "matched_summary.txt").write_text("\n".join(summ) + "\n", encoding="utf-8")
    print("\n".join(summ))

    # figura: histograma das correlações + destaque high-tier
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    allr = [r["spearman_r"] for r in rows]
    ax.hist(allr, bins=12, color="#9bbcd8", edgecolor="white")
    for r in ht:
        ax.axvline(r["spearman_r"], color="#1e8449", lw=1.4)
        ax.text(r["spearman_r"], ax.get_ylim()[1]*0.9, r["gene"], rotation=90,
                fontsize=7.5, va="top", ha="right", color="#1e6e3a")
    ax.axvline(0, color="#444", lw=1)
    ax.set_xlabel("Spearman r — promoter methylation vs expression (across tissues, same cow)")
    ax.set_ylabel("Panel genes")
    ax.set_title(f"Matched-sample validation (Holstein cow {ANIMAL}, {len(common)} tissues)\n"
                 f"{neg}/{tested} genes show the expected negative (repressive) correlation",
                 fontsize=10)
    fig.tight_layout(); fig.savefig(OUT / "matched_correlation.png", dpi=200); plt.close(fig)
    print(f"\n[atlas] saídas em {OUT}")


if __name__ == "__main__":
    main()
