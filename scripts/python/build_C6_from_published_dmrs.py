#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_C6_from_published_dmrs.py — adiciona o contraste C6 (Brahman × Angus,
fígado fetal, WGBS) ao SUB02 a partir dos DMRs PROCESSADOS e PUBLICADOS por
MacPhillamy et al. 2024 (GigaScience giae061; repo DaviesCentreInformatics/
Brahman_Angus_WGBS, DMRs/AAvsBB_DMRs_10pALL.bed, chamados com methylKit,
alinhados ao Angus/ARS-UCD).

NÃO reprocessa WGBS: usa o resultado publicado (meth.diff Angus−Brahman, qvalue)
e o cruza com as janelas dos genes-alvo do SUB02, escrevendo a anotação no MESMO
formato dos demais contrastes para a ponte multi-ômica incorporar o C6.

Filtros (iguais aos demais contrastes): qvalue < 0.05 e |meth.diff| > 15% (e 25%).
"""
from __future__ import annotations
import csv, json
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
JSON = BASE / "subprojeto" / "SUB02" / "SUB02_target_genes_coordinates.json"
COMP = BASE / "results" / "SUB02_runs" / "C6_GIGAS_BIBI_vs_BTBT_fetalLiver"
DMR = COMP / "source" / "AAvsBB_DMRs_10pALL.bed"


def load_genes():
    d = json.loads(JSON.read_text(encoding="utf-8"))["genes"]
    genes = []
    for g, info in d.items():
        try:
            chrom = str(info["chromosome"])
            start = int(info["start"]); end = int(info["end"])
            strand = int(info.get("strand", 1))
            fs = int(info.get("flanking_start", start - 10000))
            fe = int(info.get("flanking_end", end + 10000))
            genes.append((g, chrom, start, end, strand, fs, fe,
                          info.get("module", "outros")))
        except (KeyError, ValueError, TypeError):
            continue
    return genes


def feature_bucket(mid, start, end, strand):
    tss = start if strand == 1 else end
    if start <= mid <= end:
        return "gene_body"
    # distância ao TSS no sentido upstream (negativa = upstream)
    d = (mid - tss) if strand == 1 else (tss - mid)
    if d < 0:
        return "promoter" if abs(d) <= 1000 else "upstream_distal"
    return "downstream_distal"


def main():
    genes = load_genes()
    # indexa DMRs por cromossomo (apenas significativos p/ eficiência)
    sig15, sig25 = [], []
    with DMR.open(encoding="utf-8") as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            try:
                chrom = p[0]; s = int(p[1]); e = int(p[2])
                md = float(p[3]); q = float(p[4])
            except ValueError:
                continue
            if q < 0.05 and abs(md) > 15:
                sig15.append((chrom, s, e, md, q))
                if abs(md) > 25:
                    sig25.append((chrom, s, e, md, q))
    by_chr = {}
    for rec in sig15:
        by_chr.setdefault(rec[0], []).append(rec)

    # cruza com janelas dos genes
    rows15, rows25 = [], []
    for g, chrom, gs, ge, strand, fs, fe, module in genes:
        for (c, s, e, md, q) in by_chr.get(chrom, []):
            if e >= fs and s <= fe:            # sobreposição com a janela ±10 kb
                mid = (s + e) // 2
                feat = feature_bucket(mid, gs, ge, strand)
                row = {"chr": "chr" + chrom, "start": s, "end": e, "strand": "*",
                       "pvalue": "", "qvalue": q, "meth.diff": md,
                       "feature": feat, "gene": g}
                rows15.append(row)
                if abs(md) > 25:
                    rows25.append(row)

    for diff, rows in (("15", rows15), ("25", rows25)):
        outdir = COMP / f"fase2_annotation_diff{diff}"
        outdir.mkdir(parents=True, exist_ok=True)
        cols = ["chr", "start", "end", "strand", "pvalue", "qvalue",
                "meth.diff", "feature", "gene"]
        with (outdir / "annotation_table.csv").open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, quoting=csv.QUOTE_NONNUMERIC)
            w.writeheader(); w.writerows(rows)

    genes_hit = sorted({r["gene"] for r in rows15})
    print(f"[C6] DMRs significativos (q<0.05,|diff|>15%): {len(sig15)} genome-wide; "
          f">25%: {len(sig25)}")
    print(f"[C6] sobreposições DMR×gene-alvo (15%): {len(rows15)} em {len(genes_hit)} genes")
    print(f"[C6] genes: {', '.join(genes_hit)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
