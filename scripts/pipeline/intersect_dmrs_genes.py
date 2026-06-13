#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fase 2 — Interseção Genômica (BEDTools substituto puro Python)
=============================================================
Cruza uma tabela de DMRs/DMCs (saída do methylKit) com as coordenadas de
regiões regulatórias alvo (flanco ±10kb) definidas em
`all_target_genes_coordinates.json`.

Vantagens desta implementação:
    - Não depende de `bedtools` instalado (opera offline em qualquer SO).
    - Aceita DMRs em CSV, TSV ou BED.
    - Usa índice por cromossomo + interval tree simples para performance.
    - Preserva o nome da coluna original da DMR e adiciona:
        gene, module, confidence, ncbi_gene_id, dist_to_TSS

Uso:
    python intersect_dmrs_genes.py \
        --dmrs results/diff_meth.csv \
        --json all_target_genes_coordinates.json \
        --out results/dmrs_in_target_genes.csv

Colunas esperadas na DMR (case-insensitive; qualquer um funciona):
    chr | chrom | chromosome        -> cromossomo
    start | pos | position          -> posição inicial
    end                             -> posição final (opcional; se ausente, end = start)
    meth.diff | methdiff            -> diferença de metilação (opcional)
    qvalue | q.value | padj         -> q-value (opcional)

Autor: Luiz Dione Barbosa De-Melo — GenOmics Epigenoma 2026
"""
from __future__ import annotations
import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from typing import Dict, List, Tuple


# --------------------------------------------------------------------- helpers
def norm_chrom(c) -> str:
    """Normaliza cromossomo removendo 'chr' e espaços."""
    if c is None:
        return ""
    return str(c).strip().lower().removeprefix("chr")


def pick_col(header: List[str], candidates: List[str]) -> int | None:
    """Retorna o índice da primeira coluna em `header` que casa com algum alias."""
    lower = [h.strip().lower() for h in header]
    for cand in candidates:
        if cand in lower:
            return lower.index(cand)
    return None


# --------------------------------------------------------------------- core
def load_gene_regions(json_path: str) -> Dict[str, List[Tuple[int, int, dict]]]:
    """Carrega o JSON de alvos e retorna um dict chrom -> lista de (start, end, info)."""
    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    genes = data.get("genes", {})
    idx: Dict[str, List[Tuple[int, int, dict]]] = defaultdict(list)
    for name, info in genes.items():
        chrom = norm_chrom(info.get("chromosome"))
        start = info.get("flanking_start") or info.get("start")
        end = info.get("flanking_end") or info.get("end")
        if not chrom or start is None or end is None:
            continue
        payload = {
            "gene": name,
            "module": info.get("module", ""),
            "confidence": info.get("confidence", ""),
            "ensembl_id": info.get("ensembl_id", ""),
            "ncbi_gene_id": info.get("ncbi_gene_id", ""),
            "strand": info.get("strand", 1),
            "gene_start": info.get("start"),
            "gene_end": info.get("end"),
        }
        idx[chrom].append((int(start), int(end), payload))
    # ordena por start para busca mais previsível
    for c in idx:
        idx[c].sort(key=lambda t: t[0])
    return idx


def sniff_delimiter(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        sample = fh.read(2048)
    if "\t" in sample and sample.count("\t") > sample.count(","):
        return "\t"
    return ","


def dist_to_tss(pos: int, payload: dict) -> int:
    """Distância do ponto ao TSS do gene (considerando a fita)."""
    gs, ge = payload.get("gene_start"), payload.get("gene_end")
    strand = payload.get("strand", 1) or 1
    if gs is None or ge is None:
        return 0
    tss = gs if strand == 1 else ge
    return int(pos) - int(tss)


def intersect(dmrs_path: str, json_path: str, out_path: str) -> dict:
    regions = load_gene_regions(json_path)
    if not regions:
        raise SystemExit(f"[ERRO] Nenhuma região carregada de {json_path}")

    delim = sniff_delimiter(dmrs_path)
    stats = {"input_rows": 0, "hits": 0, "genes_hit": set(), "modules_hit": set()}

    with open(dmrs_path, "r", encoding="utf-8", errors="ignore", newline="") as fin:
        reader = csv.reader(fin, delimiter=delim)
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"[ERRO] Arquivo de DMRs vazio: {dmrs_path}")

        ic = pick_col(header, ["chr", "chrom", "chromosome", "#chrom"])
        istart = pick_col(header, ["start", "pos", "position"])
        iend = pick_col(header, ["end", "stop"])
        idiff = pick_col(header, ["meth.diff", "methdiff", "meth_diff", "methylation_difference", "diff"])
        iq = pick_col(header, ["qvalue", "q.value", "q_value", "padj", "fdr"])

        if ic is None or istart is None:
            raise SystemExit(
                f"[ERRO] Não foi possível identificar colunas de chr/start no header: {header}"
            )

        os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
        with open(out_path, "w", encoding="utf-8", newline="") as fout:
            writer = csv.writer(fout)
            writer.writerow(
                [
                    "chr", "start", "end",
                    "meth_diff", "qvalue",
                    "gene", "module", "confidence",
                    "ensembl_id", "ncbi_gene_id",
                    "dist_to_TSS",
                ]
            )
            for row in reader:
                if not row or len(row) <= ic:
                    continue
                stats["input_rows"] += 1
                chrom = norm_chrom(row[ic])
                try:
                    start = int(float(row[istart]))
                except (ValueError, IndexError):
                    continue
                end = start
                if iend is not None and iend < len(row):
                    try:
                        end = int(float(row[iend]))
                    except ValueError:
                        pass
                diff = row[idiff] if idiff is not None and idiff < len(row) else ""
                qval = row[iq] if iq is not None and iq < len(row) else ""

                for gs, ge, payload in regions.get(chrom, []):
                    # fim da região abaixo do começo do DMR -> sem chance (lista ordenada)
                    if ge < start:
                        continue
                    # início da região acima do fim do DMR -> já passamos
                    if gs > end:
                        break
                    # sobreposição confirmada
                    writer.writerow(
                        [
                            chrom, start, end, diff, qval,
                            payload["gene"], payload["module"], payload["confidence"],
                            payload["ensembl_id"], payload["ncbi_gene_id"],
                            dist_to_tss((start + end) // 2, payload),
                        ]
                    )
                    stats["hits"] += 1
                    stats["genes_hit"].add(payload["gene"])
                    stats["modules_hit"].add(payload["module"])

    stats["genes_hit"] = sorted(stats["genes_hit"])
    stats["modules_hit"] = sorted(stats["modules_hit"])
    return stats


# --------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dmrs", required=True, help="CSV/TSV/BED com DMRs (saída methylKit)")
    ap.add_argument("--json", default="all_target_genes_coordinates.json", help="JSON de genes alvo")
    ap.add_argument("--out", default="results/dmrs_in_target_genes.csv", help="CSV de saída")
    args = ap.parse_args()

    if not os.path.exists(args.dmrs):
        sys.exit(f"[ERRO] DMRs não encontrado: {args.dmrs}")
    if not os.path.exists(args.json):
        sys.exit(f"[ERRO] JSON não encontrado: {args.json}")

    print(f"[Info] Interseção: {args.dmrs}  ×  {args.json}")
    stats = intersect(args.dmrs, args.json, args.out)

    print(f"[Info] Linhas de entrada processadas: {stats['input_rows']}")
    print(f"[Info] Sobreposições encontradas : {stats['hits']}")
    print(f"[Info] Genes atingidos           : {len(stats['genes_hit'])}")
    if stats["genes_hit"]:
        print("       -> " + ", ".join(stats["genes_hit"][:20]) + (" ..." if len(stats["genes_hit"]) > 20 else ""))
    print(f"[Info] Módulos atingidos         : {', '.join(stats['modules_hit']) or '(nenhum)'}")
    print(f"[OK ] Resultado salvo em: {args.out}")


if __name__ == "__main__":
    main()
