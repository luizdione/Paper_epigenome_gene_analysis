#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
snps_in_cpg.py — Análise de SNPs divergentes Angus vs Nelore em sítios CpG.

OBJETIVO:
  Para cada SNP nos genes-alvo do subprojeto:
    1) Calcular freq_alelo_ALT em Angus e em Nelore (do 1000 Bull Genomes)
    2) Calcular Δfreq = |freq_Angus − freq_Nelore|
    3) Filtrar SNPs com Δfreq > limite (default 0.3)
    4) Cruzar com:
       - DMCs da Fase 1 (sítios diferencialmente metilados)
       - Motifs da Fase 4 v2 (HSF1, NFkB, SP1, etc.)
       - Features genômicas (promotor, exon, intron, intergenic)
    5) Gerar tabela ranqueada de SNPs candidatos a meQTL populacional

INPUT:
  --json                JSON dos genes-alvo do subprojeto
  --vcf                 VCF unificado (com genótipos por amostra)
  --samples-csv         CSV com colunas: sample_id, breed (Angus/Nelore/...)
  --dmcs                CSV de DMCs significativos (chr,start,end,...)
  --motifs              CSV de motifs (Fase 4 v2)
  --reference           FASTA do genoma de referência (para detectar CpG sites)
  --out                 Diretório de saída

OUTPUT:
  - {SUB}_snp_freq_diff.csv      (todos os SNPs com Δfreq calculada)
  - {SUB}_snp_candidates.csv     (filtrados com Δfreq > threshold)
  - {SUB}_snps_in_dmcs.csv       (SNPs que coincidem com DMCs)
  - {SUB}_snps_in_motifs.csv     (SNPs em sítios de motifs)
  - {SUB}_snps_summary.json      (resumo estatístico para a interface)

USO:
  python3 snps_in_cpg.py \
      --json subprojeto/SUB01/SUB01_target_genes_coordinates.json \
      --vcf subprojeto/SUB01/SUB01_variants.vcf.gz \
      --samples-csv subprojeto/SUB01/1000bg_samples_breed.csv \
      --dmcs results/SUB01_GSE161113_run04/fase1_methylkit/DMCs_significant.csv \
      --motifs results/SUB01_GSE161113_run04/fase4_homer/motifs_summary.csv \
      --reference /data/genome_reference/ARS-UCD1.3/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa \
      --out results/SUB01_GSE161113_run04/fase_snps_cpg \
      --delta-threshold 0.3
"""
from __future__ import annotations
import argparse
import csv
import datetime
import json
import os
import re
import sys
from pathlib import Path

# ============================================================
# CLI
# ============================================================
ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                             description=__doc__)
ap.add_argument("--json", required=True, help="JSON dos genes-alvo do subprojeto")
ap.add_argument("--vcf", required=True, help="VCF unificado (.vcf.gz indexado)")
ap.add_argument("--samples-csv", required=False, default=None,
                help="CSV com sample_id,breed do 1000bg. Se ausente, usa breed_field do VCF header.")
ap.add_argument("--dmcs", required=False, default=None,
                help="CSV de DMCs significativos (Fase 1)")
ap.add_argument("--motifs", required=False, default=None,
                help="CSV de motifs (Fase 4 v2)")
ap.add_argument("--reference", required=False, default=None,
                help="FASTA genoma referência (para detectar CpG sites)")
ap.add_argument("--breed-1", default="Angus",
                help="Nome da raça 1 (taurus, default: Angus)")
ap.add_argument("--breed-2", default="Nelore",
                help="Nome da raça 2 (indicus, default: Nelore)")
ap.add_argument("--delta-threshold", type=float, default=0.3,
                help="Limite mínimo de Δfreq (default: 0.3)")
ap.add_argument("--maf-min", type=float, default=0.05,
                help="MAF mínima para considerar SNP (default: 0.05)")
ap.add_argument("--cpg-window", type=int, default=2,
                help="Janela ±N pb para considerar SNP em CpG site (default: 2)")
ap.add_argument("--out", required=True, help="Diretório de saída")
args = ap.parse_args()

OUT = Path(args.out)
OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# Importar pysam (opcional — fallback se não disponível)
# ============================================================
try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False
    print("[WARN] pysam não disponível — SNPs do VCF serão lidos linearmente (lento).")

# ============================================================
# HELPERS
# ============================================================
def log(msg):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {msg}", flush=True)


def norm_chrom(c):
    return re.sub(r"^[Cc][Hh][Rr]", "", str(c))


def load_genes(json_path):
    """Carrega coordenadas dos genes-alvo do subprojeto."""
    with open(json_path) as f:
        data = json.load(f)
    genes = data.get("genes", {})
    if isinstance(genes, list):
        genes = {g.get("gene_name", f"g{i}"): g for i, g in enumerate(genes)}
    out = []
    for name, info in genes.items():
        out.append({
            "gene": name,
            "chrom": norm_chrom(info.get("chromosome", "")),
            "start": int(info.get("start", 0)),
            "end": int(info.get("end", 0)),
            "strand": info.get("strand", 1),
            "module": info.get("module", "?"),
        })
    return out


def load_dmcs(dmcs_csv):
    """Carrega DMCs significativos."""
    if dmcs_csv is None or not os.path.exists(dmcs_csv):
        return []
    out = []
    with open(dmcs_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                out.append({
                    "chrom": norm_chrom(row.get("chr", row.get("chrom", ""))),
                    "start": int(row.get("start", 0)),
                    "end": int(row.get("end", row.get("start", 0))),
                    "meth_diff": float(row.get("meth.diff", 0) or 0),
                    "qvalue": float(row.get("qvalue", 1) or 1),
                })
            except (ValueError, TypeError):
                continue
    return out


def load_motifs(motifs_csv):
    """Carrega lista de motifs enriquecidos da Fase 4 v2."""
    if motifs_csv is None or not os.path.exists(motifs_csv):
        return []
    out = []
    with open(motifs_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                enr = row.get("enrichment", "0")
                if isinstance(enr, str):
                    enr = enr.replace("x", "").strip() or "0"
                out.append({
                    "motif": row.get("motif", "?"),
                    "consensus": row.get("consensus", ""),
                    "enrichment": float(enr) if enr else 0.0,
                })
            except (ValueError, TypeError):
                continue
    return [m for m in out if m["enrichment"] >= 1.5]  # apenas motifs enriquecidos


def load_samples_breed(samples_csv):
    """Mapeia sample_id → breed do 1000bg."""
    if samples_csv is None or not os.path.exists(samples_csv):
        return None
    mapping = {}
    with open(samples_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sid = row.get("sample_id") or row.get("sampleID") or row.get("sample")
            br = row.get("breed") or row.get("Breed")
            if sid and br:
                mapping[sid.strip()] = br.strip()
    return mapping


def calc_allele_freq(genotypes, breed_samples):
    """
    Calcula frequência de alelo ALT entre amostras de uma raça.
    genotypes: dict {sample_id: '0/0' / '0/1' / '1/1' / './.'}
    breed_samples: lista de sample_ids dessa raça
    Retorna: (freq_alt, n_called, n_total)
    """
    n_total = len(breed_samples)
    n_alt = 0
    n_called = 0
    for s in breed_samples:
        gt = genotypes.get(s, "./.")
        if gt in ("./.", ".", ""):
            continue
        # Aceita 0/0, 0|0, 0/1, 1|1 etc.
        alleles = re.split(r"[|/]", gt)[:2]
        if len(alleles) < 2:
            continue
        try:
            a1, a2 = int(alleles[0]), int(alleles[1])
            n_called += 1
            n_alt += a1 + a2
        except (ValueError, TypeError):
            continue
    if n_called == 0:
        return None, 0, n_total
    return n_alt / (2 * n_called), n_called, n_total


def detect_cpg_at_position(ref_fa, chrom, pos, window=2):
    """Verifica se posição está em CpG site (CG ou CpG no ±window)."""
    if ref_fa is None:
        return False
    try:
        # Buscar nome real do contig (com ou sem prefixo chr)
        contigs = set(ref_fa.references)
        chrom_in_fa = None
        for cand in [chrom, "chr" + chrom, "Chr" + chrom]:
            if cand in contigs:
                chrom_in_fa = cand
                break
        if chrom_in_fa is None:
            return False
        # Pegar ±window pb ao redor da posição
        start = max(1, pos - window)
        end = pos + window
        seq = ref_fa.fetch(chrom_in_fa, start - 1, end).upper()
        return "CG" in seq
    except Exception:
        return False


def overlap_dmc(snp_chrom, snp_pos, dmcs):
    """Retorna o primeiro DMC que sobrepõe o SNP (±100 pb)."""
    for d in dmcs:
        if d["chrom"] == snp_chrom:
            dist = abs(snp_pos - d["start"])
            if dist <= 100:
                return d
    return None


def in_motif_consensus(seq_around_snp, motifs):
    """Retorna lista de motifs cujo consensus está em ±N pb do SNP."""
    iupac = {"A": "A", "C": "C", "G": "G", "T": "T",
             "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
             "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
             "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
    hits = []
    for m in motifs:
        try:
            pattern = "".join(iupac.get(c.upper(), c) for c in m["consensus"])
            if re.search(pattern, seq_around_snp, re.IGNORECASE):
                hits.append(m["motif"])
        except Exception:
            continue
    return hits


def feature_at_pos(snp_chrom, snp_pos, gene):
    """Categoriza posição relativa ao gene: promoter (TSS±2kb), gene_body, upstream, downstream."""
    g_start = gene["start"]
    g_end = gene["end"]
    strand = gene["strand"]
    if strand in (1, "+", "+1"):
        tss = g_start
        if abs(snp_pos - tss) <= 2000:
            return "promoter"
        elif g_start <= snp_pos <= g_end:
            return "gene_body"
        elif snp_pos < tss - 2000:
            return "upstream_distal"
        else:
            return "downstream"
    else:
        tss = g_end
        if abs(snp_pos - tss) <= 2000:
            return "promoter"
        elif g_start <= snp_pos <= g_end:
            return "gene_body"
        elif snp_pos > tss + 2000:
            return "upstream_distal"
        else:
            return "downstream"


# ============================================================
# MAIN
# ============================================================
def main():
    log("=" * 70)
    log(f"  SNPs em CpG — Análise de divergência {args.breed_1} vs {args.breed_2}")
    log("=" * 70)

    # Carregar genes-alvo
    genes = load_genes(args.json)
    log(f"Genes-alvo carregados: {len(genes)}")

    # Carregar DMCs e motifs
    dmcs = load_dmcs(args.dmcs)
    log(f"DMCs significativos: {len(dmcs)}")

    motifs = load_motifs(args.motifs)
    log(f"Motifs enriquecidos (≥1.5x): {len(motifs)} → {[m['motif'] for m in motifs]}")

    # Carregar mapping de breeds
    sample_breed = load_samples_breed(args.samples_csv)
    if sample_breed is None:
        log(f"[ERRO] CSV de samples-breed obrigatório quando o VCF não tem header anotado.")
        log(f"[INFO] Crie um CSV com colunas: sample_id, breed")
        sys.exit(1)

    breeds_set = set(sample_breed.values())
    log(f"Breeds no CSV: {sorted(breeds_set)}")
    samples_b1 = [s for s, b in sample_breed.items() if b == args.breed_1]
    samples_b2 = [s for s, b in sample_breed.items() if b == args.breed_2]
    log(f"  {args.breed_1}: {len(samples_b1)} amostras")
    log(f"  {args.breed_2}: {len(samples_b2)} amostras")

    # Abrir genoma de referência (para CpG e contexto motif)
    ref_fa = None
    if args.reference:
        if HAS_PYSAM:
            ref_fa = pysam.FastaFile(args.reference)
            log(f"Genoma de referência aberto: {args.reference}")
        else:
            log(f"[WARN] pysam ausente — análise de CpG e contexto motif desabilitada.")

    # Abrir VCF
    if not HAS_PYSAM:
        log("[ERRO] pysam é obrigatório para ler VCF.")
        sys.exit(1)

    log(f"Abrindo VCF: {args.vcf}")
    vcf = pysam.VariantFile(args.vcf)
    vcf_samples = list(vcf.header.samples)
    log(f"  Amostras no VCF: {len(vcf_samples)}")

    # Filtrar samples_b1/b2 pelos que existem no VCF
    samples_b1 = [s for s in samples_b1 if s in vcf_samples]
    samples_b2 = [s for s in samples_b2 if s in vcf_samples]
    log(f"  Após filtragem por VCF: {args.breed_1}={len(samples_b1)}, {args.breed_2}={len(samples_b2)}")
    if len(samples_b1) < 5 or len(samples_b2) < 5:
        log("[ERRO] Pelo menos 5 amostras de cada raça são necessárias.")
        sys.exit(1)

    # ── Iterar SNPs nos genes-alvo ──────────────────────────────
    log("Processando SNPs em cada gene-alvo...")
    rows_all = []
    rows_candidates = []
    rows_dmc = []
    rows_motif = []

    for gene in genes:
        chrom = gene["chrom"]
        start_q = max(1, gene["start"] - 500)
        end_q = gene["end"] + 500
        # Match contig name
        contigs = list(vcf.header.contigs)
        chrom_vcf = None
        for c in [chrom, "chr" + chrom, "Chr" + chrom]:
            if c in contigs:
                chrom_vcf = c
                break
        if chrom_vcf is None:
            log(f"  [SKIP] {gene['gene']}: chr{chrom} não no VCF")
            continue

        try:
            recs = vcf.fetch(chrom_vcf, start_q - 1, end_q)
        except Exception as e:
            log(f"  [WARN] {gene['gene']} ({chrom_vcf}:{start_q}-{end_q}): {e}")
            continue

        n_snps = 0
        n_candidates = 0
        for rec in recs:
            if rec.alts is None or len(rec.alts) == 0:
                continue
            # Apenas SNPs bialélicos
            ref = rec.ref
            alt = rec.alts[0]
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                continue
            n_snps += 1

            # Extrair genótipos para amostras de cada raça
            gt_b1 = {}
            gt_b2 = {}
            for sample, sample_data in rec.samples.items():
                gt = sample_data.get("GT", (None, None))
                if gt is None or any(g is None for g in gt):
                    gt_str = "./."
                else:
                    gt_str = "/".join(str(g) for g in gt)
                if sample in samples_b1:
                    gt_b1[sample] = gt_str
                elif sample in samples_b2:
                    gt_b2[sample] = gt_str

            f_b1, n_called_b1, _ = calc_allele_freq(gt_b1, samples_b1)
            f_b2, n_called_b2, _ = calc_allele_freq(gt_b2, samples_b2)

            if f_b1 is None or f_b2 is None:
                continue
            if n_called_b1 < 5 or n_called_b2 < 5:
                continue

            delta = abs(f_b1 - f_b2)
            maf = min(f_b1, 1 - f_b1, f_b2, 1 - f_b2)
            if maf < args.maf_min:
                continue

            # Verificar CpG
            is_cpg = detect_cpg_at_position(ref_fa, chrom, rec.pos, args.cpg_window) if ref_fa else False

            # Sobreposição com DMC
            dmc_hit = overlap_dmc(chrom, rec.pos, dmcs)

            # Contexto local da sequência (±15 pb) para motif scan
            motif_hits = []
            if ref_fa is not None:
                try:
                    contigs_fa = set(ref_fa.references)
                    fa_chrom = None
                    for c in [chrom, "chr" + chrom, "Chr" + chrom]:
                        if c in contigs_fa:
                            fa_chrom = c
                            break
                    if fa_chrom:
                        seq_local = ref_fa.fetch(fa_chrom, max(0, rec.pos - 15), rec.pos + 15).upper()
                        motif_hits = in_motif_consensus(seq_local, motifs)
                except Exception:
                    pass

            # Feature
            feat = feature_at_pos(chrom, rec.pos, gene)

            row = {
                "gene": gene["gene"],
                "module": gene["module"],
                "chrom": chrom,
                "pos": rec.pos,
                "ref": ref,
                "alt": alt,
                "rsid": rec.id or ".",
                f"freq_alt_{args.breed_1}": round(f_b1, 4),
                f"freq_alt_{args.breed_2}": round(f_b2, 4),
                f"n_called_{args.breed_1}": n_called_b1,
                f"n_called_{args.breed_2}": n_called_b2,
                "delta_freq": round(delta, 4),
                "maf": round(maf, 4),
                "feature": feat,
                "is_cpg": is_cpg,
                "dmc_overlap": "yes" if dmc_hit else "no",
                "dmc_meth_diff": round(dmc_hit["meth_diff"], 2) if dmc_hit else None,
                "motifs_in_15bp": ",".join(motif_hits) if motif_hits else "",
            }
            rows_all.append(row)

            if delta >= args.delta_threshold:
                rows_candidates.append(row)
                n_candidates += 1
                if dmc_hit:
                    rows_dmc.append(row)
                if motif_hits:
                    rows_motif.append(row)

        if n_snps > 0:
            log(f"  {gene['gene']}: {n_snps} SNPs avaliados, {n_candidates} candidatos (Δfreq>{args.delta_threshold})")

    vcf.close()
    if ref_fa:
        ref_fa.close()

    # ── Salvar CSVs ─────────────────────────────────────────────
    sub_id = Path(args.json).stem.split("_")[0]  # ex: SUB01
    log("")
    log(f"Salvando outputs em {OUT}")

    if rows_all:
        out_csv = OUT / f"{sub_id}_snp_freq_diff.csv"
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows_all[0].keys()))
            w.writeheader()
            w.writerows(rows_all)
        log(f"  [OK] Todos os SNPs: {out_csv} ({len(rows_all)} linhas)")

    if rows_candidates:
        out_csv = OUT / f"{sub_id}_snp_candidates.csv"
        rows_candidates.sort(key=lambda r: -r["delta_freq"])
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows_candidates[0].keys()))
            w.writeheader()
            w.writerows(rows_candidates)
        log(f"  [OK] Candidatos (Δfreq>{args.delta_threshold}): {out_csv} ({len(rows_candidates)} linhas)")

    if rows_dmc:
        out_csv = OUT / f"{sub_id}_snps_in_dmcs.csv"
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows_dmc[0].keys()))
            w.writeheader()
            w.writerows(rows_dmc)
        log(f"  [OK] SNPs em DMCs: {out_csv} ({len(rows_dmc)} linhas)")

    if rows_motif:
        out_csv = OUT / f"{sub_id}_snps_in_motifs.csv"
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows_motif[0].keys()))
            w.writeheader()
            w.writerows(rows_motif)
        log(f"  [OK] SNPs em motifs: {out_csv} ({len(rows_motif)} linhas)")

    # ── Resumo JSON ─────────────────────────────────────────────
    summary = {
        "subproject": sub_id,
        "date": datetime.datetime.now().isoformat(),
        "breeds": {args.breed_1: len(samples_b1), args.breed_2: len(samples_b2)},
        "thresholds": {
            "delta_freq": args.delta_threshold,
            "maf_min": args.maf_min,
            "cpg_window": args.cpg_window,
        },
        "n_snps_evaluated": len(rows_all),
        "n_candidates": len(rows_candidates),
        "n_in_dmcs": len(rows_dmc),
        "n_in_motifs": len(rows_motif),
        "n_cpg_sites": sum(1 for r in rows_all if r["is_cpg"]),
        "top_candidates": [
            {k: r[k] for k in ("gene", "module", "chrom", "pos", "ref", "alt",
                               "delta_freq", "feature", "is_cpg", "dmc_overlap",
                               "motifs_in_15bp")}
            for r in (rows_candidates[:20] if rows_candidates else [])
        ],
        "modules_affected": sorted(set(r["module"] for r in rows_candidates)),
    }
    with open(OUT / f"{sub_id}_snps_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    log(f"  [OK] Resumo: {OUT / f'{sub_id}_snps_summary.json'}")

    log("")
    log("=" * 70)
    log("  CONCLUÍDO")
    log("=" * 70)
    log(f"  SNPs avaliados:    {len(rows_all):,}")
    log(f"  Candidatos:        {len(rows_candidates):,} (Δfreq > {args.delta_threshold})")
    log(f"  SNPs em DMCs:      {len(rows_dmc):,}")
    log(f"  SNPs em motifs:    {len(rows_motif):,}")
    log(f"  SNPs em CpG:       {summary['n_cpg_sites']:,}")
    log(f"  Outputs em:        {OUT}")


if __name__ == "__main__":
    main()
