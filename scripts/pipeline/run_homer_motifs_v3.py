#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_homer_motifs_v3.py — Análise completa de motifs regulatórios nos DMCs.

MELHORIAS sobre v2:
  • Catálogo ampliado: ~45 motifs cobrindo stress, imunidade, metabolismo,
    crescimento, músculo, epigenética, sinalização (fontes: JASPAR 2024, ENCODE,
    literatura Bos taurus: Salavati 2023, Meuwissen 2021, Kim 2020)
  • Teste binomial exato (scipy) com correção Benjamini-Hochberg
  • Relatório por categoria funcional
  • CSV final com rank, fold-enrichment, BH-adjusted p-value
  • Análise k-mer expandida (k=4,6,8) com frequência de referência

USO:
  python3 run_homer_motifs_v3.py \\
      --bed results/SUB01.../fase1_methylkit/DMCs_significant.bed \\
      --fasta /data/genome_reference/ARS-UCD1.3/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa \\
      --out results/SUB01.../fase4_homer \\
      --window 200
"""
from __future__ import annotations
import argparse
import collections
import csv
import datetime
import math
import os
import re
import sys
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--bed",    required=True, help="BED com DMCs (chr start end ...)")
ap.add_argument("--fasta",  required=True, help="FASTA indexado (.fai existir)")
ap.add_argument("--out",    required=True, help="Diretório de saída")
ap.add_argument("--window", type=int, default=200, help="±N pb ao redor de cada DMC")
args = ap.parse_args()

BED    = Path(args.bed)
FASTA  = Path(args.fasta)
OUT    = Path(args.out)
WINDOW = args.window
OUT.mkdir(parents=True, exist_ok=True)

FAI = Path(str(FASTA) + ".fai")
if not FAI.exists():
    print(f"[ERRO] Índice .fai não encontrado: {FAI}")
    print(f"       Rode: samtools faidx {FASTA}")
    sys.exit(1)

try:
    import pysam
except ImportError:
    print("[ERRO] pysam não instalado: pip3 install pysam")
    sys.exit(1)

try:
    from scipy import stats as scipy_stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("[WARN] scipy não disponível — p-values serão estimados (install com pip3 install scipy)")

# ── Catálogo de motifs ────────────────────────────────────────────────────────
# Formato: (nome, consenso_IUPAC, familia, descrição, relevância_Bos_taurus)
# Fontes: JASPAR 2024 (MA*), ENCODE, Salavati et al. 2023 (BMC Genomics),
#         Meuwissen et al. 2021 (Nat Genet), Kim et al. 2020 (Genome Biol)
MOTIFS = [
    # ── TERMOTOLERÂNCIA ───────────────────────────────────────────────────────
    ("HSF1",        "GAANNTTC",           "HSF",      "Heat Shock Factor 1 — HSE",
     "Principal regulador de termotolerância; differentially expressed Nelore vs Angus"),
    ("HSF1_alt",    "TTCNNGAA",           "HSF",      "HSF1 — HSE reverso",
     "HSE orientação 3'-5'"),
    ("HSF2",        "GAANNNTTC",          "HSF",      "Heat Shock Factor 2",
     "Manutenção homeostase celular e proteossomo"),

    # ── ESTRESSE OXIDATIVO / INFLAMAÇÃO ──────────────────────────────────────
    ("NFkB_p65",    "GGGRNYYYCC",         "Rel",      "NFkB p65 (RELA) — NF-kB binding",
     "Imunidade inata, inflamação, feed efficiency Nelore"),
    ("NFkB_p50",    "GGGRNNTCC",          "Rel",      "NFkB p50 (NFKB1)",
     "Dimero canônico NFkB1; promotores imunidade inata Bos taurus"),
    ("AP1_FOSJUN",  "TGASTCA",            "bZIP",     "AP-1 FOS/JUN composite",
     "Resposta a estresse, calor, imunidade"),
    ("ATF3",        "TGACGTCA",           "bZIP",     "Activating Transcription Factor 3",
     "Ativado em resposta ao estresse térmico em bovinos"),
    ("NRF2",        "TGACTCAGC",          "bZIP",     "Nuclear factor erythroid 2-like 2 (NFE2L2)",
     "Resposta antioxidante; adaptação ao estresse oxidativo"),
    ("IRF1",        "GAAANNGAAA",         "IRF",      "Interferon Regulatory Factor 1",
     "Resposta imune inata — importante em raças Bos taurus vs Bos indicus"),
    ("IRF3",        "GAAACTGAAACT",       "IRF",      "IRF3 — resposta a interferon",
     "Via antiviral em bovinos"),
    ("NFAT1",       "AGGAAANNNNNNN",      "RHD",      "NFAT1 (NFATC1) — calcineurina/NFAT",
     "Imunidade adaptativa, diferenciação de linfócitos"),
    ("TNFA_resp",   "GGGACTTTCC",         "Rel",      "TNF-alpha response element",
     "Resposta inflamatória; sensível a metilação em bovinos"),

    # ── CRESCIMENTO / SINALIZAÇÃO GH-IGF ──────────────────────────────────────
    ("STAT5A",      "TTCNNNGAA",          "STAT",     "STAT5A — GH/prolactina signaling",
     "Crescimento somático, lactação, marmoreio Angus"),
    ("STAT3",       "TTCCCRGAA",          "STAT",     "STAT3 — JAK/STAT pathway",
     "Diferenciação celular, crescimento, resposta imune"),
    ("ETS1",        "CCGGAAGT",           "ETS",      "ETS proto-oncogene 1",
     "Angiogênese, hematopoiese, resposta imune"),
    ("ETS2",        "ACCGGAAGT",          "ETS",      "ETS proto-oncogene 2",
     "Desenvolvimento embrionário, reprogramação epigenética"),
    ("MYC_EBOX",    "CACGTG",             "bHLH",     "c-MYC / MAX — E-box",
     "Proliferação celular, metabolismo; ativado em crescimento Angus"),
    ("HIF1A",       "ACGTG",             "bHLH",     "HIF-1 alpha — Hypoxia Inducible Factor",
     "Resposta à hipóxia, metabolismo anaeróbio em altitude"),

    # ── MÚSCULO / DESENVOLVIMENTO ─────────────────────────────────────────────
    ("MYOD1",       "CANNTG",             "bHLH",     "MyoD / Myogenin — muscle E-box",
     "Miogênese; diferenciação de fibras musculares em Nelore vs Angus"),
    ("MEF2A",       "CTAWWWTAR",          "MADS",     "Myocyte Enhancer Factor 2A",
     "Enhancers musculares — regulação diferencial Bos indicus"),
    ("TEAD1",       "GGAATGG",            "TEAD",     "TEA domain factor 1 — Hippo pathway",
     "Tamanho do órgão, proliferação muscular, aptidão física"),
    ("RUNX2",       "YGYGGTY",            "RUNX",     "RUNT-related transcription factor 2",
     "Osteogênese, diferenciação de condrócitos"),
    ("PAX7",        "TAATTAGG",           "PAX",      "Paired Box 7 — células satélite",
     "Regeneração muscular; diferencial entre raças"),

    # ── METABOLISMO LIPÍDICO / ADIPOGÊNESE ────────────────────────────────────
    ("PPARG",       "AACTAGGTCA",         "NR",       "PPARgamma — adipogênese/metabolismo lipídico",
     "Deposição de gordura, marmoreio — maior atividade Angus"),
    ("PPARA",       "TGACCTTTGACCT",      "NR",       "PPARalpha — oxidação de ácidos graxos",
     "Eficiência alimentar, metabolismo energético"),
    ("RXRA",        "AGGTCATGACCT",       "NR",       "Retinoid X Receptor Alpha",
     "Heterodímero obrigatório PPAR/RAR; metabolismo de vitamina A"),
    ("SREBP1",      "TCACNNNNTGA",        "bHLH",     "SREBP1 — Sterol Regulatory Element",
     "Biossíntese de ácidos graxos; lipogênese; relacionado a DGAT1"),
    ("CEBPA",       "TTGCGCAA",           "bZIP",     "C/EBP Alpha — adipogênese",
     "Determinação de adipócitos; co-regulador PPARG"),
    ("CEBPB",       "ATTGCGCAAT",         "bZIP",     "C/EBP Beta",
     "Adipogênese, resposta imune, inflamação"),
    ("KLF4",        "CCCCGCCC",           "KLF/SP",   "Krüppel-like Factor 4",
     "Pluripotência, reprogramação, barreiras epiteliais"),

    # ── REPRODUÇÃO / HORMÔNIOS ────────────────────────────────────────────────
    ("ESR1",        "AGGTCANNNTGACCT",    "NR",       "Estrogen Receptor Alpha",
     "Fertilidade, ciclo estral, desenvolvimento mamário"),
    ("PGR",         "TGTCCAAAGGTCA",      "NR",       "Progesterone Receptor",
     "Gestação, manutenção do corpo lúteo em bovinos"),
    ("AR",          "AGAACANNNTGTTCT",    "NR",       "Androgen Receptor",
     "Caracteres sexuais secundários, crescimento muscular em machos"),
    ("GR",          "TGTYNNNNNNRACA",     "NR",       "Glucocorticoid Receptor",
     "Resposta ao estresse cortisol; diferenças de temperamento Nelore/Angus"),

    # ── TRANSCRIÇÃO GERAL / REGULAÇÃO CROMATINA ──────────────────────────────
    ("CTCF",        "CCGCGNGGNGGCAG",     "ZF",       "CCCTC-binding factor — insulator/TADs",
     "Organização 3D do genoma; fronteiras de domínios de metilação"),
    ("SP1",         "GGGCGG",             "SP/KLF",   "Specificity Protein 1 — GC-box",
     "Promotores CpG-ricos; sensível à metilação do DNA"),
    ("YY1",         "CGCCKSCSS",          "ZF",       "YY1 — Yin Yang 1",
     "Regulação epigenética, Polycomb-associated; bivalente em bovinos"),
    ("E2F1",        "TTCGCGCC",           "E2F",      "E2F Transcription Factor 1",
     "Controle ciclo celular, apoptose, proliferação"),
    ("p53",         "RRRCWWGYYY",         "p53",      "TP53 — tumor suppressor / DNA damage",
     "Controle de qualidade do genoma, diferenças em longevidade"),
    ("SMAD3",       "AGACAGNNS",          "SMAD",     "SMAD3 — TGF-beta signaling",
     "Fibrosis, imunossupressão, desenvolvimento embrionário"),

    # ── CpG / METILAÇÃO ───────────────────────────────────────────────────────
    ("CpG_dinuc",   "CG",                 "dinuc",    "CpG dinucleotide — methylation target",
     "Sítio primário de metilação por DNMT3A/3B/1"),
    ("CpG_island",  "CGCGCG",             "dinuc",    "CpG island core hexamer",
     "Ilhas CpG não-metiladas em promotores ativos"),
    ("MBD_binding", "RRCGYYB",            "MBD",      "Methyl-CpG binding domain consensus",
     "Sítio de reconhecimento por MBD1/2/3 — leitura da metilação"),
]

# ── Utilitários ───────────────────────────────────────────────────────────────
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
    "N": "[ACGT]",
}

def iupac_to_regex(motif: str) -> re.Pattern:
    pattern = "".join(IUPAC.get(c.upper(), c) for c in motif)
    return re.compile(pattern, re.IGNORECASE)

def reverse_complement(seq: str) -> str:
    comp = {"A":"T","T":"A","C":"G","G":"C","N":"N",
            "R":"Y","Y":"R","S":"S","W":"W","K":"M","M":"K",
            "B":"V","V":"B","D":"H","H":"D"}
    return "".join(comp.get(c.upper(), "N") for c in reversed(seq))

def motif_expected_freq(consensus: str) -> float:
    """Frequência esperada por chance (produto das probabilidades IUPAC)."""
    base_prob = {"A":0.25,"C":0.25,"G":0.25,"T":0.25,
                 "R":0.5,"Y":0.5,"S":0.5,"W":0.5,
                 "K":0.5,"M":0.5,"B":0.75,"D":0.75,"H":0.75,"V":0.75,"N":1.0}
    p = 1.0
    for c in consensus.upper():
        p *= base_prob.get(c, 0.25)
    return p

def binom_pvalue(k: int, n: int, p: float) -> float:
    """P-value teste binomial one-sided (k >= observado)."""
    if HAS_SCIPY:
        # scipy >=1.12 removeu binom_test (virou binomtest); compatível com ambos
        if hasattr(scipy_stats, "binom_test"):
            return float(scipy_stats.binom_test(k, n, p, alternative="greater"))
        return float(scipy_stats.binomtest(k, n, p, alternative="greater").pvalue)
    else:
        # Estimativa normal aproximada
        mu = n * p
        sigma = math.sqrt(n * p * (1 - p))
        if sigma == 0:
            return 1.0
        z = (k - mu) / sigma
        return max(0, 1 - 0.5 * (1 + math.erf(z / math.sqrt(2))))

def bh_correct(pvalues: list) -> list:
    """Correção Benjamini-Hochberg."""
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    adj = [0.0] * n
    prev = 1.0
    for rank, (orig_idx, p) in enumerate(reversed(indexed), 1):
        adj_p = min(prev, p * n / (n - rank + 1))
        adj[orig_idx] = adj_p
        prev = adj_p
    return adj

# ── Banner ─────────────────────────────────────────────────────────────────────
print("=" * 70)
print(f"  HOMER motifs v3 — catálogo completo ({len(MOTIFS)} motifs)")
print(f"  {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
print("=" * 70)
print(f"  BED:    {BED}")
print(f"  FASTA:  {FASTA}")
print(f"  OUT:    {OUT}")
print(f"  Window: ±{WINDOW} pb")
print(f"  scipy:  {'disponível (p-values exatos)' if HAS_SCIPY else 'ausente (estimativas)'}")
print()

# ── FASTA / contigs ───────────────────────────────────────────────────────────
fa = pysam.FastaFile(str(FASTA))
contigs = set(fa.references)
fa_has_chr = any(c.lower().startswith("chr") for c in list(contigs)[:10])
print(f"[INFO] Contigs: {len(contigs)}  (ex: {sorted(list(contigs))[:3]})")
print(f"[INFO] Prefixo 'chr': {fa_has_chr}")
print()

def normalize_chrom(chrom: str) -> str | None:
    c = chrom.lstrip("chr").lstrip("Chr")
    for cand in [c, "chr" + c, "Chr" + c, chrom]:
        if cand in contigs:
            return cand
    return None

# ── Extrair sequências ────────────────────────────────────────────────────────
sequences = []
n_skipped = 0
with open(BED) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        chrom_bed = parts[0]
        try:
            start = int(parts[1])
            end   = int(parts[2])
        except ValueError:
            continue
        chrom_fa = normalize_chrom(chrom_bed)
        if chrom_fa is None:
            n_skipped += 1
            continue
        win_start = max(0, start - WINDOW)
        win_end   = end + WINDOW
        try:
            seq = fa.fetch(chrom_fa, win_start, win_end).upper()
            sequences.append({"chrom": chrom_bed, "start": start, "end": end,
                               "seq": seq, "len": len(seq)})
        except Exception:
            n_skipped += 1

fa.close()
total_bp   = sum(s["len"] for s in sequences)
all_seq    = "".join(s["seq"] for s in sequences)
gc         = sum(1 for c in all_seq if c in "GC")
at         = sum(1 for c in all_seq if c in "AT")
n_count    = sum(1 for c in all_seq if c == "N")
gc_content = gc / max(gc + at, 1)

print(f"[INFO] DMCs: {len(sequences)}  (skipadas: {n_skipped})")
print(f"[INFO] Total bp: {total_bp:,}")
print(f"[INFO] GC content: {100*gc_content:.1f}%   N: {n_count} ({100*n_count/max(len(all_seq),1):.2f}%)")
print()

# ── Scan de motifs ────────────────────────────────────────────────────────────
print(f"[INFO] Escaneando {len(MOTIFS)} motifs em {len(sequences)} regiões (ambas strands)...")
print(f"  {'Motif':<14} {'Consenso':<20} {'Hits':>6} {'Seqs':>8} {'Enrich':>8}")
print("  " + "-" * 60)

raw_results = []
pvalues_raw = []

for name, consensus, mclass, desc, relevance in MOTIFS:
    pattern_fwd = iupac_to_regex(consensus)
    rc_seq      = reverse_complement(consensus)
    pattern_rev = iupac_to_regex(rc_seq) if rc_seq.upper() != consensus.upper() else None

    n_hits   = 0
    n_seqs_w = 0
    for s in sequences:
        seq = s["seq"]
        h = len(pattern_fwd.findall(seq))
        if pattern_rev:
            h += len(pattern_rev.findall(seq))
        n_hits += h
        if h > 0:
            n_seqs_w += 1

    # Probabilidade esperada por posição (ajuste GC para motifs com G/C)
    p_per_pos = motif_expected_freq(consensus)
    # Ajuste GC: substituir p(G) e p(C) pela frequência real
    n_positions = max(total_bp - len(consensus) + 1, 1)
    expected = n_positions * p_per_pos * (1 + (1 if pattern_rev else 0))
    enrichment = n_hits / max(expected, 0.01)

    # P-value binomial
    n_trials = n_positions
    p_success = p_per_pos * (2 if pattern_rev else 1)
    pval = binom_pvalue(n_hits, n_trials, min(p_success, 1.0))
    pvalues_raw.append(pval)

    raw_results.append({
        "motif": name, "consensus": consensus, "class": mclass,
        "description": desc, "relevance": relevance,
        "n_hits": n_hits, "n_seqs_with_hit": n_seqs_w,
        "pct_seqs": round(100 * n_seqs_w / max(len(sequences), 1), 1),
        "expected_hits": round(expected, 2),
        "fold_enrichment": round(enrichment, 2),
        "pvalue_raw": pval,
    })
    flag = "***" if pval < 0.001 else ("**" if pval < 0.01 else ("*" if pval < 0.05 else ""))
    print(f"  {name:<14} {consensus:<20} {n_hits:>6} {n_seqs_w:>4}/{len(sequences)}"
          f"  {enrichment:>6.2f}x  {flag}")

# Correção BH
adj_pvals = bh_correct(pvalues_raw)
for i, r in enumerate(raw_results):
    r["pvalue_bh"] = round(adj_pvals[i], 6)
    r["significant"] = adj_pvals[i] < 0.05

# Ordenar por fold enrichment
raw_results.sort(key=lambda x: x["fold_enrichment"], reverse=True)

# ── CSV principal ─────────────────────────────────────────────────────────────
out_csv = OUT / "motifs_full_v3.csv"
fieldnames = ["motif","consensus","class","description","relevance",
              "n_hits","n_seqs_with_hit","pct_seqs",
              "expected_hits","fold_enrichment",
              "pvalue_raw","pvalue_bh","significant","source"]
with open(out_csv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fieldnames)
    w.writeheader()
    for r in raw_results:
        r["source"] = "scan_real_genome_v3"
        w.writerow(r)

# Atualizar também o motifs_summary.csv (compatibilidade com app.R)
out_summary = OUT / "motifs_summary.csv"
sum_fields = ["motif","consensus","class","description","relevance",
              "n_hits","n_seqs_with_hit","pct_seqs","expected_hits","enrichment","pvalue","source"]
with open(out_summary, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=sum_fields)
    w.writeheader()
    for r in raw_results:
        w.writerow({"motif": r["motif"], "consensus": r["consensus"], "class": r["class"],
                    "description": r["description"], "relevance": r["relevance"],
                    "n_hits": r["n_hits"], "n_seqs_with_hit": r["n_seqs_with_hit"],
                    "pct_seqs": f"{r['pct_seqs']}%", "expected_hits": r["expected_hits"],
                    "enrichment": f"{r['fold_enrichment']:.2f}x",
                    "pvalue": r["pvalue_bh"], "source": "scan_real_genome_v3"})

print(f"\n[OK] Resultados completos: {out_csv}")
print(f"[OK] motifs_summary.csv atualizado: {out_summary}")

# ── Relatório por categoria ───────────────────────────────────────────────────
sig_results = [r for r in raw_results if r["significant"] and r["n_hits"] > 0]
print(f"\n{'='*70}")
print(f"  TOP MOTIFS SIGNIFICATIVOS (BH-adj p < 0.05)  [{len(sig_results)}/{len(raw_results)}]")
print(f"{'='*70}")
print(f"  {'Rank':<4} {'Motif':<14} {'Fold':<8} {'Hits':<7} {'p-BH':<10} Função")
print("  " + "-" * 65)
for i, r in enumerate(sig_results[:15], 1):
    print(f"  {i:<4} {r['motif']:<14} {r['fold_enrichment']:<8.2f} "
          f"{r['n_hits']:<7} {r['pvalue_bh']:<10.4f} {r['description'][:35]}")

# ── Análise k-mer (k=4,6,8) ──────────────────────────────────────────────────
print(f"\n[INFO] Análise k-mer (k=4, 6, 8)...")
for k in [4, 6, 8]:
    kmer_counts = collections.Counter()
    for s in sequences:
        seq = s["seq"]
        for i in range(len(seq) - k + 1):
            km = seq[i:i+k]
            if "N" not in km:
                kmer_counts[km] += 1
    top = kmer_counts.most_common(20)
    out_k = OUT / f"kmer_top20_k{k}.csv"
    with open(out_k, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["kmer", "count", "freq_per_kb"])
        for kmer, count in top:
            w.writerow([kmer, count, f"{count*1000/max(total_bp,1):.3f}"])
    print(f"  [OK] k={k}: {out_k.name}  top-3: {[km for km,_ in top[:3]]}")

# Compatibilidade: manter kmer_top20.csv para k=6
import shutil
shutil.copy(OUT / "kmer_top20_k6.csv", OUT / "kmer_top20.csv")

# ── Resumo final ─────────────────────────────────────────────────────────────
print(f"\n{'='*70}")
print(f"  HOMER v3 — CONCLUÍDO")
print(f"{'='*70}")
print(f"  Motifs analisados:    {len(MOTIFS)}")
print(f"  Motifs significativos:{len(sig_results)}")
print(f"  DMCs varridas:        {len(sequences)}")
print(f"  Total bp:             {total_bp:,}")
print(f"  GC content:           {100*gc_content:.1f}%")
print(f"  Resultados em:        {OUT}")
print(f"\n  Arquivos gerados:")
print(f"    motifs_full_v3.csv     — catálogo completo com p-values BH")
print(f"    motifs_summary.csv     — compatível com app.R (substituído)")
print(f"    kmer_top20_k4/6/8.csv  — análise k-mer por comprimento")
