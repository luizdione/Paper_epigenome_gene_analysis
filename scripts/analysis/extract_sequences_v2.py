#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_sequences_v2.py — Extração de sequências FASTA dos genes-alvo
usando GENOMA DE REFERÊNCIA + VARIANTES do 1000 Bull Genomes.

VANTAGENS sobre v1 (extract_sequences_from_vcf.py):
  • FASTA com 100% bases conhecidas (não mais 95% N)
  • Robusto a corrupção parcial dos VCFs (genoma referência preenche o resto)
  • Permite gerar consensos personalizados por raça/amostra
  • Suporta análises downstream que precisam de sequência completa
    (HOMER motif scanning, alinhamento, design de primers)

ENTRADAS:
  • --reference: genoma ARS-UCD1.3 completo (FASTA indexado via samtools faidx)
  • --vcf-dir: diretório com VCFs do 1000 Bull Genomes (opcional, para variantes)
  • --json: arquivo de coordenadas dos genes-alvo do subprojeto

SAÍDA:
  • <SUB>_reference_sequences.fa — uma sequência por gene, no formato:
        >GENE_chrN:start-end_REF
        ACGTACGT... (sequência REF do genoma)

  • Opcional: <SUB>_consensus_<RACA>.fa — consenso por raça aplicando variantes
        >GENE_chrN:start-end_Angus
        ACGTNCGT... (REF com substituições conforme genótipos)

USO:
  python3 extract_sequences_v2.py \\
      --reference /data/genome_reference/ARS-UCD1.3/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa \\
      --vcf-dir /data/1000bg_vcf \\
      --json /home/rstudio/projeto/subprojeto/SUB01/SUB01_target_genes_coordinates.json \\
      --out /home/rstudio/projeto/subprojeto/SUB01/SUB01_reference_sequences.fa \\
      --window 500
"""
from __future__ import annotations
import argparse
import datetime
import json
import re
import subprocess
import sys
from pathlib import Path

# ============================================================
# CLI
# ============================================================
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--reference", required=True,
                help="Genoma de referência ARS-UCD1.3 (FASTA, indexado .fai)")
ap.add_argument("--vcf-dir", default=None,
                help="Diretório com VCFs (opcional, para variantes/consenso)")
ap.add_argument("--json", required=True,
                help="JSON com coordenadas dos genes-alvo")
ap.add_argument("--out", required=True,
                help="Arquivo FASTA de saída")
ap.add_argument("--window", type=int, default=500,
                help="Flanking ±N pb ao redor de cada gene (default: 500)")
ap.add_argument("--integrate-variants", action="store_true",
                help="Integrar variantes do 1000bg ao FASTA (REF/ALT por raça)")
ap.add_argument("--breeds", default=None,
                help="Raças para consenso separado (ex: Angus,Nelore)")
args = ap.parse_args()

REF_FA = Path(args.reference)
VCF_DIR = Path(args.vcf_dir) if args.vcf_dir else None
JSON_PATH = Path(args.json)
OUT_FA = Path(args.out)
WINDOW = args.window

# ============================================================
# VALIDAÇÃO INICIAL
# ============================================================
if not REF_FA.exists():
    print(f"[ERRO] Genoma de referência não encontrado: {REF_FA}")
    sys.exit(1)

REF_FAI = Path(str(REF_FA) + ".fai")
if not REF_FAI.exists():
    print(f"[ERRO] Índice .fai não encontrado. Rode: samtools faidx {REF_FA}")
    sys.exit(1)

if not JSON_PATH.exists():
    print(f"[ERRO] JSON de genes não encontrado: {JSON_PATH}")
    sys.exit(1)

OUT_FA.parent.mkdir(parents=True, exist_ok=True)

# Tentar importar pysam para extração rápida
try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False

# ============================================================
# HELPERS
# ============================================================
def norm_chrom(c):
    """Remove prefixo 'chr'/'Chr' do cromossomo."""
    c = str(c)
    return re.sub(r"^[Cc][Hh][Rr]", "", c)


def find_chrom_in_fai(fai_path, chrom_norm):
    """Descobre o nome real do cromossomo no FASTA (com ou sem prefixo)."""
    with open(fai_path) as f:
        contigs = [line.split("\t")[0] for line in f]
    candidates = [chrom_norm, "chr" + chrom_norm, "Chr" + chrom_norm]
    for c in candidates:
        if c in contigs:
            return c
    return None


def extract_via_samtools(ref_fa, chrom, start, end):
    """Usa samtools faidx para extração rápida."""
    if subprocess_which("samtools") is None:
        return None
    region = f"{chrom}:{start}-{end}"
    r = subprocess.run(["samtools", "faidx", str(ref_fa), region],
                       capture_output=True, timeout=60)
    if r.returncode != 0:
        return None
    out = r.stdout.decode("ascii", errors="replace")
    lines = out.strip().split("\n")
    if len(lines) < 2:
        return None
    seq = "".join(lines[1:]).upper().replace("\n", "")
    return seq


def extract_via_pysam(ref_fa, chrom, start, end):
    """Usa pysam para extração rápida."""
    try:
        f = pysam.FastaFile(str(ref_fa))
        try:
            seq = f.fetch(chrom, start - 1, end)
            return seq.upper()
        finally:
            f.close()
    except Exception as e:
        print(f"  [WARN] pysam fetch {chrom}:{start}-{end}: {e}")
        return None


def extract_manual(ref_fa, fai_path, chrom, start, end):
    """Fallback: usa offsets do .fai para fetch manual."""
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] == chrom:
                seq_len = int(parts[1])
                offset = int(parts[2])
                line_blen = int(parts[3])
                line_len = int(parts[4])
                # Calcular byte offsets
                start_byte_in_seq = (start - 1) // line_blen * line_len + (start - 1) % line_blen
                end_byte_in_seq = (end - 1) // line_blen * line_len + (end - 1) % line_blen + 1
                with open(ref_fa, "rb") as fa:
                    fa.seek(offset + start_byte_in_seq)
                    raw = fa.read(end_byte_in_seq - start_byte_in_seq)
                seq = raw.decode("ascii", errors="replace").replace("\n", "").replace("\r", "")
                return seq.upper()[:end - start + 1]
    return None


def subprocess_which(cmd):
    """Detecta se 'cmd' está no PATH."""
    import shutil as _sh
    return _sh.which(cmd)


# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 70)
    print(f"  EXTRACT SEQUENCES v2 — Genoma de Referência + Variantes")
    print(f"  {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    print(f"  Reference: {REF_FA}")
    print(f"  Index:     {REF_FAI}")
    print(f"  JSON:      {JSON_PATH}")
    print(f"  Out:       {OUT_FA}")
    print(f"  Window:    ±{WINDOW} bp")
    print(f"  pysam:     {'SIM' if HAS_PYSAM else 'NAO (fallback manual)'}")
    print(f"  samtools:  {'SIM' if subprocess_which('samtools') else 'NAO'}")
    print()

    # Carregar genes-alvo
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    genes = data.get("genes", {})
    if isinstance(genes, list):
        genes = {g.get("gene_name", f"gene_{i}"): g for i, g in enumerate(genes)}
    print(f"[INFO] Genes-alvo carregados: {len(genes)}")

    # Detectar prefixo do FAI
    with open(REF_FAI) as f:
        first_contig = f.readline().split("\t")[0]
    print(f"[INFO] Primeiro contig no FAI: {first_contig}")
    fai_has_chr_prefix = first_contig.lower().startswith("chr")
    print(f"[INFO] FAI usa prefixo 'chr': {fai_has_chr_prefix}")
    print()

    # Extração
    extracted = 0
    skipped = []
    use_method = "pysam" if HAS_PYSAM else ("samtools" if subprocess_which("samtools") else "manual")
    print(f"[INFO] Método de extração: {use_method}")
    print()

    with open(OUT_FA, "w", encoding="ascii") as out:
        for gene_name, info in genes.items():
            chrom = norm_chrom(info.get("chromosome", ""))
            start = int(info.get("start", 0)) - WINDOW
            end = int(info.get("end", 0)) + WINDOW
            if start < 1:
                start = 1
            label = f"{gene_name}_chr{chrom}:{start}-{end}"

            # Encontrar nome do contig no FAI
            chrom_in_fai = find_chrom_in_fai(REF_FAI, chrom)
            if chrom_in_fai is None:
                print(f"  [SKIP] {gene_name}: chr{chrom} não encontrado no genoma")
                skipped.append(gene_name)
                continue

            # Extrair via método disponível
            seq = None
            if use_method == "pysam":
                seq = extract_via_pysam(REF_FA, chrom_in_fai, start, end)
            elif use_method == "samtools":
                seq = extract_via_samtools(REF_FA, chrom_in_fai, start, end)
            else:
                seq = extract_manual(REF_FA, REF_FAI, chrom_in_fai, start, end)

            if not seq or len(seq) < 10:
                print(f"  [SKIP] {gene_name}: extração falhou ({len(seq) if seq else 0} bp)")
                skipped.append(gene_name)
                continue

            # Stats
            n_known = sum(1 for c in seq if c in "ACGT")
            density = 100 * n_known / len(seq)

            # Escrever
            out.write(f">{label}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")
            extracted += 1
            if extracted <= 5 or extracted % 10 == 0:
                print(f"  [{extracted:2d}] {gene_name}: chr{chrom}:{start:,}-{end:,} "
                      f"({len(seq):,} bp, {density:.1f}% known)")

    print()
    print("=" * 70)
    print(f"  CONCLUÍDO")
    print("=" * 70)
    print(f"  Extraídos: {extracted}/{len(genes)}")
    if skipped:
        print(f"  Skipados:  {len(skipped)} → {skipped}")
    print(f"  FASTA:     {OUT_FA} ({OUT_FA.stat().st_size:,} bytes)")
    print()

    # Tentar gerar índice
    if subprocess_which("samtools"):
        print(f"[INFO] Gerando índice samtools faidx...")
        r = subprocess.run(["samtools", "faidx", str(OUT_FA)], capture_output=True)
        if r.returncode == 0:
            print(f"[OK] Índice gerado: {OUT_FA}.fai")
    elif HAS_PYSAM:
        print(f"[INFO] Gerando índice via pysam.faidx...")
        try:
            pysam.faidx(str(OUT_FA))
            print(f"[OK] Índice gerado: {OUT_FA}.fai")
        except Exception as e:
            print(f"[WARN] pysam.faidx falhou: {e}")

    if args.integrate_variants and VCF_DIR:
        print()
        print("[NOTA] Integração de variantes (--integrate-variants) ainda não implementada.")
        print("       Para consenso por raça, será necessário rodar bcftools consensus")
        print("       ou implementar lógica adicional. Por ora o FASTA tem REF do genoma.")


if __name__ == "__main__":
    main()
