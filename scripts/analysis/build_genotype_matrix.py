#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build Genotype View Matrix para MOFA2
======================================
Extrai frequências alélicas de SNPs nos genes-alvo para amostras Angus e Nelore
do 1000 Bull Genomes (TAUIND VCF), e monta uma matriz samples × SNPs para
integração como segunda view no MOFA2 (methylation + genotype).

Estratégia:
  - 20 amostras epigenoma: 10 Angus (A01-A05 × fev/jun) + 10 Nelore (N01-N05 × fev/jun)
  - VCF: 112 Angus + 385 Nelore com IDs SAMN
  - Para cada gene-alvo: extrair genótipos (0/0, 0/1, 1/1) de Angus e Nelore
  - Calcular frequência alélica alternativa por raça: freq_alt = (n_alt_alelos) / (2 × n_amostras)
  - Cada amostra epigenoma recebe o perfil de frequência da sua raça
  - Resultado: matrix 20 × N_SNPs (N_SNPs = SNPs filtrados: MAF > 1%, genotyped > 80%)

Output: genotype_matrix.csv (samples × SNPs) + genotype_snp_info.csv (SNP metadata)
"""
from __future__ import annotations
import csv
import datetime
import gzip
import json
import re
import sys
import zlib
from collections import defaultdict
from pathlib import Path

# pysam path
for _pip_path in [
    '/usr/local/lib/python3.10/dist-packages',
    '/usr/lib/python3/dist-packages',
    os.path.expanduser('~/.local/lib/python3.10/site-packages'),
]:
    if _pip_path not in sys.path:
        sys.path.insert(0, _pip_path)

# ─── Caminhos ─────────────────────────────────────────────────────────────────
VCF_DIR    = Path("/data/1000bg_vcf")
META_DIR   = Path("/data/1000bg_metadata")
JSON_PATH  = Path("./subprojeto/SUB001/SUB001_target_genes_coordinates.json")
SAMPLES_CSV = Path("./samples_GSE161113_demo.csv")
OUT_DIR    = Path("./scripts/results/SUB001_GSE161113_run_oficial/fase6_mofa2")

# Filtros de qualidade
MIN_MAF          = 0.01   # frequência alélica mínima (em qualquer raça)
MIN_GENO_RATE    = 0.80   # fração mínima de amostras genotipadas por raça
WINDOW           = 500    # flancos adicionais ao redor do gene (bp)
MAX_SNPS_PER_GENE = 50    # máximo SNPs por gene (seleciona por MAF decrescente)

# ─── Carregar metadados ────────────────────────────────────────────────────────
def load_samples():
    """Retorna dict sample_id -> breed (angus|nelore)"""
    angus_ids  = set((META_DIR / "samples_angus.txt").read_text().split())
    nelore_ids = set((META_DIR / "samples_nelore.txt").read_text().split())
    return angus_ids, nelore_ids

def load_json_genes():
    with open(JSON_PATH) as f:
        data = json.load(f)
    genes = data.get("genes", {})
    if isinstance(genes, list):
        genes = {g.get("gene_name", f"gene_{i}"): g for i, g in enumerate(genes)}
    return genes

def load_epigenoma_samples():
    """Retorna lista de (sample_id, breed) do arquivo CSV do epigenoma."""
    samples = []
    with open(SAMPLES_CSV) as f:
        for row in csv.DictReader(f):
            breed = row.get("breed", "").lower()
            if "angus" in breed:
                b = "angus"
            elif "nelore" in breed or "nelo" in breed:
                b = "nelore"
            else:
                b = breed
            samples.append((row["sampleID"], b))
    return samples

# ─── VCF helpers ──────────────────────────────────────────────────────────────
def find_vcf(chrom_norm: str):
    patterns = [
        re.compile(rf'[Cc]hr{re.escape(chrom_norm)}[-_\.]', re.IGNORECASE),
    ]
    for vcf in sorted(VCF_DIR.glob("*.vcf.gz")):
        for pat in patterns:
            if pat.search(vcf.name):
                return vcf
    return None

def norm_chrom(c):
    c = str(c)
    for pfx in ["Chr", "chr", "CHR"]:
        if c.startswith(pfx):
            c = c[len(pfx):]
    return c

def parse_gt(gt_str: str):
    """Parse GT field: returns (a1, a2) as ints or None if missing."""
    if not gt_str or gt_str in (".", "./.", ".|."):
        return None
    sep = "|" if "|" in gt_str else "/"
    parts = gt_str.split(sep)[:2]
    try:
        return int(parts[0]), int(parts[1])
    except (ValueError, IndexError):
        return None

def extract_genotypes_tabix(vcf_path: Path, regions: list[dict],
                             angus_cols: list[int], nelore_cols: list[int],
                             chrom_vcf: str) -> list[dict]:
    """
    Extrai genótipos por SNP para amostras Angus e Nelore via tabix.
    Retorna lista de SNP records com frequências alélicas.
    """
    import pysam

    snps = []
    tbx = pysam.TabixFile(str(vcf_path))

    for reg in regions:
        gene      = reg["gene"]
        start_0   = max(0, reg["start"] - 1)  # tabix 0-based
        end       = reg["end"]

        try:
            for rec in tbx.fetch(chrom_vcf, start_0, end):
                fields = rec.split("\t")
                if len(fields) < 10:
                    continue
                pos    = int(fields[1])
                ref    = fields[3]
                alt    = fields[4]
                fmt    = fields[8].split(":")

                # Apenas SNPs bialélicos simples
                if len(ref) != 1 or len(alt) != 1 or alt == ".":
                    continue
                if ref not in "ACGT" or alt not in "ACGT":
                    continue

                # Índice do campo GT no FORMAT
                try:
                    gt_idx = fmt.index("GT")
                except ValueError:
                    continue

                sample_fields = fields[9:]

                angus_ref2 = angus_alt2 = angus_miss = 0
                for col in angus_cols:
                    if col >= len(sample_fields):
                        angus_miss += 1
                        continue
                    gt_field = sample_fields[col].split(":")[gt_idx]
                    gt = parse_gt(gt_field)
                    if gt is None:
                        angus_miss += 1
                        continue
                    a1, a2 = gt
                    angus_ref2 += (a1 == 0) + (a2 == 0)
                    angus_alt2 += (a1 == 1) + (a2 == 1)

                nelore_ref2 = nelore_alt2 = nelore_miss = 0
                for col in nelore_cols:
                    if col >= len(sample_fields):
                        nelore_miss += 1
                        continue
                    gt_field = sample_fields[col].split(":")[gt_idx]
                    gt = parse_gt(gt_field)
                    if gt is None:
                        nelore_miss += 1
                        continue
                    a1, a2 = gt
                    nelore_ref2 += (a1 == 0) + (a2 == 0)
                    nelore_alt2 += (a1 == 1) + (a2 == 1)

                n_angus   = len(angus_cols)  - angus_miss
                n_nelore  = len(nelore_cols) - nelore_miss
                geno_rate_angus  = n_angus  / max(len(angus_cols),  1)
                geno_rate_nelore = n_nelore / max(len(nelore_cols), 1)

                if geno_rate_angus  < MIN_GENO_RATE: continue
                if geno_rate_nelore < MIN_GENO_RATE: continue

                freq_angus  = angus_alt2  / max(angus_ref2  + angus_alt2,  1)
                freq_nelore = nelore_alt2 / max(nelore_ref2 + nelore_alt2, 1)

                # MAF = min(freq, 1-freq) — filtro global
                maf = max(
                    min(freq_angus,  1 - freq_angus),
                    min(freq_nelore, 1 - freq_nelore)
                )
                if maf < MIN_MAF:
                    continue

                snp_id = f"{gene}_{chrom_vcf}_{pos}_{ref}_{alt}"
                snps.append({
                    "snp_id":        snp_id,
                    "gene":          gene,
                    "chrom":         chrom_vcf,
                    "pos":           pos,
                    "ref":           ref,
                    "alt":           alt,
                    "freq_angus":    round(freq_angus,  5),
                    "freq_nelore":   round(freq_nelore, 5),
                    "maf":           round(maf, 5),
                    "n_angus":       n_angus,
                    "n_nelore":      n_nelore,
                })

        except Exception as e:
            print(f"  [WARN] tabix fetch error {gene}: {e}", file=sys.stderr)

    tbx.close()
    return snps


def extract_genotypes_linear(vcf_path: Path, regions: list[dict],
                              angus_cols: list[int], nelore_cols: list[int],
                              chrom_vcf: str) -> list[dict]:
    """Fallback linear scan."""
    intervals = [(r["start"]-WINDOW, r["end"]+WINDOW, r["gene"]) for r in regions]
    gene_by_pos = {(s, e): g for s, e, g in intervals}
    min_pos = min(s for s, e, _ in intervals)
    max_pos = max(e for _, e, _ in intervals)

    snps = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    try:
        with opener(str(vcf_path), "rt", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 10:
                    continue
                try:
                    pos = int(fields[1])
                except ValueError:
                    continue
                if pos < min_pos or pos > max_pos:
                    continue

                ref = fields[3]
                alt = fields[4]
                if len(ref) != 1 or len(alt) != 1 or alt == ".":
                    continue

                # find gene
                gene = None
                for s, e, g in intervals:
                    if s <= pos <= e:
                        gene = g
                        break
                if gene is None:
                    continue

                fmt = fields[8].split(":")
                try:
                    gt_idx = fmt.index("GT")
                except ValueError:
                    continue

                sample_fields = fields[9:]
                angus_ref2 = angus_alt2 = angus_miss = 0
                for col in angus_cols:
                    if col >= len(sample_fields):
                        angus_miss += 1
                        continue
                    gt_field = sample_fields[col].split(":")[gt_idx]
                    gt = parse_gt(gt_field)
                    if gt is None:
                        angus_miss += 1
                        continue
                    a1, a2 = gt
                    angus_ref2 += (a1 == 0) + (a2 == 0)
                    angus_alt2 += (a1 == 1) + (a2 == 1)

                nelore_ref2 = nelore_alt2 = nelore_miss = 0
                for col in nelore_cols:
                    if col >= len(sample_fields):
                        nelore_miss += 1
                        continue
                    gt_field = sample_fields[col].split(":")[gt_idx]
                    gt = parse_gt(gt_field)
                    if gt is None:
                        nelore_miss += 1
                        continue
                    a1, a2 = gt
                    nelore_ref2 += (a1 == 0) + (a2 == 0)
                    nelore_alt2 += (a1 == 1) + (a2 == 1)

                n_angus  = len(angus_cols)  - angus_miss
                n_nelore = len(nelore_cols) - nelore_miss
                if n_angus  / max(len(angus_cols),  1) < MIN_GENO_RATE: continue
                if n_nelore / max(len(nelore_cols), 1) < MIN_GENO_RATE: continue

                freq_angus  = angus_alt2  / max(angus_ref2  + angus_alt2,  1)
                freq_nelore = nelore_alt2 / max(nelore_ref2 + nelore_alt2, 1)
                maf = max(
                    min(freq_angus,  1 - freq_angus),
                    min(freq_nelore, 1 - freq_nelore)
                )
                if maf < MIN_MAF:
                    continue

                snp_id = f"{gene}_{chrom_vcf}_{pos}_{ref}_{alt}"
                snps.append({
                    "snp_id":      snp_id,
                    "gene":        gene,
                    "chrom":       chrom_vcf,
                    "pos":         pos,
                    "ref":         ref,
                    "alt":         alt,
                    "freq_angus":  round(freq_angus,  5),
                    "freq_nelore": round(freq_nelore, 5),
                    "maf":         round(maf, 5),
                    "n_angus":     n_angus,
                    "n_nelore":    n_nelore,
                })
    except (OSError, EOFError, zlib.error, Exception) as e:
        print(f"  [WARN] linear scan error: {e}", file=sys.stderr)
    return snps


# ─── Main ─────────────────────────────────────────────────────────────────────
def main():
    t_start = datetime.datetime.now()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[{t_start:%H:%M:%S}] Build Genotype Matrix para MOFA2")
    print(f"  MIN_MAF={MIN_MAF}, MIN_GENO_RATE={MIN_GENO_RATE}, WINDOW={WINDOW}")

    # Carregar dados
    angus_ids, nelore_ids = load_samples()
    genes = load_json_genes()
    epi_samples = load_epigenoma_samples()

    print(f"\nAmostras VCF: {len(angus_ids)} Angus, {len(nelore_ids)} Nelore")
    print(f"Amostras epigenoma: {len(epi_samples)}")
    print(f"Genes-alvo: {len(genes)}")

    # Verificar pysam
    has_pysam = False
    try:
        import pysam
        has_pysam = True
        print("[INFO] pysam disponível → tabix ativado")
    except ImportError:
        print("[INFO] pysam não disponível → fallback linear")

    # Pegar índices de coluna para Angus/Nelore em Chr1 VCF
    chr1_vcf = next((v for v in sorted(VCF_DIR.glob("*.vcf.gz"))
                     if re.search(r'Chr1[-_.]', v.name)), None)
    all_vcf_samples = []
    with gzip.open(str(chr1_vcf), 'rt', errors='ignore') as f:
        for line in f:
            if line.startswith("#CHROM"):
                all_vcf_samples = line.strip().split("\t")[9:]
                break
            if not line.startswith("#"):
                break

    vcf_sample_idx = {s: i for i, s in enumerate(all_vcf_samples)}
    angus_cols  = [vcf_sample_idx[s] for s in angus_ids  if s in vcf_sample_idx]
    nelore_cols = [vcf_sample_idx[s] for s in nelore_ids if s in vcf_sample_idx]
    print(f"\nÍndices VCF: {len(angus_cols)} Angus, {len(nelore_cols)} Nelore")

    # Agrupar genes por cromossomo
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for gene_name, info in genes.items():
        chrom = norm_chrom(info.get("chromosome", ""))
        start = info.get("start", 0) - WINDOW
        end   = info.get("end",   0) + WINDOW
        by_chrom[chrom].append({"gene": gene_name, "start": start, "end": end})

    # Extrair SNPs por cromossomo
    all_snps: list[dict] = []
    print()
    for chrom_norm, regions in sorted(by_chrom.items(), key=lambda x: x[0].zfill(3)):
        vcf_path = find_vcf(chrom_norm)
        if vcf_path is None:
            print(f"  [skip] Chr{chrom_norm}: sem VCF", file=sys.stderr)
            continue

        t0 = datetime.datetime.now()
        print(f"[{t0:%H:%M:%S}] Chr{chrom_norm} ({len(regions)} genes) | {vcf_path.name}", flush=True)

        # Detectar formato do cromossomo no VCF
        chrom_vcf = chrom_norm
        if has_pysam:
            import pysam
            try:
                tbx = pysam.TabixFile(str(vcf_path))
                contigs = list(tbx.contigs)
                tbx.close()
                for cand in [chrom_norm, f"Chr{chrom_norm}", f"chr{chrom_norm}"]:
                    if cand in contigs:
                        chrom_vcf = cand
                        break
            except Exception:
                pass

        # Extrair genótipos
        if has_pysam:
            try:
                import pysam
                tbi = Path(str(vcf_path) + ".tbi")
                if tbi.exists():
                    snps = extract_genotypes_tabix(vcf_path, regions, angus_cols, nelore_cols, chrom_vcf)
                    print(f"  → [tabix] {len(snps)} SNPs (MAF≥{MIN_MAF:.0%})", flush=True)
                else:
                    snps = extract_genotypes_linear(vcf_path, regions, angus_cols, nelore_cols, chrom_vcf)
                    print(f"  → [linear] {len(snps)} SNPs", flush=True)
            except Exception as e:
                print(f"  [WARN] {e}", file=sys.stderr)
                snps = extract_genotypes_linear(vcf_path, regions, angus_cols, nelore_cols, chrom_vcf)
        else:
            snps = extract_genotypes_linear(vcf_path, regions, angus_cols, nelore_cols, chrom_vcf)

        # Limitar SNPs por gene (top MAF)
        gene_snps: dict[str, list] = defaultdict(list)
        for snp in snps:
            gene_snps[snp["gene"]].append(snp)

        kept = []
        for g, gs in gene_snps.items():
            gs_sorted = sorted(gs, key=lambda x: -x["maf"])
            kept.extend(gs_sorted[:MAX_SNPS_PER_GENE])

        print(f"  → {len(kept)} SNPs após filtro ({MAX_SNPS_PER_GENE}/gene)")
        all_snps.extend(kept)

    print(f"\nTotal SNPs extraídos: {len(all_snps)}")

    if not all_snps:
        print("[ERROR] Nenhum SNP extraído. Verifique VCFs e filtros.")
        sys.exit(1)

    # ─── Montar matriz samples × SNPs ─────────────────────────────────────────
    # SNP IDs (colunas)
    snp_ids = [s["snp_id"] for s in all_snps]
    # Frequências por raça
    freq_map = {
        s["snp_id"]: {"angus": s["freq_angus"], "nelore": s["freq_nelore"]}
        for s in all_snps
    }

    print(f"\n=== Montando matriz {len(epi_samples)} × {len(snp_ids)} ===")
    geno_matrix = []
    for sample_id, breed in epi_samples:
        row = {"sample": sample_id}
        for snp_id in snp_ids:
            row[snp_id] = freq_map[snp_id].get(breed, "NA")
        geno_matrix.append(row)

    # ─── Salvar outputs ────────────────────────────────────────────────────────
    out_matrix = OUT_DIR / "genotype_matrix.csv"
    out_info   = OUT_DIR / "genotype_snp_info.csv"

    with open(out_matrix, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample"] + snp_ids)
        writer.writeheader()
        writer.writerows(geno_matrix)
    print(f"✓ Matrix salva: {out_matrix}")

    snp_info_fields = ["snp_id", "gene", "chrom", "pos", "ref", "alt",
                       "freq_angus", "freq_nelore", "maf", "n_angus", "n_nelore"]
    with open(out_info, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=snp_info_fields)
        writer.writeheader()
        writer.writerows(all_snps)
    print(f"✓ SNP info salva: {out_info}")

    # Resumo
    elapsed = (datetime.datetime.now() - t_start).total_seconds()
    print(f"\n=== RESUMO ===")
    print(f"  Samples: {len(geno_matrix)}")
    print(f"  SNPs: {len(snp_ids)}")
    print(f"  Tempo: {elapsed:.1f}s")

    # Genes com SNPs
    from collections import Counter
    gene_counts = Counter(s["gene"] for s in all_snps)
    print(f"\n  SNPs por gene (top 10):")
    for g, n in gene_counts.most_common(10):
        print(f"    {g}: {n}")

    print(f"\n  Frequências alélicas amostra (primeiras 3 SNPs):")
    for row in geno_matrix[:3]:
        items = list(row.items())
        print(f"    {items[0][1]}: {items[1]}, {items[2]}, {items[3]}")


if __name__ == "__main__":
    main()
