#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extração de Sequências de Referência via VCFs do 1000 Bull Genomes
===================================================================
Integração entre os projetos Epigenoma (DMCs/genes-alvo) e Genomas (1000bg_vcf).

DESIGN:
  Este script foi projetado para rodar UMA VEZ (ex: durante a madrugada) e
  gerar um FASTA de saída reutilizável por todas as análises do subprojeto.
  O pipeline (run_pipeline.py) apenas VERIFICA se o FASTA existe — se sim,
  usa diretamente; se não, roda em modo fallback sem sequência.

ESTRATÉGIA DE LEITURA (eficiente via tabix):
  Usa pysam.TabixFile.fetch() para query regional direta nos VCFs indexados.
  Requer índice .tbi (gerado automaticamente via tabix se ausente).
  Fallback para leitura linear se pysam não disponível.

DOIS MODOS DE OPERAÇÃO:
  1. --mode dmcs  (padrão): extrai janelas ao redor dos DMCs da Fase 1
  2. --mode genes: extrai regiões completas dos genes-alvo do JSON (flancos ±window)

Uso:
  python3 scripts/extract_sequences_from_vcf.py \
      --mode genes \
      --json subprojeto/SUB001/SUB001_target_genes_coordinates.json \
      --vcf-dir /data/1000bg_vcf \
      --out subprojeto/SUB001/SUB001_reference_sequences.fa \
      --window 500

Autor: Luiz Dione Barbosa De-Melo — GenOmics Epigenoma 2026
"""
from __future__ import annotations
import argparse
import csv
import datetime
import gzip
import json
import zlib
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# Garantir que pacotes pip (pysam) sejam encontrados em diferentes instalações
for _pip_path in [
    '/usr/local/lib/python3.10/dist-packages',
    '/usr/lib/python3/dist-packages',
    os.path.expanduser('~/.local/lib/python3.10/site-packages'),
]:
    if _pip_path not in sys.path:
        sys.path.insert(0, _pip_path)


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def norm_chrom(c) -> str:
    """Normaliza cromossomo: remove 'chr' e converte para string."""
    return str(c).strip().lower().removeprefix("chr")


def vcf_chrom_number(vcf_path: str) -> str | None:
    """Extrai número/letra do cromossomo do nome do arquivo VCF."""
    m = re.search(r'[Cc]hr(\w+)', os.path.basename(vcf_path))
    return m.group(1).lower() if m else None


def find_vcf_for_chrom(vcf_dir: Path, chrom_norm: str) -> Path | None:
    """Localiza o arquivo VCF para um cromossomo."""
    for f in sorted(vcf_dir.glob("*.vcf.gz")):
        cn = vcf_chrom_number(str(f))
        if cn and cn.lower() == chrom_norm.lower():
            return f
    return None


def reconstruct_sequence(positions: dict[int, str], start: int, end: int) -> str:
    """
    Reconstrói sequência de referência para [start, end] (1-based, inclusivo).
    Posições sem SNP no VCF recebem 'N'.
    """
    length = end - start + 1
    seq = bytearray(b"N" * length)
    for pos, base in positions.items():
        idx = pos - start
        if 0 <= idx < length:
            seq[idx] = ord(base)
    return seq.decode("ascii")


def detect_chrom_prefix(vcf_path: Path) -> str:
    """Detecta se o VCF usa prefixo 'chr' nas colunas CHROM."""
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    try:
        with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                chrom = line.split("\t", 1)[0].strip()
                if chrom.lower().startswith("chr"):
                    return "chr"
                return ""
    except Exception:
        return ""
    return ""


def ensure_tabix_index(vcf_path: Path) -> bool:
    """Garante que o índice .tbi existe para o VCF. Retorna True se OK."""
    tbi_path = Path(str(vcf_path) + ".tbi")
    if tbi_path.exists():
        return True
    print(f"  [tabix] Criando índice para {vcf_path.name}...", flush=True)
    try:
        result = subprocess.run(
            ["tabix", "-p", "vcf", str(vcf_path)],
            capture_output=True, text=True, timeout=300
        )
        if result.returncode == 0:
            print(f"  [tabix] Índice criado.", flush=True)
            return True
        else:
            print(f"  [WARN] tabix falhou: {result.stderr.strip()}", file=sys.stderr)
            return False
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        print(f"  [WARN] tabix não disponível: {e}", file=sys.stderr)
        return False


# ═══════════════════════════════════════════════════════════════════════════
# Coleta de regiões de interesse
# ═══════════════════════════════════════════════════════════════════════════

def collect_regions_from_dmcs(dmcs_path: str, window: int) -> dict[str, list[dict]]:
    """Lê CSV de DMCs e retorna dict {chrom_norm -> [regiões]}."""
    regions_by_chrom: dict[str, list[dict]] = defaultdict(list)

    with open(dmcs_path, encoding="utf-8", errors="ignore", newline="") as fh:
        sample = fh.read(1024); fh.seek(0)
        delim = "\t" if sample.count("\t") > sample.count(",") else ","
        reader = csv.DictReader(fh, delimiter=delim)
        for i, row in enumerate(reader, 1):
            chrom_raw = (row.get("chr") or row.get("chrom") or "").strip()
            start_raw = row.get("start") or row.get("pos") or "0"
            end_raw   = row.get("end") or start_raw
            if not chrom_raw:
                continue
            try:
                center = (int(float(start_raw)) + int(float(end_raw))) // 2
            except ValueError:
                continue

            chrom_norm = norm_chrom(chrom_raw)
            win_start  = max(1, center - window)
            win_end    = center + window

            diff  = row.get("meth.diff", "")
            qval  = row.get("qvalue", "")
            label = f"DMC_{i}_{chrom_raw}:{win_start}-{win_end}"
            if diff:
                try: label += f"_diff={float(diff):.1f}%"
                except ValueError: pass
            if qval:
                try: label += f"_q={float(qval):.4f}"
                except ValueError: pass

            regions_by_chrom[chrom_norm].append({
                "chrom_raw": chrom_raw,
                "start": win_start,
                "end":   win_end,
                "label": label,
            })

    return regions_by_chrom


def collect_regions_from_json(json_path: str, window: int) -> dict[str, list[dict]]:
    """Lê JSON de genes-alvo e retorna dict {chrom_norm -> [regiões]}."""
    regions_by_chrom: dict[str, list[dict]] = defaultdict(list)

    with open(json_path, encoding="utf-8") as fh:
        data = json.load(fh)

    genes_raw = data.get("genes", {})
    # Normalizar lista -> dict
    if isinstance(genes_raw, list):
        genes = {}
        for i, item in enumerate(genes_raw):
            if isinstance(item, dict):
                name = (item.get("gene_name") or item.get("name") or
                        item.get("symbol") or f"gene_{i}")
                genes[name] = item
            else:
                name = str(item)
                genes[name] = {"gene_name": name}
    else:
        genes = genes_raw

    skipped = 0
    for gene_name, info in genes.items():
        chrom_raw = str(info.get("chromosome", "")).strip()
        if not chrom_raw:
            skipped += 1
            continue
        chrom_norm = norm_chrom(chrom_raw)

        gene_start = info.get("start")
        gene_end   = info.get("end")
        if gene_start is None or gene_end is None:
            skipped += 1
            continue

        reg_start = info.get("flanking_start") or max(1, int(gene_start) - window)
        reg_end   = info.get("flanking_end")   or (int(gene_end) + window)

        regions_by_chrom[chrom_norm].append({
            "chrom_raw": chrom_raw,
            "start": int(reg_start),
            "end":   int(reg_end),
            "label": f"{gene_name}_chr{chrom_raw}:{reg_start}-{reg_end}",
            "gene":  gene_name,
            "module": info.get("module", ""),
        })

    if skipped:
        print(f"  [INFO] {skipped} genes ignorados (sem cromossomo ou coordenadas).")

    return regions_by_chrom


# ═══════════════════════════════════════════════════════════════════════════
# Extração via tabix (rápida) — query regional direta
# ═══════════════════════════════════════════════════════════════════════════

def extract_positions_tabix(
    vcf_path: Path,
    regions: list[dict],
    chrom_prefix: str = ""
) -> dict[int, str]:
    """
    Usa pysam.TabixFile para query regional direta — muito mais rápido.
    Retorna dict {pos -> base_ref}.
    """
    import pysam  # type: ignore

    positions: dict[int, str] = {}
    tbx = pysam.TabixFile(str(vcf_path))

    vcf_chroms = set(tbx.contigs)

    for reg in regions:
        start = reg["start"]
        end   = reg["end"]
        chrom_norm = norm_chrom(reg["chrom_raw"])

        # Tentar vários formatos de cromossomo
        candidates = [
            chrom_prefix + chrom_norm,
            chrom_prefix + chrom_norm.upper(),
            "Chr" + chrom_norm,
            "CHR" + chrom_norm,
            chrom_norm,
        ]
        chrom_vcf = None
        for c in candidates:
            if c in vcf_chroms:
                chrom_vcf = c
                break

        if chrom_vcf is None:
            # Tentar match parcial
            for c in vcf_chroms:
                if norm_chrom(c) == chrom_norm:
                    chrom_vcf = c
                    break

        if chrom_vcf is None:
            print(f"    [WARN] Cromossomo '{chrom_norm}' não encontrado no VCF "
                  f"(contigs: {list(vcf_chroms)[:5]}...)", file=sys.stderr)
            continue

        try:
            for rec in tbx.fetch(chrom_vcf, start - 1, end):  # tabix é 0-based
                parts = rec.split("\t")
                if len(parts) < 4:
                    continue
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue
                ref = parts[3].strip()
                if len(ref) == 1 and ref.upper() in "ACGT":
                    positions[pos] = ref.upper()
        except ValueError as e:
            print(f"    [WARN] fetch erro para {chrom_vcf}:{start}-{end}: {e}",
                  file=sys.stderr)

    tbx.close()
    return positions


# ═══════════════════════════════════════════════════════════════════════════
# Extração linear (fallback sem pysam) — UMA passagem por cromossomo
# ═══════════════════════════════════════════════════════════════════════════

def extract_positions_linear(
    vcf_path: Path,
    regions: list[dict],
    progress_interval: int = 500_000
) -> dict[int, str]:
    """
    Lê o VCF UMA VEZ linearmente e coleta posições REF nas regiões de interesse.
    Mais lento que tabix mas não requer pysam.
    """
    intervals = [(r["start"], r["end"]) for r in regions]
    min_pos   = min(s for s, e in intervals)
    max_pos   = max(e for s, e in intervals)

    positions: dict[int, str] = {}
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    lines_read = 0

    print(f"  → [linear] Lendo {vcf_path.name} "
          f"(intervalo: {min_pos:,}–{max_pos:,}) ...", flush=True)

    try:
        with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                lines_read += 1

                if lines_read % progress_interval == 0:
                    print(f"    {lines_read:,} linhas lidas, "
                          f"{len(positions):,} posições coletadas ...", flush=True)

                parts = line.split("\t", 5)
                if len(parts) < 4:
                    continue
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue

                if pos < min_pos or pos > max_pos:
                    continue

                ref = parts[3].strip()
                if len(ref) != 1 or ref.upper() not in "ACGT":
                    continue

                for s, e in intervals:
                    if s <= pos <= e:
                        positions[pos] = ref.upper()
                        break

    except (OSError, EOFError, zlib.error, Exception) as ex:
        print(f"  [WARN] Erro ao ler {vcf_path.name}: {ex}", file=sys.stderr)

    print(f"  → {lines_read:,} linhas processadas, "
          f"{len(positions):,} posições REF coletadas.", flush=True)
    return positions


# ═══════════════════════════════════════════════════════════════════════════
# Pipeline principal
# ═══════════════════════════════════════════════════════════════════════════

def run(regions_by_chrom: dict[str, list[dict]],
        vcf_dir: Path,
        out_path: Path,
        min_density: float,
        use_tabix: bool = True) -> dict:

    out_path.parent.mkdir(parents=True, exist_ok=True)
    stats = {
        "chroms_processed": 0,
        "chroms_no_vcf": 0,
        "sequences_written": 0,
        "sequences_skipped_density": 0,
        "total_regions": sum(len(v) for v in regions_by_chrom.values()),
    }

    # Verificar disponibilidade de pysam
    has_pysam = False
    if use_tabix:
        try:
            import pysam  # noqa: F401
            has_pysam = True
            print("[INFO] pysam disponível — usando tabix para queries regionais rápidas.")
        except ImportError:
            print("[INFO] pysam não disponível — usando leitura linear (mais lento).")
            print("  Para instalar: pip3 install pysam --break-system-packages")

    start_time = datetime.datetime.now()
    print(f"\n[{start_time:%H:%M:%S}] Iniciando extração")
    print(f"  Cromossomos a processar : {len(regions_by_chrom)}")
    print(f"  Total de regiões        : {stats['total_regions']}")
    print(f"  Densidade mínima        : {min_density:.0%}")
    print(f"  Modo                    : {'tabix (rápido)' if has_pysam else 'linear (lento)'}")
    print(f"  Saída                   : {out_path}\n")

    with open(out_path, "w", encoding="ascii") as fa:
        for chrom_norm, regions in sorted(regions_by_chrom.items(),
                                          key=lambda x: x[0].zfill(3)):
            vcf_path = find_vcf_for_chrom(vcf_dir, chrom_norm)
            if vcf_path is None:
                print(f"[WARN] Sem VCF para chr{chrom_norm} "
                      f"({len(regions)} regiões ignoradas)", file=sys.stderr)
                stats["chroms_no_vcf"] += 1
                continue

            stats["chroms_processed"] += 1
            t0 = datetime.datetime.now()
            print(f"[{t0:%H:%M:%S}] Cromossomo {chrom_norm} "
                  f"({len(regions)} regiões) ...", flush=True)

            # Extrair posições — tabix se disponível, senão linear
            if has_pysam:
                # Detectar prefixo de cromossomo no VCF
                chrom_prefix = detect_chrom_prefix(vcf_path)
                # Garantir índice tabix
                if ensure_tabix_index(vcf_path):
                    all_positions = extract_positions_tabix(
                        vcf_path, regions, chrom_prefix
                    )
                    print(f"  → [tabix] {len(all_positions):,} posições REF coletadas.",
                          flush=True)
                else:
                    # Fallback para linear se tabix falhar
                    all_positions = extract_positions_linear(vcf_path, regions)
            else:
                all_positions = extract_positions_linear(vcf_path, regions)

            # Reconstruir e escrever cada sequência
            written_chr = 0
            skipped_density = 0
            for reg in regions:
                seq   = reconstruct_sequence(all_positions, reg["start"], reg["end"])
                known = sum(seq.count(b) for b in "ACGT")
                density = known / max(len(seq), 1)

                if density < min_density:
                    skipped_density += 1
                    stats["sequences_skipped_density"] += 1
                    if skipped_density <= 3:
                        print(f"  [skip] {reg['label'][:60]}: "
                              f"densidade {density:.1%} < {min_density:.0%} "
                              f"({known}/{len(seq)} bases conhecidas)", flush=True)
                    continue

                fa.write(f">{reg['label']}\n")
                for i in range(0, len(seq), 60):
                    fa.write(seq[i:i+60] + "\n")
                stats["sequences_written"] += 1
                written_chr += 1

            elapsed = (datetime.datetime.now() - t0).total_seconds()
            print(f"  ✓ {written_chr}/{len(regions)} escritas, "
                  f"{skipped_density} ignoradas (baixa densidade) "
                  f"({elapsed:.1f}s)\n", flush=True)

    total_elapsed = (datetime.datetime.now() - start_time).total_seconds()
    print("=" * 60)
    print(f"[CONCLUÍDO em {total_elapsed/60:.1f} min]")
    print(f"  Sequências escritas      : {stats['sequences_written']}")
    print(f"  Ignoradas (baixa dens.)  : {stats['sequences_skipped_density']}")
    print(f"  Cromossomos sem VCF      : {stats['chroms_no_vcf']}")
    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"  FASTA salvo em           : {out_path}")
        print(f"  Tamanho                  : {out_path.stat().st_size / 1024:.1f} KB")
    print("=" * 60)

    return stats


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--mode", choices=["dmcs", "genes"], default="genes",
                    help="'genes' = regiões dos genes-alvo (permanente, recomendado); "
                         "'dmcs' = janelas ao redor dos DMCs (específico por run)")
    ap.add_argument("--json", default=None,
                    help="JSON de genes-alvo (obrigatório no modo 'genes')")
    ap.add_argument("--dmcs", default=None,
                    help="CSV de DMCs da Fase 1 (obrigatório no modo 'dmcs')")
    ap.add_argument("--vcf-dir", required=True,
                    help="Diretório com VCFs do 1000 Bull Genomes (Chr*-Run8-TAUIND*.vcf.gz)")
    ap.add_argument("--out", default=None,
                    help="Caminho do FASTA de saída (padrão: deduzido do modo)")
    ap.add_argument("--window", type=int, default=500,
                    help="Janela em bp ao redor do DMC/gene quando flanco não definido. "
                         "Padrão: 500")
    ap.add_argument("--min-density", type=float, default=0.01,
                    help="Fração mínima de bases conhecidas (não-N). "
                         "Padrão: 0.01 (1%%) — VCFs contêm apenas SNPs, "
                         "não todas as bases, então densidade sempre será baixa.")
    ap.add_argument("--no-tabix", action="store_true",
                    help="Desabilitar uso de pysam/tabix e usar leitura linear.")
    args = ap.parse_args()

    vcf_dir = Path(args.vcf_dir)
    if not vcf_dir.is_dir():
        sys.exit(f"[ERRO] Diretório VCF não encontrado: {vcf_dir}")

    # Coletar regiões conforme o modo
    if args.mode == "genes":
        if not args.json:
            sys.exit("[ERRO] --json é obrigatório no modo 'genes'")
        if not Path(args.json).exists():
            sys.exit(f"[ERRO] JSON não encontrado: {args.json}")

        if args.out is None:
            json_path = Path(args.json)
            stem = json_path.stem.replace("_target_genes_coordinates", "")
            args.out = str(json_path.parent / f"{stem}_reference_sequences.fa")

        print(f"[Modo: genes] JSON: {args.json}")
        regions_by_chrom = collect_regions_from_json(args.json, args.window)

    else:  # mode == "dmcs"
        if not args.dmcs:
            sys.exit("[ERRO] --dmcs é obrigatório no modo 'dmcs'")
        if not Path(args.dmcs).exists():
            sys.exit(f"[ERRO] DMCs não encontrado: {args.dmcs}")

        if args.out is None:
            args.out = str(Path(args.dmcs).parent.parent / "fase4_homer" / "dmcs_sequences.fa")

        print(f"[Modo: dmcs] CSV: {args.dmcs}")
        regions_by_chrom = collect_regions_from_dmcs(args.dmcs, args.window)

    if not regions_by_chrom:
        sys.exit("[ERRO] Nenhuma região coletada. Verifique o arquivo de entrada.")

    out_path = Path(args.out)
    stats = run(
        regions_by_chrom,
        vcf_dir,
        out_path,
        args.min_density,
        use_tabix=not args.no_tabix
    )

    if stats["sequences_written"] == 0:
        print("\n[WARN] Nenhuma sequência gerada.", file=sys.stderr)
        print("  Causas possíveis:", file=sys.stderr)
        print("  1. Cromossomos do JSON/CSV não correspondem aos nomes nos VCFs", file=sys.stderr)
        print("     (ex: '21' vs 'Chr21' vs 'chr21')", file=sys.stderr)
        print("  2. Density < min-density: VCFs contêm apenas SNPs (~1 por kb),", file=sys.stderr)
        print("     então sequências longas terão densidade < 1%.", file=sys.stderr)
        print("     Use --min-density 0.0 para aceitar qualquer densidade.", file=sys.stderr)
        print("  3. Coordenadas dos genes fora do intervalo do VCF.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
