#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_pipeline.py — Validação end-to-end do pipeline com dados sintéticos
=============================================================================
Gera dados de metilação sintéticos realistas para os 45 genes-alvo e executa
as Fases 2 e 3 do pipeline, verificando que os scripts produzem outputs
corretos sem necessidade de dados reais do SSD.

Também gera um DMC CSV sintético (simulando saída do methylKit) para testar
a interseção genômica e perfis TSS.

Uso:
    python scripts/validate_pipeline.py [--out test_validation]
    python scripts/validate_pipeline.py --phases 2,3 --n-samples 4

Autor: Luiz Dione Barbosa De-Melo — GenOmics Epigenoma 2026
"""
from __future__ import annotations
import argparse
import csv
import gzip
import json
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
_subprojeto_json = ROOT / "subprojeto" / "all_target_genes_coordinates.json"
TARGETS_JSON = _subprojeto_json if _subprojeto_json.exists() else ROOT / "all_target_genes_coordinates.json"


def load_genes(json_path: Path) -> list[dict]:
    """Carrega genes do JSON de alvos."""
    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    genes = []
    for name, info in data.get("genes", {}).items():
        genes.append({
            "gene": name,
            "chrom": str(info.get("chromosome", "")),
            "start": info.get("start"),
            "end": info.get("end"),
            "flanking_start": info.get("flanking_start"),
            "flanking_end": info.get("flanking_end"),
            "strand": info.get("strand", 1),
            "module": info.get("module", ""),
        })
    return genes


def generate_synthetic_methylation(
    genes: list[dict],
    n_cpgs_per_gene: int = 50,
    breed: str = "Angus",
    seed: int = 42,
) -> list[list]:
    """
    Gera dados de metilação sintéticos no formato methylKit:
    chr, start, end, strand, coverage, numCs, numTs

    Angus (taurus) → metilação média ~60-75%
    Nelore (indicus) → metilação média ~45-60% (hipometilação em termotolerância)
    """
    rng = random.Random(seed)
    rows = []

    # Ajustar média de metilação por raça para simular diferenças reais
    if breed in ("Angus", "Holstein", "HF"):
        mean_meth = 68.0
    else:  # Nelore, ND, Gir
        mean_meth = 52.0

    for g in genes:
        if g["start"] is None or g["end"] is None:
            continue
        chrom = g["chrom"]
        # Gerar CpGs dentro da região flanqueada (incluindo flancos ±10kb)
        flank_start = g.get("flanking_start", g["start"] - 10000)
        flank_end = g.get("flanking_end", g["end"] + 10000)

        # Módulos de termotolerância: maior diferença entre raças
        module_boost = 0
        if g["module"] in ("thermotolerance", "thermotolerance_extra"):
            module_boost = 8 if breed in ("Angus", "Holstein", "HF") else -8

        for _ in range(n_cpgs_per_gene):
            pos = rng.randint(flank_start, flank_end)
            coverage = rng.randint(5, 80)
            # Metilação com variação por módulo
            local_meth = max(0, min(100, rng.gauss(mean_meth + module_boost, 15)))
            numC = max(0, min(coverage, round(coverage * local_meth / 100)))
            numT = coverage - numC
            freqC = round(100 * numC / coverage, 2) if coverage > 0 else 0
            strand = "F" if rng.random() > 0.5 else "R"
            rows.append([chrom, pos, pos, strand, coverage, numC, numT, freqC])

    # Ordenar por cromossomo e posição
    rows.sort(key=lambda r: (r[0], r[1]))
    return rows


def write_methylation_file(rows: list[list], path: str, compressed: bool = True):
    """Escreve arquivo de metilação no formato methylKit."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    if compressed:
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write("chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n")
            for r in rows:
                chrom, pos, _, strand, cov, numC, numT, freqC = r
                fh.write(f"{chrom}.{pos}\t{chrom}\t{pos}\t{strand}\t{cov}\t{freqC}\t{round(100-freqC,2)}\n")
    else:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n")
            for r in rows:
                chrom, pos, _, strand, cov, numC, numT, freqC = r
                fh.write(f"{chrom}.{pos}\t{chrom}\t{pos}\t{strand}\t{cov}\t{freqC}\t{round(100-freqC,2)}\n")


def generate_synthetic_dmcs(genes: list[dict], n_dmcs: int = 200, seed: int = 42) -> list[dict]:
    """
    Gera DMCs sintéticos (simulando saída methylKit) para testar Fase 2.
    Alguns DMCs caem dentro das regiões alvo, outros são intergênicos.
    """
    rng = random.Random(seed)
    dmcs = []

    # 60% dos DMCs dentro de regiões alvo
    n_in_target = int(n_dmcs * 0.6)
    for i in range(n_in_target):
        g = rng.choice(genes)
        if g["start"] is None:
            continue
        flank_start = g.get("flanking_start", g["start"] - 10000)
        flank_end = g.get("flanking_end", g["end"] + 10000)
        pos = rng.randint(flank_start, flank_end)
        diff = rng.choice([-1, 1]) * rng.uniform(15, 40)
        qval = 10 ** (-rng.uniform(2, 10))
        dmcs.append({
            "chr": g["chrom"],
            "start": pos,
            "end": pos,
            "strand": "*",
            "pvalue": qval * 10,
            "qvalue": qval,
            "meth.diff": round(diff, 2),
        })

    # 40% intergênicos (cromossomos aleatórios)
    chroms = list(set(g["chrom"] for g in genes if g["chrom"]))
    for i in range(n_dmcs - n_in_target):
        chrom = rng.choice(chroms)
        pos = rng.randint(1_000_000, 100_000_000)
        diff = rng.choice([-1, 1]) * rng.uniform(15, 35)
        qval = 10 ** (-rng.uniform(2, 8))
        dmcs.append({
            "chr": chrom,
            "start": pos,
            "end": pos,
            "strand": "*",
            "pvalue": qval * 10,
            "qvalue": qval,
            "meth.diff": round(diff, 2),
        })

    dmcs.sort(key=lambda d: (d["chr"], d["start"]))
    return dmcs


def write_dmcs_csv(dmcs: list[dict], path: str):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["chr", "start", "end", "strand", "pvalue", "qvalue", "meth.diff"])
        writer.writeheader()
        writer.writerows(dmcs)


def write_samples_csv(sample_files: list[dict], path: str):
    """Gera samples.csv para o pipeline."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sampleID", "breed", "species", "treatment", "file"])
        writer.writeheader()
        writer.writerows(sample_files)


def run_phase(name: str, cmd: list[str], cwd: str | None = None) -> bool:
    """Executa um comando e reporta sucesso/falha."""
    print(f"\n{'='*60}")
    print(f"[VALIDAÇÃO] {name}")
    print(f"  $ {' '.join(cmd)}")
    print(f"{'='*60}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd, timeout=120)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr)
        if result.returncode == 0:
            print(f"[✓] {name} — PASSOU")
            return True
        else:
            print(f"[✗] {name} — FALHOU (rc={result.returncode})")
            return False
    except Exception as e:
        print(f"[✗] {name} — ERRO: {e}")
        return False


def validate_output(path: str, min_lines: int = 1) -> bool:
    """Verifica se um arquivo de saída existe e tem conteúdo mínimo."""
    if not os.path.exists(path):
        print(f"  [✗] Arquivo esperado não encontrado: {path}")
        return False
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()
    if len(lines) < min_lines:
        print(f"  [✗] Arquivo tem apenas {len(lines)} linhas (esperado >= {min_lines}): {path}")
        return False
    print(f"  [✓] OK: {path} ({len(lines)} linhas)")
    return True


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="test_validation", help="Pasta de saída para validação")
    ap.add_argument("--phases", default="2,3", help="Fases a validar (CSV: 2,3)")
    ap.add_argument("--n-samples", type=int, default=4, help="Nº de amostras sintéticas por grupo")
    ap.add_argument("--n-cpgs", type=int, default=50, help="CpGs por gene por amostra")
    ap.add_argument("--n-dmcs", type=int, default=200, help="Nº total de DMCs sintéticos")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    phases = {p.strip() for p in args.phases.split(",")}

    print("[Info] Carregando genes de alvos...")
    genes = load_genes(TARGETS_JSON)
    print(f"[Info] {len(genes)} genes carregados")

    results = {}

    # --- Gerar DMCs sintéticos para Fase 2 ---
    if "2" in phases:
        print("\n[Geração] Criando DMCs sintéticos...")
        dmcs = generate_synthetic_dmcs(genes, n_dmcs=args.n_dmcs)
        dmc_path = str(out / "synthetic_DMCs.csv")
        write_dmcs_csv(dmcs, dmc_path)
        print(f"  → {len(dmcs)} DMCs gerados em {dmc_path}")

        # Rodar Fase 2
        intersect_out = str(out / "fase2_intersect" / "dmrs_in_target_genes.csv")
        ok = run_phase(
            "Fase 2 — Interseção genômica",
            [sys.executable, str(HERE / "intersect_dmrs_genes.py"),
             "--dmrs", dmc_path,
             "--json", str(TARGETS_JSON),
             "--out", intersect_out],
        )
        results["fase2"] = ok
        if ok:
            validate_output(intersect_out, min_lines=2)

    # --- Gerar amostras sintéticas para Fase 3 ---
    if "3" in phases:
        print("\n[Geração] Criando amostras de metilação sintéticas...")
        sample_dir = out / "synthetic_samples"
        sample_dir.mkdir(parents=True, exist_ok=True)

        sample_entries = []
        methyl_files = []
        breeds_list = []

        groups = [("Angus", args.n_samples), ("Nelore", args.n_samples)]
        for breed, n in groups:
            for i in range(n):
                seed = hash(f"{breed}_{i}") % (2**31)
                rows = generate_synthetic_methylation(genes, n_cpgs_per_gene=args.n_cpgs,
                                                      breed=breed, seed=seed)
                fname = f"SYNTH_{breed}_{i+1:02d}_CpG.txt.gz"
                fpath = str(sample_dir / fname)
                write_methylation_file(rows, fpath, compressed=True)
                methyl_files.append(fpath)
                breeds_list.append(breed)
                treatment = 0 if breed == "Angus" else 1
                sample_entries.append({
                    "sampleID": f"{breed.lower()}_{i+1:02d}",
                    "breed": breed,
                    "species": "Bos taurus" if breed == "Angus" else "Bos indicus",
                    "treatment": treatment,
                    "file": fname,
                })
                print(f"  → {fname} ({len(rows)} CpGs)")

        # samples.csv
        samples_path = str(out / "samples_synthetic.csv")
        write_samples_csv(sample_entries, samples_path)
        print(f"  → samples.csv: {samples_path}")

        # Rodar Fase 3
        fase3_out = str(out / "fase3_deeptools")
        ok = run_phase(
            "Fase 3 — Perfis TSS (deepTools-like)",
            [sys.executable, str(HERE / "tss_methylation_profile.py"),
             "--methylation", *methyl_files,
             "--breeds", ",".join(breeds_list),
             "--json", str(TARGETS_JSON),
             "--out", fase3_out],
        )
        results["fase3"] = ok
        if ok:
            validate_output(os.path.join(fase3_out, "genes_ordered.csv"), min_lines=2)
            for breed in ("Angus", "Nelore"):
                validate_output(os.path.join(fase3_out, f"matrix_{breed}.csv"), min_lines=2)

    # --- Resumo ---
    print("\n" + "=" * 60)
    print("RESUMO DA VALIDAÇÃO")
    print("=" * 60)
    all_ok = True
    for phase, ok in sorted(results.items()):
        marker = "✓" if ok else "✗"
        print(f"  {marker} {phase}: {'PASSOU' if ok else 'FALHOU'}")
        if not ok:
            all_ok = False

    if all_ok:
        print("\n[✓] Todas as fases validadas com sucesso!")
    else:
        print("\n[✗] Algumas fases falharam. Verifique os logs acima.")

    print(f"\nOutputs em: {out}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
