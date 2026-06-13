#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Orquestrador do Pipeline GenOmics Epigenoma
===========================================
Executa as 8 fases do projeto em sequência, com logs centralizados:

    Fase 0 — Pré-processamento (filtragem por coordenadas)
    Fase 1 — Análise diferencial (methylKit, R)
    Fase 2 — Interseção genômica (BEDTools-like puro Python)
    Fase 3 — Perfis TSS / heatmaps (deepTools-like)
    Fase 4 — Motifs (HOMER — requer binário externo; opcional)
    Fase 5 — Anotação ChIPseeker (R)
    Fase 6 — Integração multi-ômica (MOFA2 — R; opcional)
    Fase 7 — GLM / DSS (R)
    Fase 8 — Consolidação e relatório final

Uso:
    python run_pipeline.py \
        --metadata samples.csv \
        --data-dir /data/epigenomas/NCBI_GSE161113 \
        --out results/run_$(date +%Y%m%d)

O orquestrador detecta automaticamente quais fases podem rodar (baseado em
dependências instaladas) e pula silenciosamente as demais, registrando o
motivo. Isso garante que o pipeline **nunca trave** por falta de um pacote.

Autor: Luiz Dione Barbosa De-Melo — GenOmics Epigenoma 2026
"""
from __future__ import annotations
import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent  # scripts/pipeline/
SCRIPTS_DIR = HERE.parent                # scripts/
ROOT = SCRIPTS_DIR.parent                # raiz do projeto (epigenoma/)
R_DIR = SCRIPTS_DIR / "R"                # scripts/R/ — onde estão os .R

# ── Resolução do JSON de genes-alvo ─────────────────────────────────────────
# Prioridade (decrescente):
#   1. subprojeto/<ID>/<ID>_target_genes_coordinates.json  (subprojeto ativo)
#   2. genes_alvo/all_target_genes_coordinates.json        (cache global)
#   3. all_target_genes_coordinates.json                   (raiz, legado)
#
# O subprojeto ativo é detectado automaticamente:
#   - Lê subprojeto/subprojeto_config.json → campo "genes_file" → extrai ID
#   - Ou varre subprojeto/<ID>/<ID>_target_genes_coordinates.json
#   - Pode ser sobrescrito via --json na CLI

def _detect_subproject_json(root: Path) -> Path | None:
    """Detecta automaticamente o JSON do subprojeto ativo."""
    sp_dir = root / "subprojeto"
    if not sp_dir.exists():
        return None

    # Tentar ler subprojeto_config.json para pegar o genes_file
    cfg_path = sp_dir / "subprojeto_config.json"
    if cfg_path.exists():
        try:
            with open(cfg_path, encoding="utf-8") as fh:
                cfg = json.load(fh)
            genes_file = cfg.get("genes_file", "")
            if genes_file:
                candidate = sp_dir / genes_file if not genes_file.startswith("/") else Path(genes_file)
                # genes_file pode ser "SUB001/SUB001_target_genes_coordinates.json"
                if not candidate.exists():
                    candidate = sp_dir / genes_file.lstrip("SUB001/")
                if not candidate.exists():
                    # Tentar relativo à raiz do projeto
                    candidate = root / genes_file
                if candidate.exists():
                    return candidate
        except Exception:
            pass

    # Varrer subpastas de subprojeto/ procurando <ID>/<ID>_target_genes_coordinates.json
    best: Path | None = None
    best_count = 0
    for sub_id_dir in sorted(sp_dir.iterdir()):
        if not sub_id_dir.is_dir():
            continue
        sid = sub_id_dir.name
        candidate = sub_id_dir / f"{sid}_target_genes_coordinates.json"
        if candidate.exists():
            try:
                with open(candidate, encoding="utf-8") as fh:
                    data = json.load(fh)
                n = len(data.get("genes", {}))
                if n > best_count:  # preferir o que tem mais genes (mais completo)
                    best_count = n
                    best = candidate
            except Exception:
                pass
    return best


_subproject_json = _detect_subproject_json(ROOT)
_genes_alvo_json = ROOT / "genes_alvo" / "all_target_genes_coordinates.json"

if _subproject_json is not None:
    TARGETS_JSON = _subproject_json
elif _genes_alvo_json.exists():
    TARGETS_JSON = _genes_alvo_json
else:
    TARGETS_JSON = ROOT / "all_target_genes_coordinates.json"

# Integração Genomas: diretório dos VCFs do 1000 Bull Genomes
# Mapeado no Docker como /data/1000bg_vcf (E:\data\1000bg_vcf no host)
VCF_DIR = Path("/data/1000bg_vcf")


def run(cmd: list[str], log_file: Path, cwd: str | None = None) -> tuple[int, str]:
    """Executa comando e registra stdout+stderr em log_file. Retorna (rc, output)."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with open(log_file, "w", encoding="utf-8") as lf:
        lf.write(f"$ {' '.join(cmd)}\n\n")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                cwd=cwd, text=True, bufsize=1)
        out_buf = []
        assert proc.stdout
        for line in proc.stdout:
            lf.write(line)
            out_buf.append(line)
            print(line, end="")
        proc.wait()
        return proc.returncode, "".join(out_buf)


def has(binary: str) -> bool:
    if shutil.which(binary) is not None:
        return True
    # Paths absolutos comuns em containers rocker/rstudio e Debian/Ubuntu
    fallback_paths = [
        f"/usr/local/bin/{binary}",
        f"/usr/bin/{binary}",
        f"/usr/lib/R/bin/{binary}",
        f"/opt/R/bin/{binary}",
    ]
    return any(os.path.isfile(p) and os.access(p, os.X_OK) for p in fallback_paths)


def rscript_cmd() -> str:
    """Retorna o comando/path correto do Rscript."""
    for candidate in ["Rscript",
                      "/usr/local/bin/Rscript",
                      "/usr/bin/Rscript",
                      "/usr/lib/R/bin/Rscript"]:
        if shutil.which(candidate) or (os.path.isfile(candidate) and os.access(candidate, os.X_OK)):
            return candidate if shutil.which(candidate) else candidate
    return "Rscript"


def phase_report(name: str, status: str, detail: str = "") -> dict:
    return {"phase": name, "status": status, "detail": detail,
            "timestamp": datetime.datetime.now().isoformat(timespec="seconds")}


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", help="samples.csv: sampleID,breed,treatment,file")
    ap.add_argument("--data-dir", default="/data/epigenomas/NCBI_GSE161113",
                    help="Pasta com arquivos de metilação (.gz) — default: /data/epigenomas/NCBI_GSE161113")
    ap.add_argument("--out", default=f"results/run_{datetime.date.today():%Y%m%d}")
    ap.add_argument("--skip", default="", help="Lista CSV de fases a pular (ex.: 4,6)")
    ap.add_argument("--json", default=None,
                    help="JSON de genes-alvo (sobrescreve detecção automática do subprojeto)")
    ap.add_argument("--fasta", default=None,
                    help="FASTA de referência pré-gerado pelo extract_sequences_from_vcf.py "
                         "(opcional — se ausente, HOMER roda sem sequência)")
    args = ap.parse_args()

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    log_dir = out / "logs"; log_dir.mkdir(exist_ok=True)
    report: list[dict] = []
    skip = {s.strip() for s in args.skip.split(",") if s.strip()}

    # Resolver JSON final: argumento CLI > detecção automática
    global TARGETS_JSON
    if args.json:
        TARGETS_JSON = Path(args.json)
        if not TARGETS_JSON.exists():
            print(f"[WARN] --json especificado não encontrado: {TARGETS_JSON}")

    # Log do JSON que será usado em todas as fases
    try:
        with open(TARGETS_JSON, encoding="utf-8") as fh:
            _n_genes = len(json.load(fh).get("genes", {}))
    except Exception:
        _n_genes = "?"
    print(f"[Config] JSON de genes-alvo: {TARGETS_JSON}")
    print(f"[Config] Genes carregados   : {_n_genes}")

    # ---------------- Fase 0 — Pré-processamento
    if "0" in skip:
        report.append(phase_report("0_preprocess", "SKIPPED", "por solicitação"))
    else:
        script0 = ROOT / "baixar_e_filtrar_metilacao.py"
        if script0.exists():
            rc, _ = run([sys.executable, str(script0)], log_dir / "fase0.log")
            report.append(phase_report("0_preprocess", "OK" if rc == 0 else "FAIL", f"rc={rc}"))
        else:
            report.append(phase_report("0_preprocess", "SKIPPED", "script não encontrado"))

    # ---------------- Fase 1 — methylKit
    _rscript = rscript_cmd()
    if args.metadata and has("Rscript") and "1" not in skip:
        r1 = R_DIR / "run_methylkit_diff.R"  # scripts/R/run_methylkit_diff.R
        rc, _ = run(
            [_rscript, str(r1), "--metadata", args.metadata,
             "--data-dir", args.data_dir, "--out", str(out / "fase1_methylkit")],
            log_dir / "fase1.log",
        )
        report.append(phase_report("1_methylkit", "OK" if rc == 0 else "FAIL", f"rc={rc}"))
        dmc_csv = out / "fase1_methylkit" / "DMCs_significant.csv"
    else:
        report.append(phase_report("1_methylkit", "SKIPPED",
                                   "Rscript/metadata indisponível ou pulado"))
        # Fix: se fase1 foi pulada mas o CSV já existe de run anterior, reusar
        _existing_dmc = out / "fase1_methylkit" / "DMCs_significant.csv"
        dmc_csv = _existing_dmc if _existing_dmc.exists() else None
        if dmc_csv:
            print(f"[Info] Fase 1 pulada — reusando DMCs existentes: {dmc_csv}")

    # ---------------- Fase 2 — Interseção com alvos
    if dmc_csv and dmc_csv.exists() and "2" not in skip:
        rc, _ = run(
            [sys.executable, str(HERE / "intersect_dmrs_genes.py"),
             "--dmrs", str(dmc_csv), "--json", str(TARGETS_JSON),
             "--out", str(out / "fase2_intersect" / "dmrs_in_target_genes.csv")],
            log_dir / "fase2.log",
        )
        report.append(phase_report("2_intersect", "OK" if rc == 0 else "FAIL", f"rc={rc}"))
    else:
        report.append(phase_report("2_intersect", "SKIPPED", "sem DMCs da fase 1"))

    # ---------------- Fase 3 — Perfis TSS (precisa arquivos das amostras)
    if args.metadata and "3" not in skip:
        try:
            import csv as _csv
            with open(args.metadata, newline="") as fh:
                reader = _csv.DictReader(fh)
                rows = list(reader)
            methyl_files = [os.path.join(args.data_dir, r["file"]) for r in rows]
            breeds = ",".join(r["breed"] for r in rows)
            rc, _ = run(
                [sys.executable, str(HERE / "tss_methylation_profile.py"),
                 "--methylation", *methyl_files,
                 "--breeds", breeds, "--json", str(TARGETS_JSON),
                 "--window", "10000",
                 "--out", str(out / "fase3_deeptools")],
                log_dir / "fase3.log",
            )
            report.append(phase_report("3_tss_profile", "OK" if rc == 0 else "FAIL", f"rc={rc}"))
        except Exception as e:
            report.append(phase_report("3_tss_profile", "FAIL", str(e)))
    else:
        report.append(phase_report("3_tss_profile", "SKIPPED", "sem metadata"))

    # ---------------- Fase 4 — HOMER / Motifs (run_homer_motifs_v2 com genoma completo)
    if "4" not in skip:
        bed = out / "fase1_methylkit" / "DMCs_significant.bed"
        homer_script_v2 = HERE / "run_homer_motifs_v2.py"
        homer_script_v1 = HERE / "run_homer_motifs.py"
        homer_script = homer_script_v2 if homer_script_v2.exists() else homer_script_v1

        # FASTA do genoma completo ARS-UCD1.3 (montado em /data no container)
        GENOME_FASTA = Path("/data/genome_reference/ARS-UCD1.3/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa")

        if bed.exists() and homer_script.exists():
            fase4_out = out / "fase4_homer"

            # Sub-fase 4a: localizar FASTA — prioridade:
            #   0. Genoma completo ARS-UCD1.3 em /data (PREFERIDO — v2 usa pysam direto)
            #   1. --fasta passado na linha de comando
            #   2. FASTA permanente ao lado do JSON do subprojeto
            #   3. Sem FASTA → somente se usando v1 (modo fallback)
            dmc_fasta = None
            fasta_source = "nenhum"

            if GENOME_FASTA.exists() and GENOME_FASTA.with_suffix(".fa.fai").exists():
                dmc_fasta = str(GENOME_FASTA)
                fasta_source = "genoma completo ARS-UCD1.3"
                print(f"[Fase 4a] Usando genoma completo: {dmc_fasta}")
            elif args.fasta and Path(args.fasta).exists():
                dmc_fasta = args.fasta
                fasta_source = "argumento --fasta"
                print(f"[Fase 4a] FASTA fornecido via --fasta: {dmc_fasta}")
            else:
                # Procurar FASTA permanente ao lado do JSON do subprojeto
                json_path = Path(str(TARGETS_JSON))
                stem = json_path.stem.replace("_target_genes_coordinates", "")
                fasta_permanente = json_path.parent / f"{stem}_reference_sequences.fa"
                if fasta_permanente.exists() and fasta_permanente.stat().st_size > 100:
                    dmc_fasta = str(fasta_permanente)
                    fasta_source = "FASTA subprojeto"
                    print(f"[Fase 4a] FASTA permanente encontrado: {fasta_permanente}")
                elif GENOME_FASTA.exists() and not GENOME_FASTA.with_suffix(".fa.fai").exists():
                    print(f"[Fase 4a] Genoma encontrado mas sem índice .fai — indexando...")
                    print(f"          Execute: samtools faidx {GENOME_FASTA}")
                    dmc_fasta = str(GENOME_FASTA)
                    fasta_source = "genoma completo ARS-UCD1.3 (sem .fai)"
                else:
                    print(f"[Fase 4a] FASTA não encontrado.")
                    print(f"          Baixe: python3 scripts/data/baixar_genoma_referencia.py")
                    print(f"          Índice: samtools faidx {GENOME_FASTA}")

            # Sub-fase 4b: análise de motifs
            # v2 requer --fasta; v1 aceita sem
            use_v2 = homer_script == homer_script_v2
            if use_v2 and dmc_fasta is None:
                # Sem FASTA e v2 requer: registrar e pular
                report.append(phase_report("4_homer", "SKIPPED",
                    "run_homer_motifs_v2.py requer --fasta; genoma ARS-UCD1.3 não indexado"))
            else:
                cmd4 = [sys.executable, str(homer_script),
                        "--bed", str(bed), "--out", str(fase4_out)]
                if dmc_fasta:
                    cmd4 += ["--fasta", dmc_fasta]
                rc, _ = run(cmd4, log_dir / "fase4.log")
                detail = f"rc={rc} | FASTA: {fasta_source} | script: {homer_script.name}"
                report.append(phase_report("4_homer", "OK" if rc == 0 else "FAIL", detail))
        elif not bed.exists():
            report.append(phase_report("4_homer", "SKIPPED", "BED da fase 1 ausente"))
        else:
            report.append(phase_report("4_homer", "SKIPPED", "run_homer_motifs_v2.py não encontrado"))
    else:
        report.append(phase_report("4_homer", "SKIPPED", "por solicitação"))

    # ---------------- Fase 5 + 7 — Anotação + GLM
    if has("Rscript") and dmc_csv and dmc_csv.exists() and "5" not in skip:
        r57 = R_DIR / "run_annotation_glm.R"
        cmd = [_rscript, str(r57), "--dmrs", str(dmc_csv),
               "--json", str(TARGETS_JSON), "--out", str(out / "fase5_fase7")]
        if args.metadata:
            cmd += ["--metadata", args.metadata]
        rc, _ = run(cmd, log_dir / "fase57.log")
        report.append(phase_report("5_7_annotation_glm", "OK" if rc == 0 else "FAIL", f"rc={rc}"))
    else:
        report.append(phase_report("5_7_annotation_glm", "SKIPPED", "Rscript ou DMCs ausentes"))

    # ---------------- Fase 6 — MOFA2 (opcional, exige R + vcfs)
    if has("Rscript") and dmc_csv and dmc_csv.exists() and "6" not in skip:
        mofa_script = R_DIR / "prepare_mofa2.R"
        if mofa_script.exists() and args.metadata:
            cmd6 = [_rscript, str(mofa_script),
                    "--methylation", str(dmc_csv),
                    "--metadata", args.metadata,
                    "--json", str(TARGETS_JSON),
                    "--data-dir", args.data_dir,
                    "--out", str(out / "fase6_mofa2")]
            rc, _ = run(cmd6, log_dir / "fase6.log")
            report.append(phase_report("6_mofa2", "OK" if rc == 0 else "FAIL", f"rc={rc}"))
        else:
            report.append(phase_report("6_mofa2", "SKIPPED",
                                       "prepare_mofa2.R ou metadata ausente"))
    else:
        report.append(phase_report("6_mofa2", "SKIPPED",
                                   "Rscript/DMCs ausentes ou pulado"))

    # ---------------- Fase 8 — Relatório consolidado
    summary = {
        "run_date": datetime.datetime.now().isoformat(timespec="seconds"),
        "data_dir": args.data_dir,
        "metadata": args.metadata,
        "output_dir": str(out),
        "phases": report,
    }
    with open(out / "pipeline_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    with open(out / "pipeline_summary.txt", "w", encoding="utf-8") as fh:
        fh.write("GenOmics Epigenoma — resumo do run\n")
        fh.write("=" * 50 + "\n")
        for ph in report:
            fh.write(f"[{ph['status']:<7}] {ph['phase']:<22} {ph['detail']}\n")

    # ---------------- Relatório HTML (se generate_report.py existir)
    report_script = HERE / "generate_report.py"
    if report_script.exists():
        cmd_rep = [sys.executable, str(report_script),
                   "--results-dir", str(out),
                   "--json", str(TARGETS_JSON),
                   "--out", str(out / "report_epigenoma.html")]
        if args.metadata:
            cmd_rep += ["--metadata", args.metadata]
        rc, _ = run(cmd_rep, log_dir / "report.log")
        if rc == 0:
            print(f"\n📊 Relatório HTML: {out / 'report_epigenoma.html'}")

    print("=" * 60)
    print("Pipeline finalizado. Resumo:")
    for ph in report:
        icon = "OK" if ph["status"] == "OK" else ("FAIL" if ph["status"] == "FAIL" else "SKIP")
        print(f"  [{icon}] {ph['phase']:<24} {ph['status']:<8} {ph['detail']}")

    summary = {
        "run_date": datetime.datetime.now().isoformat(timespec="seconds"),
        "data_dir": args.data_dir,
        "metadata": args.metadata,
        "output_dir": str(out),
        "phases": report,
    }
    with open(out / "pipeline_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    with open(out / "pipeline_summary.txt", "w", encoding="utf-8") as fh:
        fh.write("GenOmics Epigenoma - resumo do run\n")
        fh.write("=" * 50 + "\n")
        for ph in report:
            fh.write(f"[{ph['status']:<7}] {ph['phase']:<22} {ph['detail']}\n")

    print(f"\nOutputs em: {out}")
    print(f"Resumo JSON: {out / 'pipeline_summary.json'}")


if __name__ == "__main__":
    main()
