#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
probe_vcfs_urls.py — Sonda HTTP HEAD em URLs candidatas para VCFs alternativos.

Não baixa nada — apenas verifica quais URLs existem (HTTP 200) e tamanho.
Útil para mapear fontes ANTES de iniciar download de GBs de dados.

Uso:
  python probe_vcfs_urls.py
"""
from __future__ import annotations
import sys
import urllib.request
import urllib.error

URLS = [
    # ── Ensembl Variation (Bos taurus, sem genótipos por animal) ────────
    ("Ensembl all",       "https://ftp.ensembl.org/pub/release-113/variation/vcf/bos_taurus/bos_taurus.vcf.gz"),
    ("Ensembl all (112)", "https://ftp.ensembl.org/pub/release-112/variation/vcf/bos_taurus/bos_taurus.vcf.gz"),
    ("Ensembl all (110)", "https://ftp.ensembl.org/pub/release-110/variation/vcf/bos_taurus/bos_taurus.vcf.gz"),
    ("Ensembl Chr16",     "https://ftp.ensembl.org/pub/release-113/variation/vcf/bos_taurus/bos_taurus_chromosome-16.vcf.gz"),
    ("Ensembl Chr22",     "https://ftp.ensembl.org/pub/release-113/variation/vcf/bos_taurus/bos_taurus_chromosome-22.vcf.gz"),
    ("Ensembl Chr25",     "https://ftp.ensembl.org/pub/release-113/variation/vcf/bos_taurus/bos_taurus_chromosome-25.vcf.gz"),
    ("Ensembl Chr1",      "https://ftp.ensembl.org/pub/release-113/variation/vcf/bos_taurus/bos_taurus_chromosome-1.vcf.gz"),  # controle

    # ── NCBI dbSNP Bos taurus archive ───────────────────────────────────
    ("NCBI dbSNP all",     "https://ftp.ncbi.nih.gov/snp/organisms/archive/bovine/VCF/00-All.vcf.gz"),
    ("NCBI bovine VCF",    "https://ftp.ncbi.nih.gov/snp/organisms/bos_taurus_9913/VCF/"),

    # ── EBI EVA (European Variation Archive) ────────────────────────────
    ("EBI EVA cattle",     "https://ftp.ebi.ac.uk/eva/datasets/cattle/"),
    ("EBI EVA Run8 base",  "https://ftp.ebi.ac.uk/eva/PRJEB42513/"),

    # ── 1000 Bull Genomes Run7 (versão anterior) ────────────────────────
    ("1000bg Run7 base",   "https://ftp.sra.ebi.ac.uk/vol1/ena_ftp/PRJEB23303/"),

    # ── Cattle Variation Database Roslin ────────────────────────────────
    ("Roslin BovReg",      "https://www.bovreg.eu/data/"),

    # ── Mirror UCSC ─────────────────────────────────────────────────────
    ("UCSC bovine SNP",    "https://hgdownload.soe.ucsc.edu/gbdb/bosTau9/"),

    # ── FAANG portal ────────────────────────────────────────────────────
    ("FAANG cattle",       "https://data.faang.org/api/file/_search?q=cattle+vcf"),
]


def probe(name, url, timeout=20):
    try:
        req = urllib.request.Request(url, method="HEAD",
                                     headers={"User-Agent": "epigenoma-probe/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            cl = r.headers.get("Content-Length", "?")
            ct = r.headers.get("Content-Type", "?")
            return f"HTTP {r.status} | size={cl} | type={ct[:30]}"
    except urllib.error.HTTPError as e:
        return f"HTTP {e.code} {e.reason}"
    except urllib.error.URLError as e:
        return f"URL_ERR: {str(e.reason)[:60]}"
    except Exception as e:
        return f"{type(e).__name__}: {str(e)[:60]}"


print("=" * 90)
print("  PROBE de URLs candidatas para VCFs alternativos (HTTP HEAD apenas)")
print("=" * 90)
print()
print(f"{'Fonte':<25} {'Status':<55}")
print("-" * 90)

for name, url in URLS:
    status = probe(name, url)
    print(f"{name:<25} {status}")

print()
print("=" * 90)
print("  Próximos passos:")
print("    Use as URLs com HTTP 200 + tamanho razoável (>100MB para genomas completos)")
print("    Para baixar uma fonte específica: python baixar_vcfs_alternativos.py --source <nome>")
print("=" * 90)
