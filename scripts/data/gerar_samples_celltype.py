#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gera os sample sheets dos contrastes resolvidos por tipo celular (C8-C14)
do dataset ZENODO19555520: Nelore (Bos indicus, treatment=0) vs
Holstein-Friesian (Bos taurus, treatment=1), n=3 vs n=3 em cada tipo celular.

Mesma estrutura/colunas do C1 (samples_C1_zenodo_NEL_vs_HF.csv), apenas
trocando LysedBlood pelo tipo celular imune. Os arquivos FILTERED_*.cov.gz
(janela +/-15 kb, painel de 73 genes) ja existem em
E:\\data\\epigenomas\\ZENODO19555520.

Saida: subprojeto/SUB02/samples/samples_C{8..14}_*.csv
"""
from pathlib import Path

HERE = Path(__file__).resolve().parents[2]            # epigenoma/
OUT = HERE / "subprojeto" / "SUB02" / "samples"

NELORE = ["ERS4410323", "ERS4410324", "ERS4410325"]   # Bos indicus, controle (0)
HF = ["ERS4410317", "ERS4410318", "ERS4410319"]       # Bos taurus, teste (1)

# tipo celular -> (id do contraste)
CELLTYPES = {
    "CD4":   "C8",
    "CD8":   "C9",
    "Bcell": "C10",
    "gdT":   "C11",
    "gran":  "C12",
    "mono":  "C13",
    "NK":    "C14",
}

HEADER = ("sampleID,breed,species,species_code,tissue_cell,treatment,"
          "condition,season,dataset,file")


def row(sample, breed, species, code, ct, treat, fname):
    return (f"{sample},{breed},{species},{code},{ct},{treat},"
            f"control,NA,ZENODO19555520,{fname}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    written = []
    for ct, cid in CELLTYPES.items():
        lines = [HEADER]
        for acc in NELORE:
            f = f"FILTERED_Nelore_{acc}_{ct}.bismark.cov.gz"
            lines.append(row(f"nelore_{acc}_{ct}", "Nelore",
                             "Bos indicus", "indicus", ct, 0, f))
        for acc in HF:
            f = f"FILTERED_HF_{acc}_{ct}.bismark.cov.gz"
            lines.append(row(f"hf_{acc}_{ct}", "Holstein-Friesian",
                             "Bos taurus", "taurus", ct, 1, f))
        name = f"samples_{cid}_zenodo_NEL_vs_HF_{ct}.csv"
        (OUT / name).write_text("\n".join(lines) + "\n", encoding="utf-8")
        written.append(name)
        print(f"[OK] {name}  ({len(lines)-1} amostras: 3 Nelore + 3 HF, {ct})")
    print(f"\nTotal: {len(written)} sample sheets em {OUT}")


if __name__ == "__main__":
    main()
