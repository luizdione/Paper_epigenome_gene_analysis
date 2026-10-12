# Thermotolerance gene set analysis

Comparative DNA-methylation analysis of indicine (*Bos taurus indicus*) versus taurine
(*Bos taurus taurus*) cattle. The pipeline reuses public reduced-representation bisulfite
sequencing (RRBS), an orthogonal whole-genome bisulfite layer, an independent heat-stress
transcriptome and a taurine-indicine ancestry-differentiation layer (Gyr versus Holstein
F<sub>ST</sub> / ancestry-informative markers computed from the 1000 Bull Genomes catalogue) to
prioritise a compact thermotolerance / immune / somatotropic candidate gene set, and tests
whether the prioritised signal is epigenetic rather than a footprint of sequence divergence.

This repository holds the **computational scripts, the run outputs, the result tables and
the result figures (with legends)** of the analysis.

## Layout
- `scripts/`
  - `R/` — methylKit differential methylation; ChIPseeker annotation; MOFA2
  - `python/` — multi-omic bridge (`bridge_deg_dmr.py`, `bridge_enrichment.py`); cis-SNP
    classification (`cissnp_annotation.py`); motif–DMC crossing (`motif_dmc_cross.py`);
    tri-layer F<sub>ST</sub>/AIM context (`triple_layer.py`); matched-atlas correlation
    (`matched_atlas_bridge.py`); orthogonal-WGBS validation (`c6_validation.py`);
    `export_online_resources.py` (supplementary tables from `results/`);
    `figures/` (result figures `Fig01`–`Fig06`, drawn at final print size)
  - `pipeline/` — DMR–gene intersection; TF-motif scan; TSS methylation profiles; orchestration
  - `data/` — retrieval/preparation of the public datasets and sample sheets
  - `analysis/`, `utils/` — sequence extraction, genotype matrices, QC
- `results/` — run outputs: `DISCOVERY_108_15kb/` (discovery run) and `SUB02_runs_108/`
  (13-contrast systematisation, multi-omic bridge, cis-SNP, motif, matched-atlas and
  fetal-liver layers; `AGGREGATE_108.json` is the machine-readable run summary);
  `ancestry/degs_x_aims_genes.tsv` (gene-level Gyr-versus-Holstein F<sub>ST</sub>/AIM status)
- `tables/` — supplementary result tables, one per Online Resource (`OR01`–`OR11`, `.tsv`)
- `figures/` — result figures `Fig01`–`Fig06` (`.png`, 600 dpi, and `.pdf`)
- `FIGURE_LEGENDS.md`, `TABLE_LEGENDS.md` — legends for the figures and tables

## Data sources (public)
- **GSE161113** — RRBS of Angus and Nellore peripheral blood (hot and thermoneutral seasons)
- **Zenodo 19555520** — Bismark cytosine reports for Holstein-Friesian, N'Dama and Nellore
  across whole blood and seven flow-sorted immune cell types
- **MacPhillamy et al. (2024), GigaScience 13:giae061** — Angus-versus-Brahman fetal-liver WGBS DMRs
- **1000 Bull Genomes Project (Run 8)** — taurine/indicine variant catalogue (cis-SNP test;
  Gyr-versus-Holstein F<sub>ST</sub> and ancestry-informative markers)
- **Zhou et al. (2020)** — matched WGBS + RNA-seq tissue atlas (one Holstein cow)
- **PRJNA885909, GSE235171, PRJNA562925** — public bovine heat-stress RNA-seq (hypothalamus and
  pituitary, granulosa cells, PBMC), reprocessed with nf-core/rnaseq + Salmon + DESeq2

## Notes on reproduction
- Breed codes and sample identifiers inside the scripts and sample sheets follow the source
  datasets as downloaded (e.g. `Nelore`, `HF`, `ND`, `A01`/`N01`); the result tables and
  figures use the spellings Nellore and Gyr.
- Download scripts in `scripts/data/` write to a data root given on the command line (`--dest`)
  or by their `BASE_DIR` constant; analysis scripts read raw data from that root and write
  to `results/`. The `high_tier` column in the run outputs under `results/` reflects the
  working prioritisation at run time; `export_online_resources.py` applies the final ten-gene
  high-tier set when it writes `tables/`.
- Figures: `cd scripts/python/figures && python make_fig01_discovery_108.py` (and the other
  `make_fig*` scripts); profile `print` (default) draws at 174 mm width.

## Licence
Code is released under the MIT licence (see `LICENSE`); result tables and figures are
released under CC BY 4.0.
