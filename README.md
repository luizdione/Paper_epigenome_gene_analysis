# Thermotolerance gene set analysis

Comparative DNA-methylation analysis of indicine (*Bos taurus indicus*) versus taurine
(*Bos taurus taurus*) cattle. The pipeline reuses public reduced-representation bisulfite
sequencing (RRBS), an orthogonal whole-genome bisulfite layer, and an independent
heat-stress transcriptome to prioritise a compact thermotolerance / immune / somatotropic
candidate gene set, and tests whether the prioritised signal is epigenetic rather than a
footprint of sequence divergence.

This repository holds the **computational scripts, the result tables, and the result
figures (with legends)** of the analysis.

## Layout
- `scripts/`
  - `R/` — methylKit differential methylation; ChIPseeker annotation; MOFA2
  - `python/` — multi-omic bridge; cis-SNP classification; motif–DMC crossing;
    tri-layer (Fst/AIM); matched-atlas correlation; orthogonal-WGBS validation; result figures
  - `pipeline/` — DMR–gene intersection; TF-motif scan; TSS methylation profiles; orchestration
  - `data/` — retrieval/preparation of the public datasets and sample sheets
  - `analysis/`, `utils/` — sequence extraction, genotype matrices, QC
- `tables/` — supplementary tables S1–S8 (`.tsv`) and a consolidated workbook (`.xlsx`)
- `figures/` — result figures `Fig01`–`Fig06` (`.png`)
- `FIGURE_LEGENDS.md`, `TABLE_LEGENDS.md` — legends for the figures and tables

## Data sources (public)
- **GSE161113** — RRBS of Angus and Nelore peripheral blood
- **Zenodo 19555520** — Bismark cytosine reports for Holstein-Friesian, N'Dama and Nelore
  across whole blood and seven flow-sorted immune cell types
- **MacPhillamy et al. (2024), GigaScience 13:giae061** — Angus-versus-Brahman fetal-liver WGBS DMRs
- **1000 Bull Genomes Project (Run 8)** — taurine/indicine variant catalogue (cis-SNP test)
- **Zhou et al. (2020), GEO GSE147087 + GSE137943** — matched WGBS + RNA-seq taurine tissue atlas
- Heat-stress transcriptome: companion RNA-seq differentially-expressed-gene sets (four tissues)

Reference assembly: **ARS-UCD1.3**.

## Methods (summary)
methylKit (logistic regression, overdispersion correction, Benjamini–Hochberg FDR) for
differential methylation; ChIPseeker for genomic annotation; a position-weight-matrix scan
for transcription-factor-motif enrichment; MOFA2 for unsupervised structure; and a gene-level
methylation–expression concordance bridge with feature-weighted prioritisation. Robustness is
assessed against an orthogonal WGBS layer, the 1000 Bull Genomes variant catalogue,
Gir-versus-Holstein differentiation, and a matched WGBS+RNA-seq atlas.

All reported values are read from versioned result files; no value is hand-entered.

## License
Released for academic reuse. Please cite the underlying public datasets listed above.
