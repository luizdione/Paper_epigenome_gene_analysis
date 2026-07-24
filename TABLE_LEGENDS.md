# Table legends

**Table 1** — Differentially methylated cytosines (DMCs) per contrast (methylKit, q < 0.05) and annotated target genes. Discovery and systematisation use the same 108-gene panel and the same uniform ±15 kb regulatory window around each gene; discovery is the June Angus-versus-Nelore comparison (equivalent to systematisation contrast C3). Per-contrast counts are in `tables/TableS2_by_comparison.tsv`.

**Table 2** — High-tier candidate genes (priority score ≥ 5) from the multi-omic bridge (ten genes), with regulatory methylation (promoter/upstream Δβ), best heat-stress expression and orthogonal fetal-liver corroboration. Full per-gene data in `tables/TableS1_bridge_gene_table.tsv`.

## Supplementary tables (`tables/`)
- **Table S1** — `TableS1_bridge_gene_table.tsv`: bridge gene set, the 47 genes differentially marked at both omic layers, with per-gene methylation × expression signals, concordance and priority score.
- **Table S1-extended** — `TableS1ext_bridge_gene_table_FULL.tsv`: complete bridge table, every gene evaluated, with the `in_methylation` / `in_expression` flags from which Table S1 is derived.
- **Table S2** — `TableS2_by_comparison.tsv`: concordant bridge genes decomposed by contributing contrast.
- **Table S3** — `TableS3_module_enrichment.tsv`: functional-module over-representation statistics (hypergeometric; universe = 108-gene panel; no module survives Benjamini–Hochberg).
- **Table S4** — `TableS4_C6_validation.tsv`: fetal-liver WGBS corroboration per candidate (25 of 47 carry a DMR; 15 in the same taurine-versus-indicine direction).
- **Table S5** — `TableS5_cissnp.tsv`: cis-genetic classification of the 2719 subspecies DMCs against the 1000 Bull Genomes catalogue (CpG-SNP / near-SNP / no-SNP).
- **Table S6** — `TableS6_genomic_context.tsv`: genomic-context (Fst / ancestry-informative-marker) status of the bridge candidates (Gir versus Holstein).
- **Table S7** — `TableS7_matched_correlation.tsv`: matched-sample promoter-methylation-versus-expression Spearman correlations per panel gene (Zhou et al. 2020 atlas; 88 genes).
- **Table S8** — `TableS8_motif_dmc.tsv`: subspecies DMCs co-localised (±25 bp) with motif-disrupting SNPs for the enriched GC-box factors SP1 and KLF4.
- **Table S10** — `TableS10_motif_enrichment.tsv`: transcription-factor-motif enrichment at the discovery DMC flanks (43 motifs tested, with fold enrichment over genomic expectation, sequence counts and nominal / BH p-values; six motifs with multi-sequence support (KLF4, NF-κB p65, NFAT1, NF-κB p50, HSF1, SP1) survive q < 0.05).
