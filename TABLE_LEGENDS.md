# Table legends

**Table 1** — Differentially methylated cytosines (DMCs) per contrast (methylKit, q < 0.05)
and annotated target genes. The discovery contrast uses the 53-gene panel with ±10 kb windows;
the systematisation contrasts use the 73-gene panel with ±15 kb windows. (Per-contrast counts
are in `tables/TableS2_by_comparison.tsv`.)

**Table 2** — High-tier candidate genes (priority score ≥ 5) from the multi-omic bridge, with
regulatory methylation, best heat-stress expression and orthogonal fetal-liver corroboration.
(Full per-gene data in `tables/TableS1_bridge_gene_table.tsv`.)

## Supplementary tables (`tables/`)
- **Table S1** — `TableS1_bridge_gene_table.tsv`: bridge gene set — the 31 genes differentially
  marked at both omic layers, with per-gene methylation × expression signals, concordance and
  priority score.
- **Table S1-extended** — `TableS1ext_bridge_gene_table_FULL.tsv`: complete bridge table — every
  gene evaluated, with the `in_methylation` / `in_expression` flags from which Table S1 is derived.
- **Table S2** — `TableS2_by_comparison.tsv`: concordant bridge genes decomposed by contributing contrast.
- **Table S3** — `TableS3_module_enrichment.tsv`: functional-module enrichment statistics.
- **Table S4** — `TableS4_C6_validation.tsv`: fetal-liver WGBS corroboration per candidate.
- **Table S5** — `TableS5_cissnp.tsv`: cis-genetic classification of subspecies DMCs against the
  1000 Bull Genomes catalogue (CpG-SNP / near-SNP / no-SNP).
- **Table S6** — `TableS6_genomic_context.tsv`: genomic-context (Fst / ancestry-informative-marker)
  status of the bridge candidates.
- **Table S7** — `TableS7_matched_correlation.tsv`: matched-sample promoter-methylation-versus-
  expression Spearman correlations per panel gene (Zhou et al. 2020 atlas).
- **Table S8** — `TableS8_motif_dmc.tsv`: subspecies DMCs co-localised (±25 bp) with
  motif-disrupting SNPs for the enriched GC-box factors SP1 and KLF4.

A consolidated workbook with all sheets is provided as `tables/Supplementary_Tables_S1-S8.xlsx`.
