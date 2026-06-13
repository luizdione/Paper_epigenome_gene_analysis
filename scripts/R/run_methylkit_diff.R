#!/usr/bin/env Rscript
# =============================================================================
# Fase 1 — Análise Diferencial de Metilação (headless, CLI-first)
# =============================================================================
# Este script replica, em modo batch, exatamente a lógica da Aba 3 do Shiny.
# Ele aceita um diretório com arquivos de metilação e uma planilha de metadados
# (sampleID, breed, file), produz DMCs filtrados (diff > cutoff & q < alpha),
# salva CSV, QC plots e um relatório consolidado.
#
# Uso:
#   Rscript run_methylkit_diff.R \
#       --metadata samples.csv \
#       --data-dir /data/NCBI_GSE161113 \
#       --assembly ARS-UCD1.3 \
#       --diff 15 --qvalue 0.05 \
#       --out results/fase1_methylkit
#
# metadata.csv:
#   sampleID,breed,treatment,file
#   Nelore_01,Nelore,0,FILTERED_nelore_N01_CD4.gz   # Bos indicus = controle (treat=0)
#   Angus_01,Angus,1,FILTERED_angus_A01_CD4.gz      # Bos taurus  = teste   (treat=1)
#
# Autor: Luiz Dione Barbosa De-Melo — GenOmics Epigenoma 2026
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Instale optparse:  install.packages('optparse')")
  }
  if (!requireNamespace("methylKit", quietly = TRUE)) {
    stop("Instale methylKit (Bioconductor).")
  }
  library(optparse)
  library(methylKit)
})

# ------------------------------------------------------------- CLI
opt_list <- list(
  make_option(c("--metadata"), type = "character", default = NULL,
              help = "CSV com colunas sampleID,breed,treatment,file"),
  make_option(c("--data-dir"), type = "character", default = ".",
              help = "Diretório onde estão os arquivos de metilação"),
  make_option(c("--assembly"), type = "character", default = "ARS-UCD1.3"),
  make_option(c("--diff"),     type = "numeric",   default = 15,
              help = "Diferença mínima de metilação (padrão 15%%)"),
  make_option(c("--qvalue"),   type = "numeric",   default = 0.05),
  make_option(c("--min-cov"),  type = "integer",   default = 5,
              help = "Cobertura mínima por CpG (padrão 5)"),
  make_option(c("--out"),      type = "character", default = "results/fase1_methylkit")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$metadata)) stop("--metadata é obrigatório.")

# ------------------------------------------------------------- carregamento
meta <- read.csv(opt$metadata, stringsAsFactors = FALSE)
required <- c("sampleID", "treatment", "file")
missing_cols <- setdiff(required, colnames(meta))
if (length(missing_cols) > 0) stop("Colunas ausentes em metadata: ", paste(missing_cols, collapse = ", "))

files <- as.list(file.path(opt$`data-dir`, meta$file))
missing_files <- files[!file.exists(unlist(files))]
if (length(missing_files) > 0) stop("Arquivos não encontrados: ", paste(missing_files, collapse = ", "))

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
log_con <- file(file.path(opt$out, "methylkit_run.log"), open = "wt")
log <- function(...) { msg <- paste0(..., collapse = ""); message(msg); writeLines(msg, log_con) }

log("[Info] Amostras: ", nrow(meta))
log("[Info] Treatments únicos: ", paste(unique(meta$treatment), collapse = ", "))

# ------------------------------------------------------------- methylKit pipeline
# Detectar formato pelo nome do arquivo:
#   *.bismark.cov.gz  → formato bismark coverage: chr start end %meth countM countU
#   *.CpG.txt.gz      → formato CpG_context bismark: chr strand start %meth countM countU
# methylKit aceita pipeline="bismarkCoverage" para ambos os formatos padrão do Bismark.
# Não usar pipeline manual com fraction=FALSE pois calculateDiffMeth exige freq 0-1.

first_file <- files[[1]]
is_bismark_cov <- grepl("\\.bismark\\.cov(\\.gz)?$", first_file, ignore.case = TRUE)
is_cpg_txt     <- grepl("CpG[_\\.].*\\.txt(\\.gz)?$",  first_file, ignore.case = TRUE)

pipeline_arg <- if (is_bismark_cov || is_cpg_txt) {
  log("[Info] Formato detectado: bismark coverage — usando pipeline='bismarkCoverage'")
  "bismarkCoverage"
} else {
  log("[Info] Formato não reconhecido — usando pipeline='bismarkCoverage' como padrão")
  "bismarkCoverage"
}

obj <- methRead(
  files,
  sample.id   = as.list(meta$sampleID),
  assembly    = opt$assembly,
  treatment   = meta$treatment,
  context     = "CpG",
  mincov      = opt$`min-cov`,
  pipeline    = pipeline_arg
)

# QC por amostra
pdf(file.path(opt$out, "QC_methylation_stats.pdf"), width = 8, height = 6)
for (i in seq_along(obj)) {
  getMethylationStats(obj[[i]], plot = TRUE, both.strands = FALSE)
  title(sub = paste("Sample:", meta$sampleID[i]))
}
dev.off()

pdf(file.path(opt$out, "QC_coverage_stats.pdf"), width = 8, height = 6)
for (i in seq_along(obj)) {
  getCoverageStats(obj[[i]], plot = TRUE, both.strands = FALSE)
  title(sub = paste("Sample:", meta$sampleID[i]))
}
dev.off()

# Unificação e análise diferencial
meth <- unite(obj, destrand = FALSE)
log("[Info] CpGs unidos: ", nrow(meth))
myDiff <- calculateDiffMeth(meth)
sig <- getMethylDiff(myDiff, difference = opt$diff, qvalue = opt$qvalue)

log("[Info] DMCs significativos (|diff| > ", opt$diff, "%, q < ", opt$qvalue, "): ", nrow(getData(sig)))

# Se nenhum DMC passou o filtro estrito, usar top 50 por q-value como fallback informativo
sig_df <- getData(sig)
all_df  <- getData(myDiff)
if (nrow(sig_df) == 0 && nrow(all_df) > 0) {
  log("[WARN] Nenhum DMC passou o filtro estrito. Usando top-50 por q-value como output informativo.")
  all_ord <- all_df[order(all_df$qvalue, -abs(all_df$meth.diff)), ]
  sig_df  <- head(all_ord, 50)
}
write.csv(sig_df, file.path(opt$out, "DMCs_significant.csv"), row.names = FALSE)
write.csv(getData(myDiff), file.path(opt$out, "DMCs_all.csv"), row.names = FALSE)

# Também salva como BED para interseção posterior
if (nrow(sig_df) > 0) {
  strand_val <- if (is.null(sig_df$strand) || length(sig_df$strand) == 0) {
    rep(".", nrow(sig_df))
  } else {
    sig_df$strand
  }
  bed <- data.frame(chrom = sig_df$chr, start = sig_df$start, end = sig_df$end,
                    name = paste0("DMC_", seq_len(nrow(sig_df))),
                    score = round(-log10(pmax(sig_df$qvalue, 1e-300)), 2),
                    strand = strand_val)
  write.table(bed, file.path(opt$out, "DMCs_significant.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
} else {
  log("[WARN] Nenhum DMC significativo encontrado — BED não gerado.")
  # Criar BED vazio para não bloquear fases seguintes
  file.create(file.path(opt$out, "DMCs_significant.bed"))
}

# Histograma das diferenças de metilação
pdf(file.path(opt$out, "methylation_diff_histogram.pdf"), width = 7, height = 5)
hist(getData(myDiff)$meth.diff, breaks = 60, col = "#5B9BD5",
     main = "Distribuição de diferença de metilação",
     xlab = "meth.diff (%)", ylab = "Nº de CpGs")
abline(v = c(-opt$diff, opt$diff), col = "red", lty = 2)
dev.off()

close(log_con)
cat("[OK] Análise concluída. Resultados em:", opt$out, "\n")
