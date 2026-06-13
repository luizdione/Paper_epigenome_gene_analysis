#!/usr/bin/env Rscript
# =============================================================================
# run_mofa2_multiview.R
# MOFA2 multi-view: methylation (Epigenoma) + genotype (1000 Bull Genomes)
# =============================================================================

# CRÍTICO: setar PYTHONPATH antes de qualquer importação do reticulate
# para que mofapy2 e dependências (scipy, numpy) sejam encontradas
Sys.setenv(PYTHONPATH = paste(
  "/usr/local/lib/python3.10/dist-packages",
  "/usr/lib/python3/dist-packages",
  sep = ":"
))
reticulate::use_python("/usr/bin/python3", required = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(MOFA2)
})

PROJ    <- "/home/rstudio/projeto"
OUT_DIR <- file.path(PROJ, "scripts/results/SUB001_GSE161113_run_oficial/fase6_mofa2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== MOFA2 Multi-view: Methylation + Genotype ===\n")
cat(format(Sys.time(), "[%H:%M:%S]"), "\n\n")

# ─── 1. Carregar methylation matrix ─────────────────────────────────────────
meth_csv <- file.path(OUT_DIR, "methylation_matrix_multisamples.csv")
cat("Carregando methylation matrix:", meth_csv, "\n")
meth_mat_raw <- read.csv(meth_csv, row.names = 1, check.names = FALSE)
cat("  Dimensões brutas:", nrow(meth_mat_raw), "amostras ×", ncol(meth_mat_raw), "genes\n")

# Transpor para features × amostras (MOFA2 espera features nas linhas)
meth_mat <- t(as.matrix(meth_mat_raw))  # genes × amostras
cat("  Para MOFA2 (features × amostras):", nrow(meth_mat), "×", ncol(meth_mat), "\n")

# ─── 2. Carregar genotype matrix ─────────────────────────────────────────────
geno_csv <- file.path(OUT_DIR, "genotype_matrix.csv")
cat("\nCarregando genotype matrix:", geno_csv, "\n")
geno_mat_raw <- read.csv(geno_csv, row.names = 1, check.names = FALSE)
cat("  Dimensões brutas:", nrow(geno_mat_raw), "amostras ×", ncol(geno_mat_raw), "SNPs\n")

# Garantir mesma ordem de amostras que meth_mat
common_samples <- intersect(colnames(meth_mat), rownames(geno_mat_raw))
cat("  Amostras comuns:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("Nenhuma amostra em comum entre methylation e genotype matrices!")
}

meth_mat  <- meth_mat[, common_samples, drop = FALSE]
geno_mat  <- t(as.matrix(geno_mat_raw[common_samples, ]))  # SNPs × amostras
cat("  Genotype para MOFA2 (SNPs × amostras):", nrow(geno_mat), "×", ncol(geno_mat), "\n")

# ─── 3. Verificar completude dos dados ───────────────────────────────────────
cat("\n--- QC das matrizes ---\n")
cat("Methylation NAs:", sum(is.na(meth_mat)), "/", length(meth_mat),
    sprintf("(%.1f%%)\n", 100*sum(is.na(meth_mat))/length(meth_mat)))
cat("Genotype NAs:  ", sum(is.na(geno_mat)), "/", length(geno_mat),
    sprintf("(%.1f%%)\n", 100*sum(is.na(geno_mat))/length(geno_mat)))

# Remover SNPs com NAs (todos são breed-level, não deve haver NAs)
geno_na_rows <- rowSums(is.na(geno_mat)) > 0
if (sum(geno_na_rows) > 0) {
  cat("  Removendo", sum(geno_na_rows), "SNPs com NA\n")
  geno_mat <- geno_mat[!geno_na_rows, , drop = FALSE]
}

# Remover SNPs com variância zero (frequência idêntica em Angus e Nelore)
geno_var <- apply(geno_mat, 1, var, na.rm = TRUE)
zero_var <- geno_var == 0
if (sum(zero_var) > 0) {
  cat("  Removendo", sum(zero_var), "SNPs com variância zero\n")
  geno_mat <- geno_mat[!zero_var, , drop = FALSE]
}
cat("  SNPs após QC:", nrow(geno_mat), "\n")

# ─── 4. Selecionar top SNPs mais informativos (reduzir dimensionalidade) ──────
# Selecionar os SNPs com maior diferença entre Angus e Nelore
# A matriz já tem freq_angus para Angus samples e freq_nelore para Nelore
# Como todos Angus têm o mesmo valor, variância entre amostras ≠ 0 mas
# variância reflete diferença Angus vs Nelore
MAX_SNPS <- 500  # máximo SNPs para MOFA2 (eficiência computacional)
if (nrow(geno_mat) > MAX_SNPS) {
  geno_var2 <- apply(geno_mat, 1, var, na.rm = TRUE)
  top_idx   <- order(geno_var2, decreasing = TRUE)[1:MAX_SNPS]
  geno_mat  <- geno_mat[top_idx, , drop = FALSE]
  cat("  Top", MAX_SNPS, "SNPs selecionados por variância\n")
}

cat("\nDimensões finais para MOFA2:\n")
cat("  methylation:", nrow(meth_mat), "features ×", ncol(meth_mat), "amostras\n")
cat("  genotype:   ", nrow(geno_mat), "features ×", ncol(geno_mat), "amostras\n")

# ─── 5. Criar objeto MOFA2 ───────────────────────────────────────────────────
mofa_data <- list(
  methylation = meth_mat,
  genotype    = geno_mat
)

cat("\nCriando objeto MOFA2 multi-view...\n")
mofa_obj <- create_mofa(mofa_data)
print(mofa_obj)

# ─── 6. Configurar opções ────────────────────────────────────────────────────
data_opts  <- get_default_data_options(mofa_obj)

model_opts <- get_default_model_options(mofa_obj)
# Número de fatores: min(5, n_samples - 1)
model_opts$num_factors <- min(5, length(common_samples) - 1)
cat("\nNúmero de fatores latentes:", model_opts$num_factors, "\n")

train_opts <- get_default_training_options(mofa_obj)
train_opts$maxiter          <- 1000
train_opts$convergence_mode <- "slow"
train_opts$seed             <- 42
train_opts$verbose          <- TRUE
train_opts$startELBO        <- 10
train_opts$freqELBO         <- 5

mofa_obj <- prepare_mofa(mofa_obj,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)

# ─── 7. Treinar ──────────────────────────────────────────────────────────────
hdf5_path <- file.path(OUT_DIR, "mofa2_multiview_model.hdf5")
cat("\nTreinando MOFA2 multi-view...\n")

tryCatch({
  sp_ver <- reticulate::import("scipy")$version$full_version
  np_ver <- reticulate::import("numpy")$version$full_version
  cat("  scipy:", sp_ver, "| numpy:", np_ver, "\n")
}, error = function(e) cat("  Versão check:", e$message, "\n"))

mofa_obj <- run_mofa(mofa_obj, outfile = hdf5_path, use_basilisk = FALSE)

# ─── 8. Resultados ───────────────────────────────────────────────────────────
cat("\n=== RESULTADOS MOFA2 MULTI-VIEW ===\n")
print(mofa_obj)

# Variância explicada por view
cat("\nVariância explicada por fator e view:\n")
var_exp <- calculate_variance_explained(mofa_obj)
print(var_exp$r2_per_factor)

cat("\nVariância total explicada por view:\n")
print(var_exp$r2_total)

# Salvar modelo
saveRDS(mofa_obj, file.path(OUT_DIR, "mofa2_multiview_trained.rds"))

# ─── 9. Plots ────────────────────────────────────────────────────────────────
# Variância explicada (2 views)
pdf(file.path(OUT_DIR, "mofa2_multiview_variance.pdf"), width = 10, height = 6)
tryCatch(
  plot_variance_explained(mofa_obj, x = "view", y = "factor"),
  error = function(e) cat("Plot variância:", e$message, "\n")
)
dev.off()

# Fatores latentes
pdf(file.path(OUT_DIR, "mofa2_multiview_factors.pdf"), width = 12, height = 8)
tryCatch({
  nf <- dim(get_factors(mofa_obj)[[1]])[2]
  plot_factors(mofa_obj, factors = 1:min(nf, 5))
}, error = function(e) cat("Plot fatores:", e$message, "\n"))
dev.off()

# Weights methylation
pdf(file.path(OUT_DIR, "mofa2_multiview_weights_meth.pdf"), width = 10, height = 6)
tryCatch(
  plot_weights(mofa_obj, view = "methylation", factor = 1, nfeatures = 10),
  error = function(e) cat("Plot weights meth:", e$message, "\n")
)
dev.off()

# Weights genotype
pdf(file.path(OUT_DIR, "mofa2_multiview_weights_geno.pdf"), width = 10, height = 6)
tryCatch(
  plot_weights(mofa_obj, view = "genotype", factor = 1, nfeatures = 15),
  error = function(e) cat("Plot weights geno:", e$message, "\n")
)
dev.off()

# ─── 10. Exportar tabelas ────────────────────────────────────────────────────
# Fatores por amostra
tryCatch({
  factors_df <- as.data.frame(get_factors(mofa_obj)[[1]])
  factors_df$sampleID <- rownames(factors_df)

  # Adicionar metadata
  meta_file <- file.path(PROJ, "samples_GSE161113_demo.csv")
  meta      <- read.csv(meta_file, stringsAsFactors = FALSE)
  factors_df <- merge(factors_df,
                      meta[, c("sampleID","breed","condition","treatment","month")],
                      by = "sampleID", all.x = TRUE)
  write.csv(factors_df,
            file.path(OUT_DIR, "mofa2_multiview_factors_metadata.csv"),
            row.names = FALSE)
  cat("Fatores com metadata salvos\n")
}, error = function(e) cat("Fatores:", e$message, "\n"))

# Weights de cada view
for (view_name in c("methylation", "genotype")) {
  tryCatch({
    w <- get_weights(mofa_obj, views = view_name, factors = "all", as.data.frame = TRUE)
    w <- w[order(abs(w$value), decreasing = TRUE), ]
    out_f <- file.path(OUT_DIR, paste0("mofa2_multiview_weights_", view_name, ".csv"))
    write.csv(w, out_f, row.names = FALSE)
    cat("Weights", view_name, "salvos:", nrow(w), "features\n")
  }, error = function(e) cat("Weights", view_name, ":", e$message, "\n"))
}

# Resumo correlações entre fatores e covariáveis
cat("\n--- Correlação fatores × breed/condition ---\n")
tryCatch({
  f <- get_factors(mofa_obj)[[1]]
  meta_sub <- meta[match(rownames(f), meta$sampleID), ]
  for (col in c("breed", "condition", "treatment")) {
    if (col %in% names(meta_sub)) {
      vals <- as.numeric(as.factor(meta_sub[[col]]))
      for (k in 1:ncol(f)) {
        r <- cor(f[, k], vals, use = "complete.obs")
        if (abs(r) > 0.3) {
          cat(sprintf("  Factor%d ~ %s: r=%.3f\n", k, col, r))
        }
      }
    }
  }
}, error = function(e) cat("Correlação:", e$message, "\n"))

cat("\n=== MOFA2 MULTI-VIEW CONCLUÍDO ===\n")
cat("Outputs em:", OUT_DIR, "\n")
cat("Arquivos gerados:\n")
for (f in list.files(OUT_DIR, pattern = "multiview")) cat(" ", f, "\n")
