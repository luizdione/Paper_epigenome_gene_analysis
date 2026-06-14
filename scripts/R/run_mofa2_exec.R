#!/usr/bin/env Rscript
# Executa MOFA2 com os dados ja preparados em fase6_mofa2/
# Roda dentro do container epigenoma_env

suppressPackageStartupMessages({
  library(MOFA2)
})

out_dir <- "./results/SUB01_GSE161113_run04/fase6_mofa2"
rds_path <- file.path(out_dir, "mofa2_input.rds")

cat("=== MOFA2 Execucao ===\n")
cat("Input RDS:", rds_path, "\n")

if (!file.exists(rds_path)) {
  stop("mofa2_input.rds nao encontrado: ", rds_path)
}

input <- readRDS(rds_path)
cat("Views disponíveis:", paste(names(input$data), collapse = ", "), "\n")
for (v in names(input$data)) {
  m <- input$data[[v]]
  cat(sprintf("  %s: %d features x %d amostras\n", v, nrow(m), ncol(m)))
}

# Criar MOFAobject
mofa_obj <- create_mofa(input$data)
cat("\nMOFA object criado\n")

# Configurar opcoes
data_opts  <- get_default_data_options(mofa_obj)
model_opts <- get_default_model_options(mofa_obj)
model_opts$num_factors <- min(10, ncol(input$data[[1]]) - 1)
train_opts <- get_default_training_options(mofa_obj)
train_opts$maxiter          <- 1000
train_opts$convergence_mode <- "medium"
train_opts$seed             <- 42
train_opts$verbose          <- TRUE

mofa_obj <- prepare_mofa(mofa_obj,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)

# Treinar modelo
hdf5_path <- file.path(out_dir, "mofa2_model.hdf5")
cat("\nTreinando MOFA2 (maxiter=1000)...\n")
# Usar Python do sistema (mofapy2 instalado via pip3) em vez de basilisk/conda
reticulate::use_python("/usr/bin/python3", required = TRUE)
mofa_obj <- run_mofa(mofa_obj, outfile = hdf5_path, use_basilisk = FALSE)

cat("\n=== RESULTADOS ===\n")
print(mofa_obj)

# Variancia explicada por fator
var_exp <- calculate_variance_explained(mofa_obj)
cat("\nVariancia explicada por view/fator:\n")
print(var_exp$r2_per_factor)

# Salvar objeto treinado
trained_rds <- file.path(out_dir, "mofa2_trained.rds")
saveRDS(mofa_obj, trained_rds)
cat("\nModelo salvo:", trained_rds, "\n")

# Graficos
pdf(file.path(out_dir, "mofa2_variance_explained.pdf"), width = 8, height = 6)
tryCatch(plot_variance_explained(mofa_obj), error = function(e) cat("Plot variancia:", e$message, "\n"))
dev.off()

pdf(file.path(out_dir, "mofa2_factors.pdf"), width = 10, height = 6)
tryCatch(plot_factors(mofa_obj, factors = 1:min(5, model_opts$num_factors)),
         error = function(e) cat("Plot fatores:", e$message, "\n"))
dev.off()

pdf(file.path(out_dir, "mofa2_weights_meth.pdf"), width = 10, height = 6)
tryCatch(plot_weights(mofa_obj, view = "methylation", factor = 1, nfeatures = 15),
         error = function(e) cat("Plot weights:", e$message, "\n"))
dev.off()

# Top weights do fator 1
cat("\nTop features - Fator 1 (metilacao):\n")
tryCatch({
  w <- get_weights(mofa_obj, views = "methylation", factors = 1, as.data.frame = TRUE)
  w <- w[order(abs(w$value), decreasing = TRUE), ]
  print(head(w, 15))
  write.csv(w, file.path(out_dir, "mofa2_weights_fator1.csv"), row.names = FALSE)
  cat("Weights salvos:", file.path(out_dir, "mofa2_weights_fator1.csv"), "\n")
}, error = function(e) cat("Erro weights:", e$message, "\n"))

cat("\n=== MOFA2 CONCLUIDO ===\n")
cat("Outputs em:", out_dir, "\n")
