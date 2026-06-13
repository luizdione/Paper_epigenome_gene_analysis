#!/usr/bin/env Rscript
# =============================================================================
# Fase 6 — Preparação de Dados para Integração Multi-ômica (MOFA2)
# =============================================================================
# Prepara as matrizes de entrada para MOFA2 a partir de:
#   1. Dados de metilação (CSV de DMCs ou matriz de metilação por amostra)
#   2. Dados genômicos (VCF filtrado do 1000 Bull Genomes)
#
# MOFA2 recebe uma lista de matrizes (views) × amostras (groups) e decompõe
# a variância em fatores latentes que explicam co-variação entre ômicas.
#
# Saídas:
#   - mofa2_input.rds      — objeto pronto para MOFA2::run_mofa()
#   - mofa2_input_summary.txt — resumo das dimensões/views
#   - methylation_matrix.csv — matriz amostras × CpGs (para auditoria)
#   - genotype_matrix.csv   — matriz amostras × variantes (se VCF fornecido)
#
# Uso:
#   Rscript prepare_mofa2.R \
#       --methylation results/fase1_methylkit/DMCs_all.csv \
#       --metadata samples.csv \
#       --vcf /data/1000bulls/filtered_target_regions.vcf.gz \
#       --json all_target_genes_coordinates.json \
#       --out results/fase6_mofa2
#
# Dependências:
#   MOFA2 (Bioconductor), VariantAnnotation ou vcfR (para VCF), jsonlite
#
# Autor: Luiz Dione Barbosa De-Melo — GenOmics Epigenoma 2026
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) stop("Instale 'optparse'.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Instale 'jsonlite'.")
  library(optparse)
  library(jsonlite)
})

# ------------------------------------------------------------- CLI
opt <- parse_args(OptionParser(option_list = list(
  make_option("--methylation", type = "character", default = NULL,
              help = "CSV de metilação (DMCs_all.csv ou matriz amostras × CpGs)"),
  make_option("--metadata", type = "character", default = NULL,
              help = "samples.csv com sampleID, breed, treatment, file"),
  make_option("--vcf", type = "character", default = NULL,
              help = "VCF filtrado (opcional — 1000 Bull Genomes)"),
  make_option("--json", type = "character", default = "all_target_genes_coordinates.json",
              help = "JSON dos genes-alvo"),
  make_option("--min-var", type = "numeric", default = 0.01,
              help = "Variância mínima por feature para inclusão (filtro de ruído)"),
  make_option("--data-dir", type = "character", default = "/data/epigenomas/NCBI_GSE161113",
              help = "Pasta com arquivos .gz de cobertura Bismark por amostra"),
  make_option("--out", type = "character", default = "results/fase6_mofa2")
)))

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(opt$out, "prepare_mofa2.log")
log_con <- file(log_path, open = "wt")
log <- function(...) { msg <- paste0(...); message(msg); writeLines(msg, log_con) }

log("[Info] Preparação MOFA2 iniciada em: ", Sys.time())

# ------------------------------------------------------------- Metadata
if (is.null(opt$metadata)) stop("--metadata é obrigatório.")
meta <- read.csv(opt$metadata, stringsAsFactors = FALSE)
log("[Info] Amostras: ", nrow(meta))
log("[Info] Breeds: ", paste(unique(meta$breed), collapse = ", "))

# Definir grupos MOFA (por breed ou breed+condition)
if ("breed" %in% colnames(meta)) {
  meta$mofa_group <- meta$breed
} else {
  meta$mofa_group <- "all"
}

# ------------------------------------------------------------- View 1: Metilação
views <- list()
if (!is.null(opt$methylation) && file.exists(opt$methylation)) {
  log("[Info] Carregando dados de metilação: ", opt$methylation)
  meth <- read.csv(opt$methylation, stringsAsFactors = FALSE)

  # Detectar formato: methylKit output (chr, start, end, meth.diff...) ou matriz
  if (nrow(meth) == 0) {
    log("[WARN] Arquivo de metilação vazio (0 DMCs) — MOFA2 requer ao menos 1 feature. Encerrando.")
    close(log_con)
    quit(save = "no", status = 0)
  }
  if (all(c("chr", "start") %in% colnames(meth))) {
    log("[Info] Formato detectado: methylKit DMC output")

    # ---------------------------------------------------------------
    # Estratégia: construir matriz amostras × CpGs-DMC a partir dos
    # arquivos Bismark Coverage individuais de cada amostra.
    # Cada feature = CpG significativo; valor = % metilação daquela amostra.
    # Isso produz uma matriz real (n_amostras × n_CpGs) para o MOFA2.
    # ---------------------------------------------------------------
    meth_matrix <- NULL

    if (!is.null(opt$metadata) && file.exists(opt$metadata)) {
      meta_local <- read.csv(opt$metadata, stringsAsFactors = FALSE)

      # Localizar pasta com os arquivos filtrados
      data_dir <- if (!is.null(opt$`data-dir`)) opt$`data-dir` else ""
      has_data_dir <- nchar(data_dir) > 0 && dir.exists(data_dir)

      if (has_data_dir && "file" %in% colnames(meta_local) && "sampleID" %in% colnames(meta_local)) {
        # Definir features = CpGs significativos como chave "chr:start"
        cpg_keys <- paste0(meth$chr, ":", meth$start)

        sample_rows <- list()
        for (i in seq_len(nrow(meta_local))) {
          sid  <- meta_local$sampleID[i]
          fname <- meta_local$file[i]
          fpath <- file.path(data_dir, fname)

          # Tentar .gz e sem .gz
          if (!file.exists(fpath) && file.exists(sub("\\.gz$", "", fpath))) {
            fpath <- sub("\\.gz$", "", fpath)
          }
          if (!file.exists(fpath)) {
            log("[WARN] Arquivo não encontrado para amostra ", sid, ": ", fpath)
            next
          }

          tryCatch({
            # Bismark coverage: chr, start, end, methylation%, count_M, count_U
            raw <- read.table(gzfile(fpath), sep="\t", header=FALSE, stringsAsFactors=FALSE,
                              col.names=c("chr","start","end","pct","count_m","count_u"))
            raw$key <- paste0(raw$chr, ":", raw$start)
            # Manter somente CpGs que são DMCs significativos
            raw_sig <- raw[raw$key %in% cpg_keys, ]
            vals <- setNames(raw_sig$pct, raw_sig$key)
            # Garantir ordem e preencher NAs onde amostra não tem cobertura
            row_vals <- vals[cpg_keys]
            names(row_vals) <- cpg_keys
            sample_rows[[sid]] <- row_vals
            log("[Info] Amostra ", sid, ": ", length(vals), " CpGs cobertos de ", length(cpg_keys), " DMCs")
          }, error = function(e) {
            log("[WARN] Erro ao ler ", sid, ": ", conditionMessage(e))
          })
        }

        if (length(sample_rows) >= 2) {
          meth_matrix <- do.call(rbind, lapply(sample_rows, as.numeric))
          rownames(meth_matrix) <- names(sample_rows)
          colnames(meth_matrix) <- cpg_keys
          # Remover CpGs com > 50% NAs
          na_frac <- colMeans(is.na(meth_matrix))
          meth_matrix <- meth_matrix[, na_frac < 0.5, drop = FALSE]
          # Imputar NAs restantes com mediana da coluna
          for (j in seq_len(ncol(meth_matrix))) {
            col_med <- median(meth_matrix[, j], na.rm = TRUE)
            meth_matrix[is.na(meth_matrix[, j]), j] <- col_med
          }
          log("[Info] Matriz de metilação por amostra: ",
              nrow(meth_matrix), " amostras × ", ncol(meth_matrix), " CpGs")
        } else {
          log("[WARN] Menos de 2 amostras carregadas dos arquivos .gz — fallback para bins genômicos")
        }
      }
    }

    # Fallback se não foi possível construir matriz por amostra
    if (is.null(meth_matrix)) {
      # Tentar recuperar a partir do methylation_matrix.csv já existente no --out
      existing_csv <- file.path(opt$out, "methylation_matrix.csv")
      if (file.exists(existing_csv)) {
        log("[Info] Fallback: reutilizando methylation_matrix.csv existente em: ", existing_csv)
        existing <- read.csv(existing_csv, row.names = 1, stringsAsFactors = FALSE)
        if (nrow(existing) >= 2) {
          meth_matrix <- as.matrix(existing)
          log("[Info] Matriz reutilizada: ", nrow(meth_matrix), " amostras × ", ncol(meth_matrix), " CpGs")
        } else {
          log("[WARN] methylation_matrix.csv existente tem menos de 2 amostras — ignorando")
        }
      }
      if (is.null(meth_matrix)) {
        log("[ERRO] Não foi possível construir matriz por amostra e nenhum fallback disponível.")
        log("       Verifique se --data-dir (", opt$`data-dir`, ") é acessível e contém os arquivos FILTERED_*.gz")
        log("[Info] Usando fallback de último recurso: bins genômicos (1Mb) — MOFA2 terá apenas 1 amostra agregada")
        meth$bin <- paste0(sub("^chr", "", meth$chr), "_", floor(meth$start / 1e6), "Mb")
        bin_meth <- tapply(meth$meth.diff, meth$bin, mean, na.rm = TRUE)
        meth_matrix <- matrix(bin_meth, nrow = 1)
        colnames(meth_matrix) <- names(bin_meth)
        rownames(meth_matrix) <- "aggregated"
        log("[WARN] Matriz de fallback (bins 1Mb): ", ncol(meth_matrix), " bins — NÃO ADEQUADA para MOFA2 real")
      }
    }

    write.csv(meth_matrix, file.path(opt$out, "methylation_matrix.csv"))
    views[["methylation"]] <- meth_matrix

  } else {
    log("[Info] Formato detectado: matriz amostras × CpGs")
    meth_matrix <- as.matrix(meth[, -1])  # primeira coluna = sampleID
    rownames(meth_matrix) <- meth[, 1]
    # Filtrar features com baixa variância
    feat_var <- apply(meth_matrix, 2, var, na.rm = TRUE)
    keep <- feat_var >= opt$`min-var`
    meth_matrix <- meth_matrix[, keep, drop = FALSE]
    log("[Info] Features mantidas após filtro de variância: ",
        ncol(meth_matrix), " / ", length(feat_var))
    write.csv(meth_matrix, file.path(opt$out, "methylation_matrix.csv"))
    views[["methylation"]] <- meth_matrix
  }
} else {
  log("[WARN] Dados de metilação não fornecidos ou não encontrados.")
}

# ------------------------------------------------------------- View 2: Genótipos (VCF)
if (!is.null(opt$vcf) && file.exists(opt$vcf)) {
  log("[Info] Carregando VCF: ", opt$vcf)

  if (requireNamespace("vcfR", quietly = TRUE)) {
    library(vcfR)
    vcf <- read.vcfR(opt$vcf, verbose = FALSE)
    gt <- extract.gt(vcf, element = "GT")

    # Converter genótipos para numérico: 0/0=0, 0/1=1, 1/1=2
    gt_num <- apply(gt, c(1, 2), function(x) {
      if (is.na(x) || x == "./.") return(NA_real_)
      alleles <- strsplit(x, "[/|]")[[1]]
      sum(as.integer(alleles))
    })

    # Transpor: amostras × variantes
    gt_matrix <- t(gt_num)
    log("[Info] Matriz de genótipos: ", nrow(gt_matrix), " amostras × ", ncol(gt_matrix), " variantes")

    # Filtrar variantes monomórficas
    var_maf <- apply(gt_matrix, 2, function(col) {
      af <- mean(col, na.rm = TRUE) / 2
      min(af, 1 - af)
    })
    keep <- var_maf >= 0.05
    gt_matrix <- gt_matrix[, keep, drop = FALSE]
    log("[Info] Variantes após filtro MAF >= 0.05: ", ncol(gt_matrix))

    write.csv(gt_matrix, file.path(opt$out, "genotype_matrix.csv"))
    views[["genotype"]] <- gt_matrix

  } else if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
    library(VariantAnnotation)
    vcf <- readVcf(opt$vcf)
    gt <- geno(vcf)$GT

    gt_num <- apply(gt, c(1, 2), function(x) {
      if (is.na(x) || x == "./.") return(NA_real_)
      alleles <- strsplit(x, "[/|]")[[1]]
      sum(as.integer(alleles))
    })
    gt_matrix <- t(gt_num)
    write.csv(gt_matrix, file.path(opt$out, "genotype_matrix.csv"))
    views[["genotype"]] <- gt_matrix
    log("[Info] VCF processado via VariantAnnotation.")

  } else {
    log("[WARN] Nem vcfR nem VariantAnnotation estão instalados. VCF não processado.")
    log("       Instale: BiocManager::install(c('vcfR', 'VariantAnnotation'))")
  }
} else {
  log("[Info] VCF não fornecido — MOFA2 rodará somente com view de metilação.")
}

# ------------------------------------------------------------- Montar input MOFA2
if (length(views) == 0) {
  log("[ERRO] Nenhuma view carregada. Verifique os dados de entrada.")
  close(log_con)
  stop("Nenhuma view disponível para MOFA2.")
}

log("[Info] Views disponíveis: ", paste(names(views), collapse = ", "))

if (requireNamespace("MOFA2", quietly = TRUE)) {
  library(MOFA2)
  log("[Info] MOFA2 disponível — criando objeto de entrada.")

  # Alinhar amostras entre views
  common_samples <- Reduce(intersect, lapply(views, rownames))
  if (length(common_samples) == 0) {
    log("[WARN] Sem amostras em comum entre views. Usando todas as amostras de cada view.")
    # MOFA2 aceita amostras parciais (missing views)
  } else {
    log("[Info] Amostras em comum entre views: ", length(common_samples))
  }

  # Criar MOFA object
  mofa_data <- list()
  for (v in names(views)) {
    mofa_data[[v]] <- t(views[[v]])  # MOFA2 espera features × amostras
  }

  # Definir grupos baseados em breed
  group_labels <- NULL
  if (nrow(meta) > 0 && "mofa_group" %in% colnames(meta)) {
    # mapear sampleID → grupo
    sample_groups <- setNames(meta$mofa_group, meta$sampleID)
    # Para cada view, determinar grupo de cada amostra
    first_view <- views[[1]]
    matched <- sample_groups[rownames(first_view)]
    if (!all(is.na(matched))) {
      group_labels <- as.character(matched)
      group_labels[is.na(group_labels)] <- "unknown"
    }
  }

  # Salvar como RDS para uso posterior
  mofa_input <- list(
    data = mofa_data,
    groups = group_labels,
    sample_metadata = meta,
    views = names(views),
    created = Sys.time()
  )
  saveRDS(mofa_input, file.path(opt$out, "mofa2_input.rds"))
  log("[OK] Input MOFA2 salvo em: ", file.path(opt$out, "mofa2_input.rds"))

  # Gerar script auxiliar para rodar MOFA2
  mofa_script <- file.path(opt$out, "run_mofa2.R")
  out_dir_str  <- opt$out
  rds_path_str <- file.path(opt$out, "mofa2_input.rds")
  hdf5_str     <- file.path(opt$out, "mofa2_model.hdf5")
  writeLines(c(
    "#!/usr/bin/env Rscript",
    "# Script para rodar MOFA2 — Projeto Epigenoma GenOmics",
    paste0("# Gerado por prepare_mofa2.R | ", Sys.time()),
    "",
    "# CRITICO: setar PYTHONPATH e Python ANTES de qualquer import do reticulate/MOFA2",
    "Sys.setenv(PYTHONPATH = paste(",
    "  '/usr/local/lib/python3.10/dist-packages',",
    "  '/usr/lib/python3/dist-packages',",
    "  sep = ':'",
    "))",
    "Sys.setenv(RETICULATE_PYTHON = '/usr/bin/python3')",
    "reticulate::use_python('/usr/bin/python3', required = TRUE)",
    "",
    "library(MOFA2)",
    "",
    paste0("OUT_DIR  <- '", out_dir_str, "'"),
    paste0("RDS_PATH <- '", rds_path_str, "'"),
    "",
    "# Carregar input preparado",
    "input <- readRDS(RDS_PATH)",
    "",
    "# Validar número de amostras antes de treinar",
    "n_samples <- ncol(input$data[[1]])",
    "message('[Info] Amostras no input MOFA2: ', n_samples)",
    "if (n_samples < 2) {",
    "  stop('[ERRO] MOFA2 requer >= 2 amostras. Input tem apenas ', n_samples,",
    "       '. Rode prepare_mofa2.R novamente com --data-dir acessivel.')",
    "}",
    "message('[Info] Views: ', paste(names(input$data), collapse = ', '))",
    "message('[Info] Dimensoes view 1: ',",
    "        nrow(input$data[[1]]), ' features x ', ncol(input$data[[1]]), ' amostras')",
    "",
    "# Criar MOFAobject",
    "mofa_obj <- create_mofa(input$data)",
    "",
    "# Configurar opcoes",
    "data_opts  <- get_default_data_options(mofa_obj)",
    "model_opts <- get_default_model_options(mofa_obj)",
    "model_opts$num_factors <- 2",
    "train_opts <- get_default_training_options(mofa_obj)",
    "train_opts$maxiter <- 1000",
    "train_opts$convergence_mode <- 'slow'",
    "train_opts$seed <- 42",
    "",
    "mofa_obj <- prepare_mofa(mofa_obj,",
    "  data_options  = data_opts,",
    "  model_options = model_opts,",
    "  training_options = train_opts",
    ")",
    "",
    "# Treinar modelo",
    paste0("mofa_obj <- run_mofa(mofa_obj, outfile = '", hdf5_str, "', use_basilisk = FALSE)"),
    "",
    "# Resultados — Variancia explicada",
    paste0("pdf(file.path(OUT_DIR, 'mofa2_variance_explained.pdf'), width = 8, height = 6)"),
    "plot_variance_explained(mofa_obj)",
    "dev.off()",
    "",
    "# Resultados — Fatores latentes",
    "n_factors <- min(model_opts$num_factors, get_dimensions(mofa_obj)[['K']])",
    paste0("pdf(file.path(OUT_DIR, 'mofa2_factors.pdf'), width = 10, height = 6)"),
    "tryCatch(plot_factors(mofa_obj, factors = seq_len(n_factors)),",
    "         error = function(e) message('[WARN] plot_factors: ', conditionMessage(e)))",
    "dev.off()",
    "",
    "# Resultados — Pesos (loadings) da primeira view, fator 1",
    "first_view <- names(input$data)[1]",
    paste0("pdf(file.path(OUT_DIR, 'mofa2_weights.pdf'), width = 10, height = 6)"),
    "tryCatch(plot_weights(mofa_obj, view = first_view, factor = 1, nfeatures = 20),",
    "         error = function(e) message('[WARN] plot_weights: ', conditionMessage(e)))",
    "dev.off()",
    "",
    "# Salvar objeto treinado como RDS para analise downstream",
    "saveRDS(mofa_obj, file.path(OUT_DIR, 'mofa2_model.rds'))",
    "",
    "message('[OK] MOFA2 concluido. Resultados em: ', OUT_DIR)"
  ), mofa_script)
  log("[OK] Script de execução MOFA2 gerado: ", mofa_script)

} else {
  log("[WARN] MOFA2 não instalado. Salvando dados preparados para uso posterior.")
  log("       Instale: BiocManager::install('MOFA2')")
  saveRDS(views, file.path(opt$out, "mofa2_input.rds"))
}

# ------------------------------------------------------------- Resumo
summary_path <- file.path(opt$out, "mofa2_input_summary.txt")
sink(summary_path)
cat("GenOmics Epigenoma — Resumo da preparação MOFA2\n")
cat("=", rep("=", 49), "\n", sep = "")
cat("Data:", format(Sys.time()), "\n")
cat("Amostras (metadata):", nrow(meta), "\n")
cat("Views preparadas:", length(views), "\n")
for (v in names(views)) {
  m <- views[[v]]
  cat(sprintf("  %s: %d amostras × %d features\n", v, nrow(m), ncol(m)))
}
cat("\nArquivos gerados:\n")
for (f in list.files(opt$out, full.names = FALSE)) {
  cat(sprintf("  - %s\n", f))
}
sink()

log("[FIM] Preparação MOFA2 concluída. Resumo em: ", summary_path)
close(log_con)
cat("[OK] Resultados em:", opt$out, "\n")
