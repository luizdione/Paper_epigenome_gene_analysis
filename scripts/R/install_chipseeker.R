if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
cat("=== Instalando ChIPseeker e dependencias bovinas ===\n")
BiocManager::install(
  c("ChIPseeker", "TxDb.Btaurus.UCSC.bosTau9.refGene", "org.Bt.eg.db"),
  ask = FALSE,
  update = FALSE
)
cat("=== Verificando instalacao ===\n")
pkgs <- c("ChIPseeker", "TxDb.Btaurus.UCSC.bosTau9.refGene", "org.Bt.eg.db")
for (pkg in pkgs) {
  ok <- requireNamespace(pkg, quietly=TRUE)
  cat(sprintf("  %s: %s\n", pkg, if(ok) "OK" else "FALHOU"))
}
