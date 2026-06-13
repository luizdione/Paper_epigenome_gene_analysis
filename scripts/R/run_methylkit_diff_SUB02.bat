@echo off
REM =============================================================================
REM SUB02 - Pipeline Fase 1 methylKit em loop pelas 5 comparacoes (Windows)
REM Equivalente ao .sh; assume Docker rstudio rodando com mounts padrao do projeto.
REM =============================================================================

docker exec -i rstudio bash -lc "bash /home/rstudio/projeto/scripts/R/run_methylkit_diff_SUB02.sh"
