#!/usr/bin/env python3
"""Valida outputs da Fase 3 TSS ±10 kb."""
import csv, os

OUT = "./results/SUB01_GSE161113_run04/fase3_deeptools_10kb"

# Verificar arquivos
files = ["matrix_Angus.csv", "matrix_Nelore.csv",
         "heatmap_Angus.png", "heatmap_Nelore.png",
         "profile_comparison.png", "genes_ordered.csv"]
print("=== Arquivos gerados ===")
for f in files:
    path = os.path.join(OUT, f)
    size = os.path.getsize(path) if os.path.exists(path) else -1
    status = f"{size/1024:.1f} KB" if size > 0 else "AUSENTE ❌"
    print(f"  {f}: {status}")

# Top genes por delta de metilação
print("\n=== Top 10 genes por |delta_pct| (Angus - Nelore) ===")
with open(os.path.join(OUT, "genes_ordered.csv")) as fh:
    rows = list(csv.DictReader(fh))

with_data = [r for r in rows
             if r.get("Angus_mean_pct","0") != "0.0000"
             or r.get("Nelore_mean_pct","0") != "0.0000"]
print(f"Genes com dados: {len(with_data)} / {len(rows)} total")

top = sorted(with_data, key=lambda x: abs(float(x.get("delta_pct", 0))), reverse=True)[:10]
print(f"{'Gene':<20} {'Angus%':>8} {'Nelore%':>8} {'Delta':>8}")
print("-" * 48)
for r in top:
    print(f"{r['gene']:<20} {float(r.get('Angus_mean_pct',0)):>8.2f} "
          f"{float(r.get('Nelore_mean_pct',0)):>8.2f} "
          f"{float(r.get('delta_pct',0)):>+8.2f}")
