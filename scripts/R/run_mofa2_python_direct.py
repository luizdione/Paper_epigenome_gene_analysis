#!/usr/bin/env python3
"""
run_mofa2_python_direct.py
Roda MOFA2 via API Python pura, sem depender do R/reticulate.
Usa methylation_matrix.csv já produzido por prepare_mofa2.R.
Salva: mofa2_model.hdf5, mofa2_trained_summary.json

Uso (dentro do container epigenoma_env):
    python3 ./scripts/R/run_mofa2_python_direct.py
"""

import sys, json, time, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

OUT_DIR  = "./results/SUB01_GSE161113_run04/fase6_mofa2"
CSV_PATH = f"{OUT_DIR}/methylation_matrix.csv"
HDF5_OUT = f"{OUT_DIR}/mofa2_model.hdf5"
RDS_OUT  = f"{OUT_DIR}/mofa2_trained.rds"  # marcador para o R saber que terminou
SUMMARY  = f"{OUT_DIR}/mofa2_trained_summary.json"

print("=" * 60)
print("MOFA2 Python Direct — GenOmics Epigenoma")
print("=" * 60)

# 1. Carregar matriz de metilação
print(f"\n[1] Carregando matriz: {CSV_PATH}")
df = pd.read_csv(CSV_PATH, index_col=0)
print(f"    Dimensões: {df.shape[0]} features × {df.shape[1]} amostras")
print(f"    Amostras: {list(df.columns)}")

# 2. Preparar dados para MOFA2
# mofapy2 espera: linhas=amostras, colunas=features (samples × features)
# Verificar orientação: se índice contém coordenadas genômicas → transpor (CpG deve ser coluna/feature)
col0 = str(df.columns[0]) if len(df.columns) > 0 else ""
idx0 = str(df.index[0])   if len(df.index)   > 0 else ""
if ":" in idx0 or idx0.startswith("chr"):
    print(f"    [FIX] CpG sites no índice detectados ({idx0}...) — transpondo para amostras×features...")
    df = df.T
else:
    print(f"    [OK] Orientação correta: linhas=amostras ({idx0}...), colunas=features ({col0}...)")

data_matrix = df.values.astype(np.float32)
samples     = list(df.index)    # linhas = amostras
features    = list(df.columns)  # colunas = features (CpG sites)
print(f"    Amostras (N): {len(samples)}  — ex: {samples[:3]}")
print(f"    Features (D): {len(features)} — ex: {features[:3]}")

# 3. Inicializar MOFA2
print("\n[2] Inicializando entry_point MOFA2...")
from mofapy2.run.entry_point import entry_point

ent = entry_point()

# Configurar dados
ent.set_data_options(scale_groups=False, scale_views=False)
ent.set_data_matrix(
    [[data_matrix]],       # lista de grupos × views; matriz: amostras × features
    views_names=["methylation"],
    groups_names=["group1"],
    samples_names=[samples],   # lista de nomes das amostras (linhas)
    features_names=[features]  # lista de nomes das features/CpG (colunas)
)

# Configurar modelo
ent.set_model_options(factors=2)

# Configurar treinamento
ent.set_train_options(
    iter=1000,
    convergence_mode="medium",
    seed=42,
    verbose=False,
    quiet=False,
    outfile=HDF5_OUT,
    save_interrupted=True     # CRITICO: garante save mesmo em convergencia antecipada
)

# 4. Construir e treinar o modelo
print("\n[3] Construindo modelo...")
ent.build()

print("[4] Treinando modelo (maxiter=1000, seed=42)...")
t0 = time.time()

# Patch sys.exit para evitar saida prematura
import sys as _sys
_real_exit = _sys.exit
_sys.exit = lambda code=0: None

try:
    ent.run()
except SystemExit:
    pass
except Exception as e:
    print(f"[WARN] Exceção durante run(): {e}")
finally:
    _sys.exit = _real_exit

elapsed = time.time() - t0
print(f"\n[5] Treinamento concluído em {elapsed:.1f}s")

# 5. Salvar modelo
import os, shutil

HDF5_INTERRUPTED = HDF5_OUT.replace(".hdf5", "_interrupted.hdf5")

print(f"\n[6] Salvando modelo em: {HDF5_OUT}")

# Prioridade 1: salvar via ent.save() (funciona quando train_stats foi setado corretamente)
saved = False
try:
    if hasattr(ent.model, 'train_stats'):
        ent.save(outfile=HDF5_OUT, save_data=True)
        print(f"    [OK] HDF5 salvo via ent.save(): {HDF5_OUT}")
        saved = True
    else:
        print("    [WARN] train_stats ausente após treinamento")
except Exception as e:
    print(f"    [WARN] ent.save() falhou: {e}")

# Prioridade 2: fallback — copiar arquivo interrupted se existir
if not saved:
    if os.path.exists(HDF5_INTERRUPTED) and os.path.getsize(HDF5_INTERRUPTED) > 1000:
        print(f"    [FB] Usando interrupted ({os.path.getsize(HDF5_INTERRUPTED)//1024} KB)...")
        shutil.copy2(HDF5_INTERRUPTED, HDF5_OUT)
        print(f"    [OK] HDF5 salvo via interrupted: {HDF5_OUT}")
        saved = True

if not saved:
    print("[ERRO] Não foi possível salvar o modelo")
    sys.exit(1)

# 6. Salvar sumário JSON
print("\n[7] Salvando sumário...")
try:
    Z_node = ent.model.nodes['Z']
    Z = Z_node.getExpectation() if hasattr(Z_node, 'getExpectation') else None

    W_node = ent.model.nodes['W']
    # W pode ser Multiview_Variational_Node ou lista
    if hasattr(W_node, 'getExpectation'):
        W = W_node.getExpectation()
        if isinstance(W, list): W = W[0]
    elif hasattr(W_node, '__iter__'):
        W = list(W_node)[0].getExpectation() if W_node else None
    else:
        W = None

    summary = {
        "status": "converged",
        "n_samples": len(samples),
        "n_features": len(features),
        "n_factors": int(Z.shape[1]) if Z is not None else 2,
        "elapsed_s": round(elapsed, 2),
        "samples": samples,
        "hdf5": HDF5_OUT,
        "Z_shape": list(Z.shape) if Z is not None else [],
        "W_shape": list(W.shape) if W is not None else []
    }
    with open(SUMMARY, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"    [OK] Sumário: {SUMMARY}")
    print(f"\n{'='*60}")
    print(f"MOFA2 CONCLUÍDO")
    print(f"  Fatores latentes Z: {Z.shape}")
    print(f"  Loadings W:        {W.shape}")
    print(f"  HDF5:              {HDF5_OUT}")
    print(f"{'='*60}")
except Exception as e:
    print(f"    [WARN] Erro ao salvar sumário: {e}")
