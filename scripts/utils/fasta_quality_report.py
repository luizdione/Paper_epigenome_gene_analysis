#!/usr/bin/env python3
"""Relatório de qualidade do FASTA do subprojeto."""
import json
import sys
from pathlib import Path

base = Path(".")
sub = sys.argv[1] if len(sys.argv) > 1 else "SUB01"
fa = base / "subprojeto" / sub / f"{sub}_reference_sequences.fa"
js = base / "subprojeto" / sub / f"{sub}_target_genes_coordinates.json"

target = json.loads(js.read_text(encoding="utf-8")).get("genes", {})

cur, buf, data = None, [], {}
def flush():
    global cur, buf, data
    if cur:
        s = "".join(buf)
        n = sum(1 for c in s.upper() if c in "ACGT")
        data[cur.split("_chr")[0]] = (len(s), n, 100*n/max(len(s),1))
    cur, buf = None, []

for line in open(fa):
    line = line.rstrip()
    if line.startswith(">"):
        flush()
        cur = line[1:]
    else:
        buf.append(line)
flush()

print(f"=== FASTA Quality Report — {sub} ===")
print(f"JSON target genes: {len(target)}")
print(f"FASTA sequences:   {len(data)}")
print(f"Match: {'OK' if len(data)==len(target) else 'MISMATCH'}")
print()
print(f"{'Gene':<12}{'bp':>10}{'known':>10}{'%known':>10}")
print("-" * 42)
# Por qualidade crescente para destacar os ruins primeiro
for g in sorted(data, key=lambda x: data[x][2]):
    bp, k, p = data[g]
    flag = "  *all-N*" if p == 0 else ("  low" if p < 1 else "")
    print(f"{g:<12}{bp:>10,}{k:>10,}{p:>9.2f}%{flag}")

# Genes no JSON mas faltando no FASTA
missing = sorted(set(target) - set(data))
if missing:
    print()
    print(f"AUSENTES no FASTA ({len(missing)}): {missing}")
