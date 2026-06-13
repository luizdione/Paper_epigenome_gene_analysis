import json
import os

# Caminhos
JSON_PATH = "all_target_genes_coordinates.json"
BED_PATH = r"E:\data\epigenomas\targets.bed"


def json_to_bed():
    with open(JSON_PATH, 'r', encoding='utf-8') as f:
        data = json.load(f)

    genes = data.get("genes", {})

    with open(BED_PATH, 'w') as bed:
        for gene_name, info in genes.items():
            chrom = info['chromosome']
            # O BED usa coordenadas 0-based
            start = info['flanking_start']
            end = info['flanking_end']

            # Formato BED: chr  start  end  name
            line = f"{chrom}\t{start}\t{end}\t{gene_name}\n"
            bed.write(line)

    print(f"[Sucesso] Arquivo BED gerado em: {BED_PATH}")


if __name__ == "_