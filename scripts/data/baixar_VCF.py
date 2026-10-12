import json
import requests
import os
import time

# --- CONFIGURAÇÕES ---
JSON_PATH = "all_target_genes_coordinates.json"
# Diretório base no seu SSD (Windows)
OUTPUT_BASE_DIR = r"E:\data\epigenomas"
VCF_DIR = os.path.join(OUTPUT_BASE_DIR, "vcfs")

# Servidor da Ensembl API
SERVER = "https://rest.ensembl.org"


def criar_diretorio(caminho):
    if not os.path.exists(caminho):
        os.makedirs(caminho)
        print(f"[Info] Diretório criado: {caminho}")


def baixar_vcf_ensembl(gene, chrom, start, end):
    """
    Faz a requisição na API da Ensembl para buscar variantes (SNPs/Indels)
    em uma região genômica específica no formato VCF.
    """
    # Endpoint de sobreposição de variantes para a região no Bos taurus
    endpoint = f"/overlap/region/bos_taurus/{chrom}:{start}-{end}?feature=variation"
    url = f"{SERVER}{endpoint}"

    # Cabeçalho exigindo a resposta em formato VCF
    headers = {"Content-Type": "text/x-vcf"}

    print(f"[{gene}] Baixando VCF da região chr{chrom}:{start}-{end}...")

    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        caminho_arquivo = os.path.join(VCF_DIR, f"{gene}_variants.vcf")
        with open(caminho_arquivo, "w") as f:
            f.write(response.text)
        print(f"  └─ Sucesso! Salvo em: {caminho_arquivo}")
    else:
        print(f"  └─ [Erro] Falha ao baixar {gene}. Status: {response.status_code}")
        print(f"     Detalhes: {response.text}")


def main():
    # 1. Preparar o diretório de saída no SSD
    criar_diretorio(VCF_DIR)

    # 2. Carregar o arquivo JSON
    try:
        with open(JSON_PATH, "r", encoding="utf-8") as f:
            dados = json.load(f)
    except FileNotFoundError:
        print(f"[Erro] Arquivo {JSON_PATH} não encontrado. Coloque-o na mesma pasta do script.")
        return

    genes = dados.get("genes", {})
    total = len(genes)

    print(f"Iniciando o download de VCFs para {total} genes alvo...\n")

    # 3. Iterar sobre cada gene e baixar o VCF da região com flanking
    for i, (nome_gene, info) in enumerate(genes.items(), 1):
        chrom = info.get("chromosome")
        start = info.get("flanking_start")
        end = info.get("flanking_end")

        # O HspB11 no JSON está com chrom = "3" mas verifique casos de valor None
        if not chrom or not start or not end:
            print(f"[{nome_gene}] Coordenadas incompletas. Pulando...")
            continue

        print(f"--- Processando {i}/{total} ---")
        baixar_vcf_ensembl(nome_gene, chrom, start, end)

        # Pausa de 0.5s para respeitar o limite de requisições da API da Ensembl (15 req/seg)
        time.sleep(0.5)

    print("\n[Concluído] Download de VCFs finalizado.")


if __name__ == "__main__":
    main()