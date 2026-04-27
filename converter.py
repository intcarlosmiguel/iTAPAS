import sys
import os
import numpy as np

def converter_para_tntp(input_file, output_file):
    print(f"Lendo arquivo: {input_file}...")

    try:
        # Lê o arquivo inteiro para uma matriz NumPy
        # Colunas Entrada: 0=Init, 1=Term, 2=Cap, 3=FFT, 4=Len
        raw_data = np.loadtxt(input_file)
        
        # Verifica se o arquivo não está vazio
        if raw_data.size == 0:
            print("Aviso: Arquivo vazio ignorado.")
            return

        # Se houver apenas uma linha, o numpy cria um array 1D, precisamos garantir que seja 2D
        if raw_data.ndim == 1:
            raw_data = raw_data.reshape(1, -1)

    except Exception as e:
        print(f"Erro ao ler arquivo: {e}")
        return

    # --- Estatísticas para o Metadata ---
    # Pega o maior valor entre a coluna 0 (Init) e 1 (Term) convertendo para inteiro
    num_nodes = int(np.max(raw_data[:, 0:2])) 
    num_zones = num_nodes # Assumindo Zonas == Nós
    num_links = raw_data.shape[0] # Número de linhas

    # --- Construção da Matriz de Saída ---
    
    # Cria vetores constantes para as colunas novas
    ones = np.ones(num_links)
    zeros = np.zeros(num_links)
    
    # Valores padrão
    col_b     = ones * 0.15
    col_power = ones * 4
    col_speed = zeros
    col_toll  = zeros
    col_type  = ones

    # Monta a nova matriz empilhando as colunas na ordem correta do TNTP:
    # Ordem TNTP: Init, Term, Cap, Length, FFT, B, Power, Speed, Toll, Type
    
    output_data = np.column_stack((
        raw_data[:, 0], # Init
        raw_data[:, 1], # Term
        raw_data[:, 2], # Capacity
        raw_data[:, 4], # Length (Era a coluna 4 na entrada)
        raw_data[:, 3], # FFT    (Era a coluna 3 na entrada) -> TROCA FEITA AQUI
        col_b,
        col_power,
        col_speed,
        col_toll,
        col_type
    ))

    print(f"Escrevendo arquivo: {output_file}...")
    print(f"Estatísticas: {num_nodes} nós, {num_links} links.")

    # Garante que a pasta de destino existe
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w') as f:
        # 1. Metadados
        f.write(f"<NUMBER OF ZONES> {num_zones}\n")
        f.write(f"<NUMBER OF NODES> {num_nodes}\n")
        f.write(f"<FIRST THRU NODE> 1\n")
        f.write(f"<NUMBER OF LINKS> {num_links}\n")
        f.write(f"<END OF METADATA>\n\n\n")

        # 2. Cabeçalho das colunas
        f.write("~ \tInit node\tTerm node\tCapacity\tLength\tFree Flow Time\tB\tPower\tSpeed limit\tToll\tType\t;\n")

        # 3. Dados
        # Formatamos usando np.savetxt para velocidade.
        # fmt define o tipo: %d (int), %.6g (float genérico), e adicionamos o ';' no final
        np.savetxt(f, output_data, 
                   fmt='\t%d\t%d\t%.6g\t%.6g\t%.6g\t%.2f\t%d\t%d\t%d\t%d\t;', 
                   delimiter='\t')

    print("Conversão concluída com sucesso!")

# --- Configuração de execução ---

input_dir = './fortaleza/'
output_dir = './output/edges_tntp/'

# Garante que o diretório de entrada existe antes de listar
if os.path.exists(input_dir):
    files = os.listdir(input_dir)
    files = [f for f in files if f.startswith('edges_')]
    for file in files:
        if file.endswith('.txt'): # Filtro simples de segurança
            input_path = os.path.join(input_dir, file)
            
            # Gera nome de saída
            suffix = file.split("_")[-1].replace(".txt","")
            output_name = f'edges_{suffix}_tntp.txt'
            output_path = os.path.join(output_dir, output_name)

            converter_para_tntp(input_path, output_path)
else:
    print(f"Diretório não encontrado: {input_dir}")