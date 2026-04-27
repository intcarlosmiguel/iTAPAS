import sys
import os
import numpy as np
from collections import defaultdict

def converter_od_para_tntp(input_file, output_file, max_node_id=None):
    """
    Converte arquivo OD simples para formato TNTP.
    :param max_node_id: Se fornecido, força o número de zonas. 
                        Se None, usa o maior ID encontrado no arquivo.
    """
    print(f"Lendo arquivo OD: {input_file}...")

    try:
        # Lê o arquivo. Espera-se 3 colunas: Origin, Destination, Flow
        raw_data = np.loadtxt(input_file)

        # Verifica se o arquivo está vazio
        if raw_data.size == 0:
            print("Aviso: Arquivo vazio ignorado.")
            return

        # Garante que seja 2D mesmo se tiver apenas uma linha
        if raw_data.ndim == 1:
            raw_data = raw_data.reshape(1, -1)

    except Exception as e:
        print(f"Erro ao ler arquivo: {e}")
        return

    # --- Processamento dos Dados ---

    # Converte colunas para tipos apropriados
    origins = raw_data[:, 0].astype(int)
    dests   = raw_data[:, 1].astype(int)
    flows   = raw_data[:, 2] # Mantém float para o fluxo

    # --- Lógica do Número de Zonas ---
    max_id_encontrado = int(np.max(raw_data[:, 0:2]))

    if max_node_id is None:
        num_zones = max_id_encontrado
    else:
        num_zones = int(max_node_id)
        # Aviso de segurança caso o arquivo tenha IDs maiores que o limite forçado
        if max_id_encontrado > num_zones:
            print(f"AVISO: O arquivo contém o nó {max_id_encontrado}, "
                  f"mas você forçou o limite de {num_zones} zonas. "
                  "Alguns dados podem ser ignorados ou ficar fora do intervalo.")

    total_flow = np.sum(flows)

    # Agrupa os dados por Origem usando um dicionário para acesso rápido
    od_map = defaultdict(list)
    for o, d, flow in zip(origins, dests, flows):
        od_map[o].append((d, flow))

    # --- Escrita do Arquivo ---
    
    print(f"Escrevendo arquivo: {output_file}...")
    print(f"Estatísticas: {num_zones} zonas (Max ID nos dados: {max_id_encontrado}), Fluxo Total: {total_flow:.2f}")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w') as f:
        # 1. Metadados
        f.write(f"<NUMBER OF ZONES> {num_zones}\n")
        f.write(f"<TOTAL OD FLOW> {total_flow}\n")
        f.write(f"<END OF METADATA>\n\n")

        # 2. Dados por Origem
        # Iteramos de 1 até num_zones definido
        for zone_id in range(1, num_zones + 1):
            f.write(f"Origin {zone_id}\n")
            
            if zone_id in od_map:
                destinations = od_map[zone_id]
                
                count = 0
                for dest, flow in destinations:
                    # Se forçamos um limite de zonas e o destino é maior que esse limite,
                    # teoricamente não deveríamos escrever, mas o TNTP aceita se for apenas passagem.
                    # Mantemos a escrita normal.
                    
                    f.write(f"\t{dest} : {flow:.6g} ;")
                    
                    count += 1
                    if count % 5 == 0:
                        f.write("\n")
                
                if count % 5 != 0:
                    f.write("\n")
            else:
                f.write("\n")

    print("Conversão de OD concluída com sucesso!")

# --- Configuração de execução ---

input_dir = './od_outputs/OD_10_300/'
output_dir = './input/od_tntp/'

# SE VOCÊ QUISER FORÇAR UM NÚMERO DE ZONES, COLOQUE AQUI:
# Exemplo: FORCED_ZONES = 1072
FORCED_ZONES = None 

if os.path.exists(input_dir):
    files = os.listdir(input_dir)
    for file in files:
        if file.endswith('.txt'):
            input_path = os.path.join(input_dir, file)
            
            suffix = file.split("_")[-1].replace(".txt","")
            if suffix == file.replace(".txt", ""):
                 output_name = f'{file.replace(".txt", "")}_tntp.txt'
            else:
                 output_name = f'od_{suffix}_tntp.txt'
            
            output_path = os.path.join(output_dir, output_name)

            # Passando o argumento opcional aqui
            converter_od_para_tntp(input_path, output_path, max_node_id=36296)
else:
    print(f"Diretório de entrada não encontrado: {input_dir}")