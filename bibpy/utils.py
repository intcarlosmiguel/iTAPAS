import numpy as np
import networkx as nx
from collections import defaultdict
BPR_ALPHA = 0.15
BPR_BETA = 4.0
TOLERANCIA_FLUXO = 1e-10
TOLERANCIA_CUSTO = 1e-10
MAX_ITERACOES = 5000

# =============================================================================
# 2. CARREGAMENTO DE DADOS (chamadas por executar_itapas)
# =============================================================================

def carregar_rede(caminho_arquivo: str) -> nx.DiGraph:
    """Carrega a topologia da rede a partir de um arquivo de texto."""
    grafo = nx.DiGraph()
    dados_rede = np.loadtxt(caminho_arquivo, comments='#')
    for linha in dados_rede:
        u, v, capacidade, tempo_livre = int(linha[0]), int(linha[1]), linha[2], linha[3]
        comprimento = linha[4] if len(linha) >= 5 else 0.0
        grafo.add_edge(
            u, v,
            capacidade=capacidade,
            tempo_fluxo_livre=tempo_livre,
            comprimento=comprimento,
            fluxo=0.0,
            custo=tempo_livre,
            fluxos_por_origem=defaultdict(float)
        )
    print(f"Rede carregada com {grafo.number_of_nodes()} nós e {grafo.number_of_edges()} arcos.")
    return grafo


def carregar_viagens(caminho_arquivo: str) -> dict:
    """Carrega a matriz de viagens (Origem-Destino) a partir de um arquivo."""
    dados_viagens = np.loadtxt(caminho_arquivo, comments='#')
    # Garante que os dados sejam tratados como 2D mesmo com uma única linha
    if dados_viagens.ndim == 1:
        dados_viagens = np.array([dados_viagens])
    viagens = {(int(linha[0]), int(linha[1])): linha[2] for linha in dados_viagens}
    print(f"Matriz OD carregada com {len(viagens)} pares OD.")
    return viagens

def gerar_viagens_aleatorias(grafo: nx.DiGraph, num_pares_od: int, volume_por_par: float) -> dict:
    """
    Gera uma matriz OD com pares aleatórios de nós pertencentes ao grafo.
    
    Args:
        grafo: Grafo direcional da rede topológica.
        num_pares_od: Quantidade de pares Origem-Destino a serem criados.
        volume_por_par: Volume/vazão padronizado atribuído a cada par OD.
        
    Returns:
        Um dicionário mapeando (origem, destino) para seu respectivo volume.
    """
    import random
    
    nos = list(grafo.nodes())
    viagens = {}
    
    # Previne loop infinito caso o grafo seja muito pequeno
    max_possiveis = len(nos) * (len(nos) - 1)
    if num_pares_od > max_possiveis:
        num_pares_od = max_possiveis
        
    while len(viagens) < num_pares_od:
        origem = random.choice(nos)
        destino = random.choice(nos)
        
        if origem != destino and (origem, destino) not in viagens:
            viagens[(origem, destino)] = volume_por_par
            
    print(f"Matriz OD aleatória gerada com {len(viagens)} pares OD.")
    return viagens

def calcular_gap_relativo(grafo: nx.DiGraph, viagens: dict, origens: list) -> float:
    """Calcula o 'Relative Gap', a métrica de convergência padrão."""
    # CORREÇÃO: Usando as chaves corretas
    tempo_total_viagem = sum(d['fluxo'] * d['custo'] for _, _, d in grafo.edges(data=True))
    if abs(tempo_total_viagem) < TOLERANCIA_FLUXO: 
        return 0.0
    
    tempo_viagem_spt = 0
    for origem in origens:
        # CORREÇÃO: Usando a chave 'custo'
        _, custos = nx.dijkstra_predecessor_and_distance(grafo, source=origem, weight='custo')
        for (o, d), demanda in viagens.items():
            if o == origem and d in custos:
                tempo_viagem_spt += demanda * custos[d]
    
    # Prevenção de divisão por zero caso o tempo de viagem spt seja maior
    if tempo_total_viagem <= 0: return float('inf')
    
    return (tempo_total_viagem - tempo_viagem_spt) / tempo_total_viagem