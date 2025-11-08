"""
Implementação do algoritmo iTAPAS em estilo funcional para o Problema de
Alocação de Tráfego (TAP).

Este script evita o uso de classes, preferindo funções que recebem o estado
da rede (grafo, viagens, etc.) como argumento e retornam um novo estado
modificado. Isso torna o fluxo de dados explícito e o código mais modular.

ORGANIZAÇÃO: Funções ordenadas pela ordem de execução no algoritmo.
"""

import networkx as nx # type: ignore
import numpy as np
from collections import defaultdict


# =============================================================================
# PARÂMETROS GLOBAIS DO ALGORITMO
# =============================================================================

BPR_ALPHA = 0.15
BPR_BETA = 4.0
TOLERANCIA_FLUXO = 1e-9
TOLERANCIA_CUSTO = 1e-9
MAX_ITERACOES = 100


# =============================================================================
# 1. FUNÇÃO PRINCIPAL DE EXECUÇÃO (PONTO DE ENTRADA)
# =============================================================================

def executar_itapas(arquivo_rede: str, arquivo_viagens: str, max_iter: int, gap_convergencia: float):
    """Função principal que orquestra a execução do algoritmo iTAPAS."""
    grafo = carregar_rede(arquivo_rede)
    viagens = carregar_viagens(arquivo_viagens)
    origens = []
    for (o,d) in viagens.keys():
        if o not in origens:
            origens.append(o)
    #origens = [34219,26251,3140]
    conjunto_pas = []
    grafo = atribuicao_inicial(grafo, viagens)

    for i in range(MAX_ITERACOES):
        for origem in origens:
            grafo, conjunto_pas = processar_origem(grafo, origem, conjunto_pas)
        grafo, conjunto_pas = deslocamento_global_pas(grafo, conjunto_pas)
        gap = calcular_gap_relativo(grafo, viagens, origens)
        print(f"Gap Relativo: {gap:.8e} | Iteração: {i+1} | PAS ativos: {len(conjunto_pas)}")
        if gap < gap_convergencia:
            print("Convergência atingida!")
            break

    print("\nAlocação finalizada.")
    return grafo


# =============================================================================
# 2. CARREGAMENTO DE DADOS (chamadas por executar_itapas)
# =============================================================================

def carregar_rede(caminho_arquivo: str) -> nx.DiGraph:
    """Carrega a topologia da rede a partir de um arquivo de texto."""
    grafo = nx.DiGraph()
    dados_rede = np.loadtxt(caminho_arquivo, comments='#')
    for linha in dados_rede:
        u, v, capacidade, tempo_livre = int(linha[0]), int(linha[1]), linha[2], linha[3]
        grafo.add_edge(
            u, v,
            capacidade=capacidade,
            tempo_fluxo_livre=tempo_livre,
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


# =============================================================================
# 3. ALOCAÇÃO INICIAL (chamada por executar_itapas)
# =============================================================================

def atribuicao_inicial(grafo: nx.DiGraph, viagens: dict) -> nx.DiGraph:
    """Realiza a alocação inicial 'tudo-ou-nada'."""
    print("Realizando alocação inicial 'tudo-ou-nada'...")
    for (o, d), demanda in viagens.items():
        # CORREÇÃO: Usa a chave 'tempo_fluxo_livre' para o peso
        caminho = nx.shortest_path(grafo, source=o, target=d, weight='tempo_fluxo_livre')
        for i in range(len(caminho) - 1):
            u, v = caminho[i], caminho[i+1]
            grafo[u][v]['fluxo'] += demanda
            grafo[u][v]['fluxos_por_origem'][o] += demanda
    
    atualizar_todos_custos(grafo)
    return grafo


def atualizar_todos_custos(grafo: nx.DiGraph):
    """Atualiza os custos de todos os arcos no grafo (modifica in-place)."""
    for u, v in grafo.edges():
        atualizar_custo_arco(grafo, u, v)


def atualizar_custo_arco(grafo: nx.DiGraph, u: int, v: int):
    """Atualiza o custo de um único arco no grafo (modifica o grafo in-place)."""
    dados_arco = grafo[u][v]
    fluxo = dados_arco['fluxo']
    capacidade = dados_arco['capacidade']
    tempo_livre = dados_arco['tempo_fluxo_livre']
    
    custo = tempo_livre * (1 + BPR_ALPHA * (fluxo / capacidade) ** BPR_BETA)
    grafo[u][v]['custo'] = custo


# =============================================================================
# 4. PROCESSAMENTO POR ORIGEM (chamada no loop principal)
# =============================================================================

def processar_origem(grafo: nx.DiGraph, origem: int, conjunto_pas: list) -> (nx.DiGraph, list): # type: ignore
    """Executa uma iteração de equilíbrio para uma única origem."""
    # CORREÇÃO: Usa a chave 'custo' para o peso do Dijkstra
    preds, custos_spt = nx.dijkstra_predecessor_and_distance(grafo, source=origem, weight='custo')
    arcos_desequilibrados,ids = identificar_arcos_desequilibrados(grafo, origem, custos_spt)
    #print(origem)
    #print(arcos_desequilibrados)
    # Ordena arcos_desequilibrados pela primeira coluna e depois pela segunda.
    # Reordena também `ids` para manter o paralelismo entre as listas.
    if arcos_desequilibrados:
        # Parear arcos com seus ids e ordenar: ascendente por u, descendente por v
        pareados = list(zip(arcos_desequilibrados, ids))
        pareados.sort(key=lambda item: (item[0][0], -item[0][1]))
        arcos_desequilibrados, ids = zip(*pareados)
        arcos_desequilibrados = list(arcos_desequilibrados)
        ids = list(ids)
    while arcos_desequilibrados:
        u, v = arcos_desequilibrados.pop(0)
        pas = identificar_pas_fluxo_maximo(grafo, u, v, origem, preds)
        if not pas: 
            continue
        
        fluxo_deslocado, grafo = deslocar_fluxo_no_pas(grafo, pas)
        if fluxo_deslocado > TOLERANCIA_FLUXO:
            conjunto_pas,entrou = adicionar_pas_ao_conjunto(conjunto_pas, pas)
    return grafo, conjunto_pas


def identificar_arcos_desequilibrados(grafo: nx.DiGraph, origem: int, custos_spt: dict) -> list:
    """Encontra arcos com fluxo positivo da origem que não estão na SPT."""
    desequilibrados = []
    ids_desequilibrados = []
    i = 0
    for u, v, dados in grafo.edges(data=True):
        if dados['fluxos_por_origem'][origem] > TOLERANCIA_FLUXO:
            if u not in custos_spt or v not in custos_spt: 
                continue
            # CORREÇÃO: Acessando a chave 'custo' em português.
            custo_reduzido = dados['custo'] + custos_spt[u] - custos_spt[v]
            if custo_reduzido > TOLERANCIA_CUSTO:
                desequilibrados.append((u, v))
                ids_desequilibrados.append(i)
        i += 1
    return desequilibrados, ids_desequilibrados


def identificar_pas_fluxo_maximo(grafo: nx.DiGraph, u_deseq: int, v_deseq: int, origem: int, preds: dict) -> dict:
    """Identifica um PAS usando a busca retroativa de fluxo máximo."""
    cabeca = v_deseq
    
    nos_caminho_s1 = [cabeca]
    no_atual = cabeca
    while no_atual in preds and origem not in nos_caminho_s1:
        no_anterior = preds[no_atual][0]
        nos_caminho_s1.insert(0, no_anterior)
        no_atual = no_anterior
    
    conjunto_nos_s1 = set(nos_caminho_s1)
    caminho_retroativo_s2 = []
    no_atual = u_deseq
    nos_visitados_s2 = {no_atual}
    
    while no_atual not in conjunto_nos_s1:
        predecessores = list(grafo.predecessors(no_atual))
        if not predecessores: 
            return None

        melhor_pred = max(
            predecessores, 
            key=lambda p: grafo[p][no_atual]['fluxos_por_origem'][origem],
            default=None
        )
        
        if melhor_pred is None or grafo[melhor_pred][no_atual]['fluxos_por_origem'][origem] < TOLERANCIA_FLUXO:
            return None
        if melhor_pred in nos_visitados_s2:
            return None

        caminho_retroativo_s2.insert(0, (melhor_pred, no_atual))
        nos_visitados_s2.add(melhor_pred)
        no_atual = melhor_pred

    cauda = no_atual
    
    idx_inicio_s1 = nos_caminho_s1.index(cauda)
    s1 = [(nos_caminho_s1[i], nos_caminho_s1[i+1]) for i in range(idx_inicio_s1, len(nos_caminho_s1)-1)]
    s2 = caminho_retroativo_s2 + [(u_deseq, v_deseq)]
    
    return {'s1': s1, 's2': s2, 'origem': origem, 'cabeca': cabeca, 'cauda': cauda}


def deslocar_fluxo_no_pas(grafo: nx.DiGraph, pas: dict) -> (float, nx.DiGraph): # type: ignore
    """Calcula e aplica o deslocamento de fluxo ótimo em um PAS."""
    s1, s2, origem = pas['s1'], pas['s2'], pas['origem']
    
    # CORREÇÃO: Acessando a chave 'custo'
    custo_s1 = sum(grafo[u][v]['custo'] for u, v in s1)
    custo_s2 = sum(grafo[u][v]['custo'] for u, v in s2)
    if custo_s1 > custo_s2:
        s1, s2, custo_s1, custo_s2 = s2, s1, custo_s2, custo_s1

    if (custo_s2 - custo_s1) < TOLERANCIA_CUSTO:
        return 0.0, grafo

    def derivada_bpr(u, v):
        dados = grafo[u][v]
        if dados['fluxo'] < TOLERANCIA_FLUXO: 
            return 0.0
        # CORREÇÃO: Usando as chaves corretas
        termo = (BPR_ALPHA * BPR_BETA * dados['tempo_fluxo_livre'] / 
                 (dados['capacidade'] ** BPR_BETA))
        return termo * (dados['fluxo'] ** (BPR_BETA - 1))

    derivada_s1 = sum(derivada_bpr(u, v) for u, v in s1)
    derivada_s2 = sum(derivada_bpr(u, v) for u, v in s2)
    
    denominador = derivada_s1 + derivada_s2
    if denominador < TOLERANCIA_FLUXO: 
        return 0.0, grafo

    fluxo_maximo_deslocavel = min(grafo[u][v]['fluxos_por_origem'][origem] for u, v in s2)
    delta_newton = (custo_s2 - custo_s1) / denominador
    delta = max(0, min(delta_newton, fluxo_maximo_deslocavel))

    if delta < TOLERANCIA_FLUXO: 
        return 0.0, grafo

    for u, v in s1:
        grafo[u][v]['fluxo'] += delta
        grafo[u][v]['fluxos_por_origem'][origem] += delta
        atualizar_custo_arco(grafo, u, v)
    for u, v in s2:
        grafo[u][v]['fluxo'] -= delta
        grafo[u][v]['fluxos_por_origem'][origem] -= delta
        atualizar_custo_arco(grafo, u, v)
        
    return delta, grafo


def adicionar_pas_ao_conjunto(conjunto_pas: list, novo_pas: dict) ->  (list,bool): # type: ignore
    """Adiciona um novo PAS ao conjunto global, evitando duplicatas topológicas."""
    sig_s1 = tuple(sorted(novo_pas['s1']))
    sig_s2 = tuple(sorted(novo_pas['s2']))
    nova_assinatura = tuple(sorted((sig_s1, sig_s2)))
    for pas_existente in conjunto_pas:
        ex_sig_s1 = tuple(sorted(pas_existente['s1']))
        ex_sig_s2 = tuple(sorted(pas_existente['s2']))
        assinatura_existente = tuple(sorted((ex_sig_s1, ex_sig_s2)))
        if nova_assinatura == assinatura_existente:
            return conjunto_pas,False
    
    return conjunto_pas + [novo_pas],True


# =============================================================================
# 5. DESLOCAMENTO GLOBAL DE PAS (chamada no loop principal)
# =============================================================================

def deslocamento_global_pas(grafo: nx.DiGraph, conjunto_pas: list, num_deslocamentos=20) -> (nx.DiGraph, list): # type: ignore
    """Reequilibra repetidamente todos os PAS no conjunto global."""
    if not conjunto_pas: 
        return grafo, conjunto_pas
    
    pas_ativos = list(conjunto_pas)
    for _ in range(num_deslocamentos):
        pas_a_remover = []
        for pas in pas_ativos:
            origem = pas['origem']
            fluxo_min_s2 = min(grafo[u][v]['fluxos_por_origem'][origem] for u, v in pas['s2'])
            
            if fluxo_min_s2 < TOLERANCIA_FLUXO:
                pas_a_remover.append(pas)
            else:
                dx, grafo = deslocar_fluxo_no_pas(grafo, pas)
        
        if pas_a_remover:
            pas_ativos = [p for p in pas_ativos if p not in pas_a_remover]

    return grafo, pas_ativos


# =============================================================================
# 6. CÁLCULO DE CONVERGÊNCIA (chamada no loop principal)
# =============================================================================

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


# =============================================================================
# EXECUÇÃO
# =============================================================================

if __name__ == '__main__':

    grafo_final = executar_itapas(
        arquivo_rede='./fortaleza/edges.txt',
        arquivo_viagens='./od_outputs/OD_10_300/OD_0.txt',
        max_iter=100,
        gap_convergencia=1e-10
    )

    # Imprime os resultados finais
    """ print("\n--- Resultados Finais ---")
    print("Fluxo e Custo por Arco:")
    for u, v, dados in sorted(grafo_final.edges(data=True)):
        print(f"Arco ({u}->{v}): Fluxo = {dados['fluxo']:.2f}, Custo = {dados['custo']:.2f}") """