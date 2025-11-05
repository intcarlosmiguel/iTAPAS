# solucionador_admm_funcional.py

import networkx as nx
import numpy as np
from collections import defaultdict
import time

# --- Constantes do Modelo ---
# Parâmetros da função de custo BPR: t = t_livre * (1 + ALFA * (fluxo/cap)^BETA)
ALFA = 0.15
BETA = 4.0

# --- Funções Puras e de Utilitário ---

def carregar_rede_de_arquivo(caminho_arquivo):
    """
    Carrega uma rede a partir de um arquivo de texto.

    Args:
        caminho_arquivo (str): O caminho para o arquivo da rede.
                                Formato: nó_origem nó_destino capacidade tempo_livre

    Returns:
        nx.DiGraph: Um objeto de grafo NetworkX com os dados da rede.
    """
    grafo = nx.DiGraph()
    with open(caminho_arquivo, 'r') as f:
        for linha in f:
            if linha.strip().startswith('#') or not linha.strip():
                continue
            partes = linha.split()
            u, v = int(partes[0]), int(partes[1])
            capacidade = float(partes[2])
            tempo_livre = float(partes[3])
            grafo.add_edge(u, v, capacidade=capacidade, tempo_livre=tempo_livre)
    return grafo

def carregar_viagens_de_arquivo(caminho_arquivo):
    """
    Carrega a demanda de viagens (pares O-D) de um arquivo de texto.

    Args:
        caminho_arquivo (str): O caminho para o arquivo de viagens.
                                Formato: origem destino demanda

    Returns:
        list: Uma lista de tuplas (origem, destino, demanda).
    """
    viagens = []
    with open(caminho_arquivo, 'r') as f:
        for linha in f:
            if linha.strip().startswith('#') or not linha.strip():
                continue
            partes = linha.split()
            o, d, demanda = int(partes[0]), int(partes[1]), float(partes[2])
            viagens.append((o, d, demanda))
    return viagens

def calcular_custo_bpr(fluxo, dados_arco):
    """Calcula o custo de um arco usando a função BPR."""
    capacidade = dados_arco['capacidade']
    if capacidade <= 0: return dados_arco['tempo_livre']
    
    return dados_arco['tempo_livre'] * (
        1 + ALFA * (fluxo / capacidade) ** BETA
    )

def calcular_derivada_bpr(fluxo, dados_arco):
    """Calcula a derivada da função de custo BPR."""
    capacidade = dados_arco['capacidade']
    if capacidade <= 0 or fluxo <= 0: return 0.0

    return (dados_arco['tempo_livre'] * ALFA * BETA * 
            (fluxo ** (BETA - 1)) / (capacidade ** BETA))

# --- Funções de Preparação e Estruturação ---

def agrupar_arcos_em_blocos(grafo):
    """
    Implementa um algoritmo guloso de coloração de arestas para agrupar
    os arcos em blocos independentes. Arcos no mesmo bloco não compartilham nós.

    Args:
        grafo (nx.DiGraph): O grafo da rede.

    Returns:
        list: Uma lista de listas, onde cada sublista contém os índices
              dos arcos que pertencem a um bloco.
    """
    cores_arcos = {}
    arcos_incidentes_no = defaultdict(list)
    for u, v in grafo.edges():
        arcos_incidentes_no[u].append((u, v))
        arcos_incidentes_no[v].append((u, v))

    mapa_arco_para_idx = {arco: i for i, arco in enumerate(grafo.edges())}

    for arco in grafo.edges():
        u, v = arco
        cores_usadas = set()
        
        for arco_vizinho in arcos_incidentes_no[u] + arcos_incidentes_no[v]:
            if arco_vizinho in cores_arcos:
                cores_usadas.add(cores_arcos[arco_vizinho])
        
        cor = 1
        while cor in cores_usadas:
            cor += 1
        cores_arcos[arco] = cor

    num_cores = max(cores_arcos.values()) if cores_arcos else 0
    blocos = [[] for _ in range(num_cores)]
    for arco, cor in cores_arcos.items():
        idx_arco = mapa_arco_para_idx[arco]
        blocos[cor - 1].append(idx_arco)
        
    return blocos

def preparar_estruturas_de_dados(grafo, viagens):
    """
    Inicializa todas as estruturas de dados e o estado inicial do solver.

    Args:
        grafo (nx.DiGraph): O grafo da rede.
        viagens (list): A lista de viagens O-D.

    Returns:
        dict: Um dicionário contendo o estado completo do solver.
    """
    # Mapeamentos para usar índices numéricos
    nos = sorted(grafo.nodes())
    mapa_no_para_idx = {no: i for i, no in enumerate(nos)}
    
    arcos = list(grafo.edges(data=True))
    mapa_arco_para_idx = {(u, v): i for i, (u, v, _) in enumerate(arcos)}
    
    origens = sorted(list(set(o for o, d, dem in viagens if dem > 0)))
    mapa_origem_para_idx = {origem: i for i, origem in enumerate(origens)}

    num_nos = len(nos)
    num_arcos = len(arcos)
    num_origens = len(origens)
    
    # Vetor de demanda/geração (g_n^o)
    vetor_g = np.zeros((num_origens, num_nos))
    for o, d, demanda in viagens:
        if demanda > 0 and o in mapa_origem_para_idx:
            idx_o = mapa_origem_para_idx[o]
            idx_d = mapa_no_para_idx[d]
            idx_no_o = mapa_no_para_idx[o]
            vetor_g[idx_o, idx_no_o] += demanda
            vetor_g[idx_o, idx_d] -= demanda
    
    print("2. Agrupando arcos (Edge Coloring)...")
    blocos = agrupar_arcos_em_blocos(grafo)
    print(f"   -> A rede foi dividida em {len(blocos)} blocos.")
    
    # O estado do solver é um dicionário que será passado entre as funções
    estado_solver = {
        "grafo": grafo,
        "arcos": arcos,
        "mapa_no_para_idx": mapa_no_para_idx,
        "mapa_origem_para_idx": mapa_origem_para_idx,
        "num_nos": num_nos,
        "num_arcos": num_arcos,
        "num_origens": num_origens,
        "blocos": blocos,
        "vetor_g": vetor_g,
        
        # Variáveis que mudam a cada iteração
        "fluxos_primarios": np.zeros((num_origens, num_arcos)),
        "variaveis_duais": np.zeros((num_origens, num_nos)),
    }
    return estado_solver

# --- Funções do Núcleo do Algoritmo ADMM ---

def resolver_subproblema_do_arco(idx_arco, estado, fluxos_saida, fluxos_entrada, rho):
    """
    Resolve o subproblema para um único arco usando Projeção de Gradiente.
    Esta função modifica os arrays de fluxo diretamente para eficiência.

    Args:
        idx_arco (int): O índice do arco a ser otimizado.
        estado (dict): O estado atual do solver.
        fluxos_saida (np.array): Fluxo saindo de cada nó, por origem.
        fluxos_entrada (np.array): Fluxo entrando em cada nó, por origem.
        rho (float): O parâmetro de penalidade do ADMM.
    """
    # Desempacota dados necessários do estado
    u, v, dados_arco = estado["arcos"][idx_arco]
    idx_u, idx_v = estado["mapa_no_para_idx"][u], estado["mapa_no_para_idx"][v]
    fluxos_primarios = estado["fluxos_primarios"]

    fluxos_antigos = fluxos_primarios[:, idx_arco].copy()
    fluxo_total_no_arco = np.sum(fluxos_antigos)
    
    # Derivadas da função BPR
    custo_atual = calcular_custo_bpr(fluxo_total_no_arco, dados_arco)
    derivada_custo = calcular_derivada_bpr(fluxo_total_no_arco, dados_arco)
    
    # Derivada de segunda ordem (s_a^o) da Lagrangiana Aumentada (Eq. 48)
    derivada_segunda_ordem = derivada_custo + 2 * rho
    if derivada_segunda_ordem < 1e-9: return

    # Termos constantes (e_n^o) do desbalanceamento de fluxo
    e_cauda = fluxos_saida[:, idx_u] - fluxos_antigos - fluxos_entrada[:, idx_u] - estado["vetor_g"][:, idx_u]
    e_cabeca = fluxos_saida[:, idx_v] - (fluxos_entrada[:, idx_v] - fluxos_antigos) - estado["vetor_g"][:, idx_v]

    # Derivada de primeira ordem (d_a^o) da Lagrangiana Aumentada (Eq. 47)
    termo_lambda = estado["variaveis_duais"][:, idx_u] - estado["variaveis_duais"][:, idx_v]
    termo_rho = rho * (e_cauda - e_cabeca)
    derivada_primeira_ordem = custo_atual + 2 * rho * fluxos_antigos + termo_lambda + termo_rho

    # Atualização do fluxo via Projeção de Gradiente (Eq. 46) com passo alpha=1
    fluxos_novos = fluxos_antigos - (derivada_primeira_ordem / derivada_segunda_ordem)
    fluxos_primarios[:, idx_arco] = np.maximum(0, fluxos_novos)
    
    # Atualiza os fluxos de entrada/saída para que o próximo arco no loop use valores corretos
    mudanca_fluxo = fluxos_primarios[:, idx_arco] - fluxos_antigos
    fluxos_saida[:, idx_u] += mudanca_fluxo
    fluxos_entrada[:, idx_v] += mudanca_fluxo

def atualizar_variaveis_primarias(estado, rho):
    """
    Itera sobre os blocos e atualiza os fluxos primarios (v_a^o).

    Args:
        estado (dict): O estado atual do solver.
        rho (float): O parâmetro de penalidade.

    Returns:
        dict: O estado do solver com os fluxos primários atualizados.
    """
    # Pré-calcula os fluxos de entrada/saída para eficiência
    fluxos_saida = np.zeros((estado["num_origens"], estado["num_nos"]))
    fluxos_entrada = np.zeros((estado["num_origens"], estado["num_nos"]))
    for idx_arco, (u, v, _) in enumerate(estado["arcos"]):
        idx_u, idx_v = estado["mapa_no_para_idx"][u], estado["mapa_no_para_idx"][v]
        fluxo = estado["fluxos_primarios"][:, idx_arco]
        fluxos_saida[:, idx_u] += fluxo
        fluxos_entrada[:, idx_v] += fluxo

    for bloco in estado["blocos"]:
        # Em uma implementação real, este loop interno seria paralelizado
        for idx_arco in bloco:
            resolver_subproblema_do_arco(idx_arco, estado, fluxos_saida, fluxos_entrada, rho)
            
    return estado # O estado foi modificado in-place, mas retornamos por consistência

def atualizar_variaveis_duais(estado, rho):
    """
    Atualiza as variáveis duais (lambda_n^o) com base no resíduo primal.

    Args:
        estado (dict): O estado atual do solver.
        rho (float): O parâmetro de penalidade.

    Returns:
        dict: O estado do solver com as variáveis duais atualizadas.
    """
    fluxos_saida = np.zeros((estado["num_origens"], estado["num_nos"]))
    fluxos_entrada = np.zeros((estado["num_origens"], estado["num_nos"]))
    for idx_arco, (u, v, _) in enumerate(estado["arcos"]):
        idx_u, idx_v = estado["mapa_no_para_idx"][u], estado["mapa_no_para_idx"][v]
        fluxos = estado["fluxos_primarios"][:, idx_arco]
        fluxos_saida[:, idx_u] += fluxos
        fluxos_entrada[:, idx_v] += fluxos
    
    residuo_primal = fluxos_saida - fluxos_entrada - estado["vetor_g"]
    estado["variaveis_duais"] += rho * residuo_primal
    return estado

def verificar_convergencia(estado, viagens):
    """
    Calcula o Relative Gap (RG) para verificar a convergência.

    Args:
        estado (dict): O estado atual do solver.
        viagens (list): A lista de viagens O-D.

    Returns:
        float: O valor do Relative Gap.
    """
    fluxos_totais_arcos = np.sum(estado["fluxos_primarios"], axis=0)
    
    custos_atuais = np.array([
        calcular_custo_bpr(fluxos_totais_arcos[i], estado["arcos"][i][2]) 
        for i in range(estado["num_arcos"])
    ])
    
    nx.set_edge_attributes(estado["grafo"], {
        (u, v): custo for (u, v, _), custo in zip(estado["arcos"], custos_atuais)
    }, name='weight')

    custo_total_sp = 0.0
    for o, d, demanda in viagens:
        if demanda > 0:
            try:
                custo_sp = nx.shortest_path_length(estado["grafo"], source=o, target=d, weight='weight')
                custo_total_sp += custo_sp * demanda
            except nx.NetworkXNoPath:
                return float('inf')
    
    custo_total_sistema = np.sum(fluxos_totais_arcos * custos_atuais)
    if custo_total_sistema < 1e-9: return 0.0
        
    relative_gap = abs(1.0 - custo_total_sp / custo_total_sistema)
    return relative_gap

# --- Função Principal de Execução ---

def resolver_equilibrio_usuario(grafo, viagens, rho, max_iteracoes, tolerancia):
    """
    Executa o loop principal do algoritmo ADMM para encontrar o Equilíbrio de Usuário.

    Args:
        grafo (nx.DiGraph): O grafo da rede.
        viagens (list): A lista de viagens O-D.
        rho (float): Parâmetro de penalidade da Lagrangiana Aumentada.
        max_iteracoes (int): Número máximo de iterações.
        tolerancia (float): Critério de parada para o Relative Gap.

    Returns:
        dict: Um dicionário mapeando cada arco ao seu fluxo de equilíbrio final.
    """
    print("1. Preparando estruturas de dados...")
    estado_solver = preparar_estruturas_de_dados(grafo, viagens)
    
    print("\n3. Iniciando o processo iterativo do ADMM...")
    tempo_inicio = time.time()
    
    for i in range(max_iteracoes):
        estado_solver = atualizar_variaveis_primarias(estado_solver, rho)
        estado_solver = atualizar_variaveis_duais(estado_solver, rho)
        
        relative_gap = verificar_convergencia(estado_solver, viagens)
        
        print(f"   Iteração {i+1:2d}: Relative Gap = {relative_gap:.6e}")
        
        if relative_gap < tolerancia:
            print(f"\nConvergência atingida em {i+1} iterações.")
            break
    else:
        print("\nNúmero máximo de iterações atingido.")

    tempo_fim = time.time()
    print(f"Tempo total de execução: {tempo_fim - tempo_inicio:.4f} segundos.")
    
    # Extrai e formata o resultado final
    fluxos_finais = np.sum(estado_solver["fluxos_primarios"], axis=0)
    return {(u, v): fluxo for (u, v, _), fluxo in zip(estado_solver["arcos"], fluxos_finais)}

# --- Ponto de Entrada Principal ---
if __name__ == "__main__":
    ARQUIVO_REDE = './example/edges.txt'
    ARQUIVO_VIAGENS = './example/od.txt'

    # Carregar dados
    grafo_rede = carregar_rede_de_arquivo(ARQUIVO_REDE)
    lista_viagens = carregar_viagens_de_arquivo(ARQUIVO_VIAGENS)
    
    # Executar o solver
    fluxos_equilibrio = resolver_equilibrio_usuario(
        grafo=grafo_rede, 
        viagens=lista_viagens, 
        rho=0.01, 
        max_iteracoes=500, 
        tolerancia=1e-10
    )
    
    # Exibir resultados
    print("\n--- Fluxos de Equilíbrio nos Arcos ---")
    for (u, v), fluxo in fluxos_equilibrio.items():
        print(f"  Arco ({u:2d} -> {v:2d}): {fluxo:8.2f}")