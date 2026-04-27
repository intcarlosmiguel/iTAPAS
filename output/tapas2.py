import networkx as nx
import numpy as np
from collections import defaultdict
import heapq
import os
import sys

# --- Parâmetros do Modelo ---
ALPHA = 0.15
BETA = 4.0
MAX_ITERACOES = 50       # Aumente se necessário para redes grandes
EPSILON_GAP = 1e-4       # Critério de convergência
PASSO_SHIFT = 0.25       # Passo de suavização (0 < n <= 1). Reduzi para 0.25 para estabilidade.

# --- 1. Funções de Carregamento de Arquivos ---

def carregar_rede(caminho_arquivo: str) -> nx.DiGraph:
    """Carrega a topologia da rede a partir de um arquivo de texto real."""
    print(f"Carregando rede de: {caminho_arquivo}")
    if not os.path.exists(caminho_arquivo):
        print(f"ERRO: Arquivo {caminho_arquivo} não encontrado.")
        sys.exit(1)

    grafo = nx.DiGraph()
    try:
        # Tenta carregar ignorando comentários (#)
        dados_rede = np.loadtxt(caminho_arquivo, comments='#')
        
        # Garante array 2D mesmo se tiver apenas uma linha
        if dados_rede.ndim == 1:
            dados_rede = np.array([dados_rede])
            
        count = 0
        for linha in dados_rede:
            # Assumindo colunas: 0=u, 1=v, 2=capacidade, 3=tempo_livre
            u, v = int(linha[0]), int(linha[1])
            capacidade = float(linha[2])
            tempo_livre = float(linha[3])
            
            grafo.add_edge(
                u, v,
                capacidade=capacidade,
                t0=tempo_livre,
                fluxo=0.0,
                custo=tempo_livre,
                # Estrutura fundamental do TAPAS (Origem -> fluxo)
                fluxos_origem=defaultdict(float)
            )
            count += 1
        print(f"-> Rede carregada com {grafo.number_of_nodes()} nós e {count} arcos.")
        return grafo
    except Exception as e:
        print(f"ERRO CRÍTICO ao ler rede: {e}")
        sys.exit(1)

def carregar_viagens(caminho_arquivo: str) -> dict:
    """Carrega a matriz OD a partir de um arquivo."""
    print(f"Carregando viagens de: {caminho_arquivo}")
    if not os.path.exists(caminho_arquivo):
        print(f"ERRO: Arquivo {caminho_arquivo} não encontrado.")
        sys.exit(1)

    viagens = {}
    try:
        dados_viagens = np.loadtxt(caminho_arquivo, comments='#')
        if dados_viagens.ndim == 1:
            dados_viagens = np.array([dados_viagens])
            
        for linha in dados_viagens:
            o, d, demanda = int(linha[0]), int(linha[1]), float(linha[2])
            if demanda > 0:
                viagens[(o, d)] = demanda
                
        print(f"-> Matriz OD carregada com {len(viagens)} pares O-D.")
        return viagens
    except Exception as e:
        print(f"ERRO CRÍTICO ao ler viagens: {e}")
        sys.exit(1)

# --- 2. Funções de Custo e Caminho Mínimo ---

def custo_bpr(fluxo, capacidade, t0):
    """Função de custo BPR padrão."""
    # Evita divisão por zero caso capacidade seja 0 (proteção)
    if capacidade <= 0: return t0 * 100 
    return t0 * (1.0 + ALPHA * ((fluxo / capacidade) ** BETA))

def atualizar_custos_rede(G):
    """Atualiza o custo de viagem em todos os arcos baseado no fluxo atual."""
    for u, v, dados in G.edges(data=True):
        dados['custo'] = custo_bpr(dados['fluxo'], dados['capacidade'], dados['t0'])

def dijkstra_origem(G, origem):
    """
    Calcula caminhos mínimos de uma origem para todos os nós.
    Retorna: (distancias, predecessores)
    """
    pq = [(0, origem)]
    dist = {n: float('inf') for n in G.nodes}
    dist[origem] = 0
    pred = {n: None for n in G.nodes}
    
    while pq:
        d, u = heapq.heappop(pq)
        
        # Otimização: se já achamos um caminho menor antes, ignora
        if d > dist[u]: continue
        
        for v in G.successors(u):
            custo_arco = G[u][v]['custo']
            nova_dist = d + custo_arco
            
            if nova_dist < dist[v]:
                dist[v] = nova_dist
                pred[v] = u
                heapq.heappush(pq, (nova_dist, v))
    return dist, pred

# --- 3. Métricas de Convergência (GAP) ---

def calcular_gap(G, viagens):
    """Calcula o Average Excess Cost (AEC)."""
    # Custo total atual do sistema (TSTt)
    custo_total_sistema = sum(d['fluxo'] * d['custo'] for u, v, d in G.edges(data=True))
    
    # Custo total se todos viajassem pelo caminho mínimo atual (SPTt)
    custo_minimo_virtual = 0.0
    
    # Agrupa destinos por origem para otimizar chamadas de Dijkstra
    ods_por_origem = defaultdict(list)
    for (o, d), dem in viagens.items():
        ods_por_origem[o].append((d, dem))
    
    total_demanda = 0
    
    for origem, destinos in ods_por_origem.items():
        # Verifica se a origem existe no grafo (pode haver OD para nós que não estão na edge list)
        if origem not in G: continue
        
        dist, _ = dijkstra_origem(G, origem)
        
        for destino, demanda in destinos:
            if destino in dist and dist[destino] != float('inf'):
                custo_minimo_virtual += dist[destino] * demanda
                total_demanda += demanda
            else:
                # Log opcional para debug de conexidade
                pass 
                
    if total_demanda == 0: return 0.0
    
    # AEC = (Custo Atual - Custo Mínimo Possível) / Total Viagens
    return 1 - (custo_minimo_virtual / custo_total_sistema)

# --- 4. Algoritmo Principal ---

def alocacao_inicial_aon(G, viagens):
    """Realiza a alocação Tudo-ou-Nada inicial."""
    print("Iniciando alocação Tudo-ou-Nada...")
    
    # Limpa fluxos anteriores
    for u, v, d in G.edges(data=True):
        d['fluxo'] = 0.0
        d['fluxos_origem'] = defaultdict(float) # Reinicia dicionário
        
    ods_por_origem = defaultdict(list)
    for (o, d), dem in viagens.items():
        ods_por_origem[o].append((d, dem))
        
    for origem, destinos in ods_por_origem.items():
        if origem not in G: continue
        
        _, pred = dijkstra_origem(G, origem)
        
        for destino, demanda in destinos:
            if destino not in pred: continue
            
            # Reconstrói caminho (backtracking do destino até origem)
            curr = destino
            caminho_valido = True
            edges_path = []
            
            while curr != origem:
                pai = pred[curr]
                if pai is None: 
                    caminho_valido = False
                    break
                edges_path.append((pai, curr))
                curr = pai
            
            if caminho_valido:
                for u, v in edges_path:
                    G[u][v]['fluxo'] += demanda
                    G[u][v]['fluxos_origem'][origem] += demanda
    
    atualizar_custos_rede(G)

def iteracao_tapas_simplificada(G, viagens):
    """
    Realiza o shift de fluxo baseado em origens.
    Retorna True se houve mudança significativa nos fluxos.
    """
    fluxo_mudou = False
    origens_ativas = set(o for o, d in viagens.keys())
    
    # Itera sobre cada origem (Conceito Origin-Based)
    for origem in origens_ativas:
        if origem not in G: continue
        
        # 1. Calcula a árvore ótima atual (T*)
        _, pred_otimo = dijkstra_origem(G, origem)
        
        # 2. Verifica violações de fluxo nos arcos
        # Para cada nó v, verifica quem chega nele
        for v in G.nodes():
            if v == origem: continue
            
            u_otimo = pred_otimo[v] # De onde eu deveria vir se fosse ótimo?
            if u_otimo is None: continue 
            
            # Verifica todos os arcos chegando em v
            for u in list(G.predecessors(v)):
                # Se o arco (u,v) NÃO é o arco ótimo para esta origem...
                if u != u_otimo:
                    dados_arco_ruim = G[u][v]
                    fluxo_origem_aqui = dados_arco_ruim['fluxos_origem'].get(origem, 0.0)
                    
                    # ...mas tem fluxo passando por ele
                    if fluxo_origem_aqui > 1e-6:
                        # Identificamos um PAS local:
                        # Rota ruim: chega em v via u
                        # Rota boa:  chega em v via u_otimo
                        
                        # Calcula Delta de fluxo
                        delta = fluxo_origem_aqui * PASSO_SHIFT
                        
                        # Executa o Shift
                        # 1. Tira do ruim
                        dados_arco_ruim['fluxos_origem'][origem] -= delta
                        dados_arco_ruim['fluxo'] -= delta
                        
                        # 2. Põe no bom
                        if G.has_edge(u_otimo, v):
                            dados_arco_bom = G[u_otimo][v]
                            dados_arco_bom['fluxos_origem'][origem] += delta
                            dados_arco_bom['fluxo'] += delta
                            fluxo_mudou = True
                            
    return fluxo_mudou

# --- Execução Main ---

def main():
    # Definição dos caminhos dos arquivos
    arquivo_rede = './fortaleza/edges.txt'
    arquivo_viagens = './od_outputs/OD_10_300/OD_1.txt'
    
    print("--- INICIANDO ALGORITMO DE ALOCAÇÃO DE TRÁFEGO ---")
    
    # 1. Carregar Dados
    G = carregar_rede(arquivo_rede)
    viagens = carregar_viagens(arquivo_viagens)
    
    if G.number_of_nodes() == 0 or len(viagens) == 0:
        print("Interrompendo: Rede ou Viagens vazias.")
        return

    # 2. Solução Inicial
    alocacao_inicial_aon(G, viagens)
    gap = calcular_gap(G, viagens)
    print(f"\n[Inicial] Alocação Tudo-ou-Nada completa.")
    print(f"Gap Inicial: {gap:.8f}\n")
    
    # 3. Loop de Equilíbrio
    print(f"Iniciando iterações (Max: {MAX_ITERACOES})...")
    
    for i in range(1, MAX_ITERACOES + 1):
        mudou = iteracao_tapas_simplificada(G, viagens)
        atualizar_custos_rede(G)
        
        gap = calcular_gap(G, viagens)
        
        print(f"Iter {i:03d}: Gap Relativo (AEC) = {gap:.8f}")
        
        if gap < EPSILON_GAP:
            print(f"\nCONVERGÊNCIA ATINGIDA na iteração {i} (Gap < {EPSILON_GAP})")
            break
        
        if not mudou and i > 1:
            print(f"\nEstagnação de fluxos na iteração {i}. Terminando.")
            break
            
    print("\n--- Processo Finalizado ---")

if __name__ == "__main__":
    main()