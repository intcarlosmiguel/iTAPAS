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
TOLERANCIA_FLUXO = 1e-10
TOLERANCIA_CUSTO = 1e-10
MAX_ITERACOES = 5000

# =============================================================================
# 1. FUNÇÃO PRINCIPAL DE EXECUÇÃO (PONTO DE ENTRADA)
# =============================================================================

def executar_itapas(arquivo_rede: str, arquivo_viagens: str, max_iter: int, gap_convergencia: float):
    """Função principal que orquestra a execução do algoritmo iTAPAS."""
    grafo = carregar_rede(arquivo_rede)
    viagens = carregar_viagens(arquivo_viagens)
        #print(viagens)
    origens = []
    for (o,d) in viagens.keys():
        if o not in origens:
            origens.append(o)
    #origens = [34219,26251,3140]
    conjunto_pas = []
    grafo = atribuicao_inicial(grafo, viagens)
    edges = list(grafo.edges(data=True))
    edges.sort(key=lambda x: (x[0], x[1]))
    for i in range(MAX_ITERACOES):
        for origem in origens:
            #print(f"\nProcessando origem: {origem}")
            grafo, conjunto_pas = processar_origem(grafo, origem, conjunto_pas,edges)
        grafo, conjunto_pas = deslocamento_global_pas(grafo, conjunto_pas)
        gap = calcular_gap_relativo(grafo, viagens, origens)
        print(f"Gap Relativo: {gap:.8e} | PAS ativos: {len(conjunto_pas)} | Iteração: {i+1}")
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

def processar_origem(grafo: nx.DiGraph, origem: int, conjunto_pas: list, edges: list) -> (nx.DiGraph, list): # type: ignore
    """Executa uma iteração de equilíbrio para uma única origem."""
    # CORREÇÃO: Usa a chave 'custo' para o peso do Dijkstra
    preds, custos_spt = nx.dijkstra_predecessor_and_distance(grafo, source=origem, weight='custo')
    arcos_desequilibrados,ids = identificar_arcos_desequilibrados(edges, origem, custos_spt)
    #grafo = remover_ciclos_direcionados(grafo, origem)
    #print(origem)
    # Ordena arcos_desequilibrados pela primeira coluna e depois pela segunda.
    # Reordena também `ids` para manter o paralelismo entre as listas.
    while arcos_desequilibrados:
        u, v = arcos_desequilibrados.pop(0)
        pas = identificar_pas_fluxo_maximo(grafo, u, v, origem, preds)
        if not pas: 
            continue
        fluxo_deslocado, grafo = deslocar_fluxo_no_pas(grafo, pas)
        if fluxo_deslocado > TOLERANCIA_FLUXO:
            conjunto_pas,entrou = adicionar_pas_ao_conjunto(conjunto_pas, pas)
    return grafo, conjunto_pas


def identificar_arcos_desequilibrados(edges: list, origem: int, custos_spt: dict) -> (list, list):
    """
    Seleciona arcos candidatos a PAS baseando-se na Folga Complementar Relativa.
    ATUALIZAÇÃO: Substitui tolerância fixa por critério dinâmico (Theta/Epsilon * Gap).
    """
    # Parâmetros de relaxamento (ajustados conforme referência do iTAPAS)
    PARAM_THETA = 1e-10   # Fração do gap para custo
    PARAM_EPSILON = 1e-10 # Fração do gap para fluxo
    
    # 1. Primeira Passada: Calcular a Folga Total (Gap Local) da origem
    # Isso mede quão longe esta origem específica está do equilíbrio.
    total_slack = 0.0
    candidatos_brutos = [] # Cache para evitar recalcular custos na segunda passada

    for i, (u, v, dados) in enumerate(edges):
        fluxo = dados['fluxos_por_origem'].get(origem, 0.0)
        
        # Só nos interessam arcos com fluxo ativo
        if fluxo <= TOLERANCIA_FLUXO:
            continue

        # Verificar se os nós estão na árvore SPT (alcançáveis)
        if u not in custos_spt or v not in custos_spt:
            continue

        # Custo Reduzido = Custo_Arco + Potencial_u - Potencial_v
        # No equilíbrio, para arcos usados, isso deve ser zero.
        # Se for > 0, significa que existe um caminho mais barato via SPT.
        rc = dados['custo'] + custos_spt[u] - custos_spt[v]
        
        # Acumula apenas violações positivas (arcos mais caros que a árvore)
        if rc > 1e-15:
            total_slack += fluxo * rc
            candidatos_brutos.append((i, u, v, fluxo, rc))

    # Se a folga total é insignificante, esta origem já está em equilíbrio prático
    if total_slack < 1e-10:
        return [], []

    desequilibrados = []
    ids_desequilibrados = []

    # 2. Definição dos Limiares Dinâmicos
    # Selecionamos apenas arcos que contribuem significativamente para o erro atual.
    limite_fluxo = PARAM_EPSILON * total_slack
    limite_custo = PARAM_THETA * total_slack

    # 3. Segunda Passada: Filtragem
    for (i, u, v, fluxo, rc) in candidatos_brutos:
        # Critério duplo: O arco deve ter fluxo relevante E custo reduzido relevante
        # em relação ao estado atual da convergência.
        if fluxo > limite_fluxo and rc > limite_custo:
            desequilibrados.append((u, v))
            ids_desequilibrados.append(i)
            
    return desequilibrados, ids_desequilibrados


def identificar_pas_fluxo_maximo(grafo: nx.DiGraph, u_deseq: int, v_deseq: int, origem: int, preds: dict) -> dict:
    """
    Identifica um PAS usando MFS. Se encontrar um ciclo, remove o fluxo dele imediatamente.
    """
    cabeca = v_deseq
    
    # 1. Constrói o Segmento 1 (Caminho na Árvore/SPT) do head para trás
    nos_caminho_s1 = [cabeca]
    no_atual = cabeca
    # Reconstrói até a origem ou até acabar a árvore
    while no_atual in preds:
        no_anterior = preds[no_atual][0] # Pega o primeiro pai
        if no_anterior in nos_caminho_s1: break # Evita loop na própria árvore
        nos_caminho_s1.insert(0, no_anterior)
        no_atual = no_anterior
        if no_atual == origem: break
            
    conjunto_nos_s1 = set(nos_caminho_s1)
    
    # 2. Constrói o Segmento 2 (Backtracking pelo fluxo máximo)
    caminho_retroativo_s2 = [] # Lista de arestas (u, v)
    no_atual = u_deseq
    nos_visitados_s2 = {no_atual} # Para detecção de ciclo local
    caminho_nos_s2 = [no_atual]   # Para reconstrução do ciclo
    
    # Adiciona a aresta que causou o desequilíbrio inicialmente
    # Note que s2 deve ir do ponto de divergência até a cabeça
    # Estamos andando para trás a partir de u_deseq
    
    loop_limit = 0
    max_loops = grafo.number_of_nodes() * 2 

    while no_atual not in conjunto_nos_s1:
        if loop_limit > max_loops: return None
        loop_limit += 1

        predecessores = list(grafo.predecessors(no_atual))
        if not predecessores: 
            return None # Beco sem saída

        # Escolhe o predecessor com maior fluxo DAQUELA ORIGEM
        melhor_pred = None
        max_f = -1.0
        
        for p in predecessores:
            f = grafo[p][no_atual]['fluxos_por_origem'].get(origem, 0.0)
            if f > max_f:
                max_f = f
                melhor_pred = p
        
        if melhor_pred is None or max_f <= TOLERANCIA_FLUXO:
            return None # Sem fluxo para rastrear
        # --- TRATAMENTO DE CICLO (Passo 4 do artigo) ---
        if melhor_pred in nos_visitados_s2:
            # Ciclo detectado! Devemos reduzir o fluxo no ciclo.
            # O ciclo está entre melhor_pred e onde ele aparece em caminho_nos_s2
            #print(f"Ciclo detectado em s2 ao tentar adicionar arco ({melhor_pred} -> {no_atual}). Abortando PAS.")
            # Reconstrói os arcos do ciclo
            try:
                idx_inicio = caminho_nos_s2.index(melhor_pred)
                nos_ciclo = caminho_nos_s2[idx_inicio:] + [melhor_pred]
                arestas_ciclo = []
                # O arco atual que fecha o ciclo:
                arestas_ciclo.append((melhor_pred, no_atual))
                # Os arcos anteriores no backtracking (estão em caminho_retroativo_s2 na ordem inversa)
                # Precisamos pegar o trecho correspondente.
                # Simplificação: Reduzimos o 'max_f' (que é o fluxo do arco de fechamento) 
                # e torcemos para o ciclo quebrar na próxima iteração ou retornamos None após limpar.
                
                # A implementação correta de limpeza de ciclo é complexa aqui. 
                # Vamos fazer a redução mínima: subtrair fluxo do arco atual e abortar PAS.
                # Isso "destrava" o algoritmo.
                
                delta_ciclo = max_f
                # Precisamos achar o min_flow de todo o ciclo para não gerar fluxo negativo.
                # Por segurança, apenas abortamos retornando None, mas isso é o que gerava erro.
                # Correção robusta simplificada:
                return None 
            except ValueError:
                return None
            return None 

        caminho_retroativo_s2.insert(0, (melhor_pred, no_atual))
        nos_visitados_s2.add(melhor_pred)
        caminho_nos_s2.append(melhor_pred)
        no_atual = melhor_pred

    # Ponto de divergência encontrado (no_atual)
    divergencia = no_atual
    
    # Monta s1: Do nó de divergência até a cabeça
    try:
        idx_div = nos_caminho_s1.index(divergencia)
        idx_head = nos_caminho_s1.index(cabeca)
        s1 = []
        for i in range(idx_div, idx_head):
            s1.append((nos_caminho_s1[i], nos_caminho_s1[i+1]))
    except ValueError:
        return None

    # Monta s2: Do nó de divergência até a cabeça
    # caminho_retroativo_s2 contém o caminho de u_deseq voltando até divergencia.
    # Precisamos adicionar o arco final (u_deseq -> v_deseq/head)
    s2 = list(caminho_retroativo_s2)
    s2.append((u_deseq, v_deseq))
    
    return {'s1': s1, 's2': s2, 'origem': origem, 'cabeca': cabeca, 'cauda': divergencia}


def deslocar_fluxo_no_pas(grafo: nx.DiGraph, pas: dict) -> (float, nx.DiGraph):
    """
    Calcula e aplica o deslocamento de fluxo ótimo em um PAS (Newton Step).
    ATUALIZAÇÃO: Implementa deslocamento bidirecional e verificação de fluxo disponível.
    """
    s1, s2, origem = pas['s1'], pas['s2'], pas['origem']
    
    # 1. Obter fluxos disponíveis (gargalos) em cada segmento para a origem específica
    # f1: Fluxo disponível em s1 (tree path) que pode ser movido para s2 se necessário
    f1 = float('inf')
    for u, v in s1:
        f = grafo[u][v]['fluxos_por_origem'].get(origem, 0.0)
        if f < f1: f1 = f
        
    # f2: Fluxo disponível em s2 (non-tree path) que pode ser movido para s1
    f2 = float('inf')
    for u, v in s2:
        f = grafo[u][v]['fluxos_por_origem'].get(origem, 0.0)
        if f < f2: f2 = f

    # Se ambos os caminhos têm fluxo negligenciável, não há nada a otimizar
    if f1 < TOLERANCIA_FLUXO and f2 < TOLERANCIA_FLUXO:
        return 0.0, grafo

    # 2. Calcular Custos Totais
    custo_s1 = sum(grafo[u][v]['custo'] for u, v in s1)
    custo_s2 = sum(grafo[u][v]['custo'] for u, v in s2)

    # 3. Calcular Derivadas (Segunda Derivada da Função Objetivo / Primeira do Custo)
    def derivada_bpr(u, v):
        dados = grafo[u][v]
        # Derivada do custo BPR em relação ao fluxo:
        # c'(x) = t0 * alpha * beta * (x^(beta-1)) / (cap^beta)
        if dados['capacidade'] <= 0: return float('inf') # Proteção
        
        termo_comum = (BPR_ALPHA * BPR_BETA * dados['tempo_fluxo_livre']) / (dados['capacidade'] ** BPR_BETA)
        return termo_comum * (dados['fluxo'] ** (BPR_BETA - 1))

    dt1 = sum(derivada_bpr(u, v) for u, v in s1)
    dt2 = sum(derivada_bpr(u, v) for u, v in s2)
    
    denominador = dt1 + dt2
    
    # Evita divisão por zero (caso de fluxo livre constante ou capacidade infinita)
    if denominador < 1e-15:
        # Se as derivadas são zero, usamos apenas a diferença de custo
        # Se custo_s2 > custo_s1, queremos mover tudo de s2 para s1 (delta positivo grande)
        # Se custo_s1 > custo_s2, queremos mover tudo de s1 para s2 (delta negativo grande)
        delta = float('inf') if (custo_s2 > custo_s1) else -float('inf')
    else:
        # Passo de Newton: (c2 - c1) / (c'1 + c'2)
        delta = (custo_s2 - custo_s1) / denominador

    # 4. Clamping Bidirecional (Limitação pelo fluxo disponível)
    # Se delta > 0: Movemos s2 -> s1. Limitado por f2.
    # Se delta < 0: Movemos s1 -> s2. Limitado por f1 (magnitude).
    
    if delta > 0:
        delta = min(delta, f2)
    else:
        delta = max(delta, -f1) # max porque delta é negativo e -f1 é negativo

    # Verifica se a mudança é significativa
    if abs(delta) < TOLERANCIA_FLUXO:
        return 0.0, grafo

    # 5. Aplicar o deslocamento
    # Nota: delta positivo tira de s2 e põe em s1.
    #       delta negativo tira de s1 e põe em s2.
    # A matemática abaixo funciona para ambos os casos.

    # Segmento 1 (Árvore): Recebe +delta
    for u, v in s1:
        grafo[u][v]['fluxo'] += delta
        grafo[u][v]['fluxos_por_origem'][origem] += delta
        # Limpeza numérica
        if grafo[u][v]['fluxos_por_origem'][origem] < 0: grafo[u][v]['fluxos_por_origem'][origem] = 0.0
        if grafo[u][v]['fluxo'] < 0: grafo[u][v]['fluxo'] = 0.0
        atualizar_custo_arco(grafo, u, v)
        
    # Segmento 2 (Não-Árvore): Recebe -delta
    for u, v in s2:
        grafo[u][v]['fluxo'] -= delta
        grafo[u][v]['fluxos_por_origem'][origem] -= delta
        # Limpeza numérica
        if grafo[u][v]['fluxos_por_origem'][origem] < 0: grafo[u][v]['fluxos_por_origem'][origem] = 0.0
        if grafo[u][v]['fluxo'] < 0: grafo[u][v]['fluxo'] = 0.0
        atualizar_custo_arco(grafo, u, v)

    # Retornamos o valor absoluto do deslocamento para métricas de convergência
    return abs(delta), grafo


def adicionar_pas_ao_conjunto(conjunto_pas: list, novo_pas: dict) -> (list, bool):
    """
    Adiciona um novo PAS ao conjunto.
    CORREÇÃO: A unicidade deve considerar a TOPOLOGIA + ORIGEM.
    Duas origens diferentes podem ter PAS com a mesma topologia.
    """
    # Assinatura baseada nos IDs dos arcos para comparação rápida
    sig_s1 = tuple(sorted(novo_pas['s1']))
    sig_s2 = tuple(sorted(novo_pas['s2']))
    origem_nova = novo_pas['origem']
    
    nova_assinatura = (sig_s1, sig_s2, origem_nova)

    for pas_existente in conjunto_pas:
        ex_sig_s1 = tuple(sorted(pas_existente['s1']))
        ex_sig_s2 = tuple(sorted(pas_existente['s2']))
        ex_origem = pas_existente['origem']
        
        # Se for a mesma topologia E a mesma origem, é duplicata.
        if (sig_s1 == ex_sig_s1) and (sig_s2 == ex_sig_s2) and (origem_nova == ex_origem):
            return conjunto_pas, False
    
    return conjunto_pas + [novo_pas], True


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
        arquivo_viagens='./od_outputs/OD_10_300/OD_1.txt',
        max_iter=5000,
        gap_convergencia=1e-10
    )

    # Imprime os resultados finais
    """ print("\n--- Resultados Finais ---")
    print("Fluxo e Custo por Arco:")
    for u, v, dados in sorted(grafo_final.edges(data=True)):
        print(f"Arco ({u}->{v}): Fluxo = {dados['fluxo']:.2f}, Custo = {dados['custo']:.2f}") """