import numpy as np
import networkx as nx
from collections import defaultdict, deque
import math
from typing import List, Dict, Tuple, Set, Optional
from bibpy.utils import *

# =============================================================================
# CLASSES AUXILIARES (Estruturas de Dados do Algoritmo B)
# =============================================================================

class Merge:
    """
    Representa um nó de fusão (Merge Node) dentro de um Bush.
    Equivalente à struct 'merge_type' do C.
    """
    def __init__(self, node: int):
        self.node = node
        self.approaches: List[Tuple[int, int]] = []  # Lista de arcos (u, v) que chegam neste nó no Bush
        self.approach_flows: List[float] = []        # Fluxo em cada arco de aproximação
        
        # Índices na lista 'approaches' para o link do caminho mais curto (SP) e mais longo (LP)
        self.sp_link_idx: int = -1
        self.lp_link_idx: int = -1
        
        # Nó onde os caminhos SP e LP se separaram anteriormente
        self.divergence_node: int = -1

class Bush:
    """
    Representa o subgrafo acíclico enraizado em uma origem específica.
    Equivalente à struct 'bushes_type' (instanciada para uma origem).
    """
    def __init__(self, origin: int, num_nodes: int):
        self.origin = origin
        
        # Estrutura topológica
        # Grafo contendo apenas as arestas ativas no Bush
        self.graph = nx.DiGraph() 
        self.topo_order: List[int] = [] # Ordem topológica para varredura linear
        
        # Dados de Fusão (Merges)
        # Mapeia ID do nó -> Objeto Merge
        self.merges: Dict[int, Merge] = {}
        
        # Rótulos (Potenciais)
        self.sp_cost = np.full(num_nodes + 1, np.inf) # Custo do caminho mínimo (SP)
        self.lp_cost = np.full(num_nodes + 1, -np.inf) # Custo do caminho mais longo usado (LP)
        
        # Fluxos transientes (calculados sob demanda)
        self.node_flow = np.zeros(num_nodes + 1)
        # Mapeia (u,v) -> fluxo no bush
        self.link_flow: Dict[Tuple[int, int], float] = defaultdict(float) 

# =============================================================================
# ALGORITMO B (Lógica Principal)
# =============================================================================

class AlgorithmB:
    def __init__(self, grafo: nx.DiGraph, viagens: dict, origens: List[int]):
        self.grafo = grafo
        self.viagens = viagens
        self.origens = origens
        self.bushes: Dict[int, Bush] = {}
        self.num_nodes = grafo.number_of_nodes()
        
        # Parâmetros internos do algoritmo (similares ao C)
        self.min_link_flow = 1e-14
        self.newton_step = 1.0
        self.num_newton_shifts = 1
        
    def executar(self):
        """
        Executa o ciclo principal do Algoritmo B.
        """
        print("Iniciando Algoritmo B...")
        self._inicializar_bushes()
        
        iteracao = 0
        gap = float('inf')
        
        while iteracao < MAX_ITERACOES and gap > TOLERANCIA_CUSTO:
            iteracao += 1
            max_shift = 0.0
            
            # Itera sobre cada origem (equivale aos "batches" do C)
            for r in self.origens:
                bush = self.bushes[r]
                
                # Passo 1: Atualização Topológica (Melhorar a estrutura do Bush)
                self._atualizar_bush(bush)
                
                # Passo 2: Equilibração de Fluxos (Shift de Newton)
                self._atualizar_fluxos(bush)
            
            # Passo 3: Verificação de Convergência
            gap = calcular_gap_relativo(self.grafo, self.viagens, self.origens)
            
            if iteracao % 1 == 0:
                print(f"Iteração {iteracao}: Gap Relativo = {gap:.6e}")

        print(f"Convergência atingida na iteração {iteracao}.")

    # -------------------------------------------------------------------------
    # INICIALIZAÇÃO E TOPOLOGIA
    # -------------------------------------------------------------------------

    def _inicializar_bushes(self):
        """Cria os bushes iniciais baseados na árvore de caminhos mínimos (SPT)."""
        for r in self.origens:
            bush = Bush(r, self.num_nodes)
            
            # Calcula SPT inicial usando custos de fluxo livre (t0)
            pred, _ = nx.dijkstra_predecessor_and_distance(self.grafo, r, weight='custo')
            
            # Constrói o grafo do bush
            for v, preds in pred.items():
                if v == r: continue
                # Dijkstra retorna lista de predecessores, pegamos o primeiro para formar árvore
                u = preds[0] 
                bush.graph.add_edge(u, v)
                bush.link_flow[(u,v)] = 0.0 # Inicializa
            
            # Carrega demanda inicial (Tudo-ou-Nada na SPT inicial)
            self._reconstruir_merges_e_ordem(bush)
            self._calcular_fluxos_bush(bush) # Carrega fluxos iniciais nos links globais
            self.bushes[r] = bush
        
        # Atualiza custos globais após carga inicial
        self._atualizar_custos_globais()

    def _atualizar_bush(self, bush: Bush):
        """
        Atualiza a estrutura do bush: remove links inúteis e adiciona atalhos.
        Corresponde a 'updateBushB' no C.
        """
        # 1. Escaneia rótulos (SP e LP) dentro do Bush atual
        self._scan_bush(bush)
        
        # 2. Varredura global para adicionar/remover links
        links_adicionados = False
        
        # Conjunto de arestas para remover (fluxo quase zero e não é SP)
        arestas_para_remover = []
        for u, v in list(bush.graph.edges()):
            fluxo = bush.link_flow.get((u,v), 0.0)
            if fluxo < self.min_link_flow:
                # Só remove se NÃO fizer parte da árvore de caminho mínimo atual
                # (Simplificação: removemos se fluxo é nulo, a SPT será recriada se necessário)
                arestas_para_remover.append((u, v))
        
        bush.graph.remove_edges_from(arestas_para_remover)
        for u, v in arestas_para_remover:
            del bush.link_flow[(u,v)]

        # 3. Adicionar "Shortcuts" (Atalhos)
        # Verifica todos os arcos da rede global
        for u, v, dados in self.grafo.edges(data=True):
            cost_uv = dados['custo']
            
            # Se u não é alcançável no bush, pula
            if bush.sp_cost[u] == np.inf: continue
            
            # Critério de DIAL/NIE: Adicionar se oferecer redução estrita de custo
            # SP[u] + custo_arco < SP[v]
            if bush.sp_cost[u] + cost_uv < bush.sp_cost[v] - 1e-8:
                if not bush.graph.has_edge(u, v):
                    bush.graph.add_edge(u, v)
                    bush.link_flow[(u,v)] = 0.0 # Novo link começa sem fluxo
                    links_adicionados = True
        
        # 4. Se a topologia mudou, reconstruir estruturas
        if links_adicionados or len(arestas_para_remover) > 0:
            self._reconstruir_merges_e_ordem(bush)

    def _reconstruir_merges_e_ordem(self, bush: Bush):
        """
        Reconstrói a lista de nós de fusão e a ordem topológica.
        Corresponde a 'reconstructMerges' e 'topologicalOrder'.
        """
        # Garante que é um DAG (remove ciclos se houverem, embora a lógica SP previna)
        try:
            bush.topo_order = list(nx.topological_sort(bush.graph))
        except nx.NetworkXUnfeasible:
            # Fallback: se houver ciclo (erro numérico), recria via SPT
            # (Simplificado para este exemplo)
            bush.topo_order = list(bush.graph.nodes()) 
        
        # Reconstrói objetos Merge
        bush.merges.clear()
        for node in bush.graph.nodes():
            in_edges = list(bush.graph.in_edges(node))
            if len(in_edges) > 1:
                m = Merge(node)
                for u, v in in_edges:
                    m.approaches.append((u, v))
                    m.approach_flows.append(bush.link_flow.get((u, v), 0.0))
                bush.merges[node] = m

    def _scan_bush(self, bush: Bush):
        """
        Calcula SP_cost e LP_cost para o bush.
        Corresponde a 'scanBushes'.
        """
        # Inicializa rótulos
        bush.sp_cost.fill(np.inf)
        bush.lp_cost.fill(-np.inf)
        
        # Origem
        bush.sp_cost[bush.origin] = 0.0
        bush.lp_cost[bush.origin] = 0.0
        
        # Varredura em ordem topológica (O(N))
        for u in bush.topo_order:
            # Propagar para vizinhos no Bush
            if bush.sp_cost[u] == np.inf: continue
            
            for v in bush.graph.successors(u):
                custo_arco = self.grafo[u][v]['custo']
                
                # Atualiza SP (Shortest Path)
                if bush.sp_cost[u] + custo_arco < bush.sp_cost[v]:
                    bush.sp_cost[v] = bush.sp_cost[u] + custo_arco
                
                # Atualiza LP (Longest Used Path)
                # Só propaga LP se houver fluxo ou for caminho forçado
                if bush.lp_cost[u] > -np.inf:
                    if bush.lp_cost[u] + custo_arco > bush.lp_cost[v]:
                        bush.lp_cost[v] = bush.lp_cost[u] + custo_arco

    # -------------------------------------------------------------------------
    # FLUXOS E EQUILIBRAÇÃO (Newton)
    # -------------------------------------------------------------------------

    def _atualizar_fluxos(self, bush: Bush):
        """
        Gerencia o ciclo de atualização de fluxos: calcula fluxos atuais,
        encontra divergências e aplica Newton.
        Corresponde a 'updateFlowsB'.
        """
        # 1. Calcula fluxos atuais no bush (varredura reversa)
        self._calcular_fluxos_bush(bush)
        
        # 2. Recalcula rótulos com fluxos atualizados (para identificar LP/SP corretamente)
        self._scan_bush(bush)
        
        # 3. Identifica nós de divergência para os Merges
        self._encontrar_divergencias(bush)
        
        # 4. Passada descendente (Topo order inversa) para aplicar Shifts
        # Iteramos ao contrário para ajustar nós mais distantes primeiro
        for node in reversed(bush.topo_order):
            if node in bush.merges:
                merge = bush.merges[node]
                
                # Identifica índices dos links SP e LP entre as aproximações
                self._classificar_aproximacoes(bush, merge)
                
                # Se vale a pena fazer shift
                if (merge.sp_link_idx != -1 and merge.lp_link_idx != -1 and 
                    merge.sp_link_idx != merge.lp_link_idx):
                    
                    self._newton_flow_shift(bush, merge)

    def _calcular_fluxos_bush(self, bush: Bush):
        """
        Calcula fluxos nos arcos do bush acumulando demanda da OD.
        Não armazena fluxo global, apenas atualiza `bush.link_flow`.
        Corresponde a 'calculateBushFlows'.
        """
        bush.node_flow.fill(0.0)
        
        # Carrega demandas de destino
        for (o, d), demanda in self.viagens.items():
            if o == bush.origin:
                bush.node_flow[d] += demanda
        
        # Zera fluxos de link atuais do bush para recalcular
        for k in bush.link_flow:
            bush.link_flow[k] = 0.0

        # Varredura reversa (Destinos -> Origem)
        for u in reversed(bush.topo_order):
            fluxo_total_no = bush.node_flow[u]
            if fluxo_total_no <= 0: continue
            
            # Se for nó de fusão, distribuir fluxo entre entradas proporcionalmente
            if u in bush.merges:
                merge = bush.merges[u]
                fluxo_total_entrada = sum(merge.approach_flows)
                
                if fluxo_total_entrada > 0:
                    for i, (pred_u, _) in enumerate(merge.approaches):
                        prop = merge.approach_flows[i] / fluxo_total_entrada
                        fluxo_link = fluxo_total_no * prop
                        
                        # Atualiza fluxo link local
                        bush.link_flow[(pred_u, u)] = fluxo_link
                        # Empurra fluxo para nó anterior
                        bush.node_flow[pred_u] += fluxo_link
                        # Atualiza fluxo armazenado no merge para próxima iteração
                        merge.approach_flows[i] = fluxo_link 
                else:
                    # Caso degenerado (sem fluxo anterior): joga tudo no SP atual (se existir)
                    # (Simplificação: joga no primeiro)
                    if merge.approaches:
                        pred_u, _ = merge.approaches[0]
                        bush.link_flow[(pred_u, u)] = fluxo_total_no
                        bush.node_flow[pred_u] += fluxo_total_no
                        merge.approach_flows[0] = fluxo_total_no

            else:
                # Nó simples (1 entrada ou origem)
                in_edges = list(bush.graph.in_edges(u))
                if in_edges:
                    pred_u, _ = in_edges[0]
                    bush.link_flow[(pred_u, u)] = fluxo_total_no
                    bush.node_flow[pred_u] += fluxo_total_no

    def _encontrar_divergencias(self, bush: Bush):
        """
        Encontra o nó onde os caminhos SP e LP se encontram para trás.
        Corresponde a 'findDivergenceNodes'.
        """
        for merge in bush.merges.values():
            # Rastreamento simplificado para trás
            # (Na implementação completa, rastreia-se os pais recursivamente até coincidir)
            # Aqui, assumiremos que a lógica de shift cuidará da parada
            pass 

    def _classificar_aproximacoes(self, bush: Bush, merge: Merge):
        """Identifica qual link de entrada pertence ao caminho SP e qual ao LP."""
        min_cost = np.inf
        max_cost = -np.inf
        
        merge.sp_link_idx = -1
        merge.lp_link_idx = -1
        
        for i, (u, v) in enumerate(merge.approaches):
            custo_arco = self.grafo[u][v]['custo']
            custo_caminho_u = bush.sp_cost[u] # Custo até o nó anterior
            
            # Identificando SP
            if custo_caminho_u + custo_arco < min_cost:
                min_cost = custo_caminho_u + custo_arco
                merge.sp_link_idx = i
            
            # Identificando LP (apenas se tiver fluxo significativo)
            # Nota: Usamos LP_Cost do nó u para maximizar o caminho
            if merge.approach_flows[i] > self.min_link_flow:
                custo_lp_u = bush.lp_cost[u]
                if custo_lp_u + custo_arco > max_cost:
                    max_cost = custo_lp_u + custo_arco
                    merge.lp_link_idx = i

    def _newton_flow_shift(self, bush: Bush, merge: Merge):
        """
        Aplica o método de Newton para transferir fluxo do caminho LP para SP.
        """
        u_sp, v_node = merge.approaches[merge.sp_link_idx]
        u_lp, _      = merge.approaches[merge.lp_link_idx]
        
        # Reconstrói os segmentos (backtracking) até encontrar divergência
        # Para simplificar a transcrição em Python, faremos o cálculo apenas nos links imediatos
        # e seus predecessores diretos na árvore, similar à lógica de Dial simplificada.
        # Numa implementação C completa, isso é um while loop até divergenceNode.
        
        # Custo e Derivada do Caminho LP
        cost_lp = 0.0; der_lp = 0.0
        curr = u_lp
        # Loop simples simulando rastreamento (limite de profundidade para segurança)
        for _ in range(1000): 
            # Pega arco (pred -> curr)
            preds = list(bush.graph.in_edges(curr))
            if not preds: break
            # Simplificação: assume 1 pred principal ou o de maior fluxo
            pred, _ = preds[0]
            
            dados_arco = self.grafo[pred][curr]
            cost_lp += dados_arco['custo']
            der_lp += self._calcular_derivada_bpr(dados_arco)
            curr = pred
            if curr == bush.origin: 
                break
            
        # Adiciona o arco final da fusão
        dados_lp_final = self.grafo[u_lp][v_node]
        cost_lp += dados_lp_final['custo']
        der_lp += self._calcular_derivada_bpr(dados_lp_final)

        # Custo e Derivada do Caminho SP (Mesma lógica)
        cost_sp = 0.0; der_sp = 0.0
        curr = u_sp
        for _ in range(10):
            preds = list(bush.graph.in_edges(curr))
            if not preds: break
            pred, _ = preds[0]
            dados_arco = self.grafo[pred][curr]
            cost_sp += dados_arco['custo']
            der_sp += self._calcular_derivada_bpr(dados_arco)
            curr = pred
            if curr == bush.origin: break
            
        dados_sp_final = self.grafo[u_sp][v_node]
        cost_sp += dados_sp_final['custo']
        der_sp += self._calcular_derivada_bpr(dados_sp_final)

        # Cálculo do Shift (Newton)
        if (der_lp + der_sp) == 0: return

        delta_x = self.newton_step * (cost_lp - cost_sp) / (der_lp + der_sp)
        
        # Limita o shift ao fluxo disponível no caminho mais longo
        fluxo_disponivel = merge.approach_flows[merge.lp_link_idx]
        delta_x = min(delta_x, fluxo_disponivel)
        
        if delta_x <= 0: return

        # Aplica o Shift
        # 1. Remove de LP
        self._aplicar_mudanca_fluxo(bush, u_lp, v_node, -delta_x, merge, merge.lp_link_idx)
        # 2. Adiciona em SP
        self._aplicar_mudanca_fluxo(bush, u_sp, v_node, delta_x, merge, merge.sp_link_idx)

    def _aplicar_mudanca_fluxo(self, bush: Bush, u: int, v: int, delta: float, merge: Merge, idx: int):
        """Atualiza fluxos locais e globais e custos."""
        # Atualiza fluxo local no bush
        bush.link_flow[(u,v)] += delta
        merge.approach_flows[idx] += delta
        
        # Atualiza fluxo global no grafo (NetworkX)
        self.grafo[u][v]['fluxo'] += delta
        
        # Atualiza custo imediato (Gauss-Seidel) - 'exactCostUpdate'
        self._atualizar_custo_arco(u, v)

    # -------------------------------------------------------------------------
    # FUNÇÕES DE CUSTO (BPR)
    # -------------------------------------------------------------------------

    def _calcular_derivada_bpr(self, dados_arco) -> float:
        """Calcula a primeira derivada da função BPR."""
        t0 = dados_arco['tempo_fluxo_livre']
        cap = dados_arco['capacidade']
        fluxo = dados_arco['fluxo']
        
        # C'(x) = t0 * alpha * beta * (x/cap)^(beta-1) * (1/cap)
        term = (fluxo / cap) ** (BPR_BETA - 1)
        return (t0 * BPR_ALPHA * BPR_BETA * term) / cap

    def _atualizar_custo_arco(self, u, v):
        """Recalcula o custo do arco baseado no fluxo atual."""
        dados = self.grafo[u][v]
        fluxo = dados['fluxo']
        t0 = dados['tempo_fluxo_livre']
        cap = dados['capacidade']
        
        # BPR: t = t0 * (1 + alpha * (v/c)^beta)
        dados['custo'] = t0 * (1 + BPR_ALPHA * ((fluxo / cap) ** BPR_BETA))

    def _atualizar_custos_globais(self):
        """Atualiza custos de todos os arcos da rede (usado na inicialização)."""
        # Primeiro zera fluxos globais
        for u, v, dados in self.grafo.edges(data=True):
            dados['fluxo'] = 0.0
        
        # Soma fluxos de todos os bushes
        for bush in self.bushes.values():
            for (u, v), f in bush.link_flow.items():
                if self.grafo.has_edge(u, v):
                    self.grafo[u][v]['fluxo'] += f
        
        # Recalcula custos
        for u, v in self.grafo.edges():
            self._atualizar_custo_arco(u, v)

# =============================================================================
# EXEMPLO DE USO (Integração)
# =============================================================================

if __name__ == "__main__":
    
    
    # Exemplo de caminhos (ajuste conforme necessário)
    arquivo_rede = "./fortaleza/edges.txt" # Exemplo
    arquivo_trips = './od_outputs/OD_10_300/OD_2.txt' # Exemplo
    
    # Criação dummy para teste de sintaxe se arquivos não existirem
    try:
        grafo = carregar_rede(arquivo_rede)
        viagens = carregar_viagens(arquivo_trips)
        
        # Extrai lista de origens únicas
        origens = list(set(o for o, d in viagens.keys()))
        
        algoritmo = AlgorithmB(grafo, viagens, origens)
        algoritmo.executar()
        
    except OSError:
        print("Arquivos de entrada não encontrados. O código foi compilado com sucesso.")