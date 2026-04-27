"""
ADMM-JOR: Parallel Computing Framework for Traffic Assignment Problem

Implementation of the algorithm described in:
"A novel parallel computing framework for traffic assignment problem: Integrating 
alternating direction method of multipliers with Jacobi over relaxation method"
(Liu et al., 2024)

This module solves the deterministic user equilibrium (DUE) traffic assignment problem
using the ADMM-JOR algorithm with adaptive relaxation factor.

Time Complexity: O(P * I * |A|) where P is number of blocks, I is iterations, |A| is edges
Space Complexity: O(|O| * |A|) for storing origin-based link flows
"""

import numpy as np
import networkx as nx
from typing import Dict, Tuple, List
from collections import defaultdict
import time

# Parameters
BPR_ALPHA = 0.15
BPR_BETA = 4.0
TOLERANCIA_GAP = 1e-6
MAX_ITERACOES = 500
RHO_PENALTY = 1.0
TOLERANCIA_DUAL = 1e-8


# =============================================================================
# 1. DATA STRUCTURES
# =============================================================================

class NetworkState:
    """Stores the state of network flow and costs."""
    
    def __init__(self, graph: nx.DiGraph, origins: List[int], destinations: List[int]):
        self.graph = graph
        self.origins = origins
        self.destinations = destinations
        self.num_nodes = graph.number_of_nodes()
        self.num_edges = graph.number_of_edges()
        self.num_origins = len(origins)
        
        # Origin-based link flows: v_a^o for each link a and origin o
        self.flows = {}  # (arco, origem) -> fluxo
        self.flows_iter = {}  # flows from ADMM iteration
        
        # Dual variables (Lagrange multipliers)
        self.lambda_dual = {}  # (nó, origem) -> multiplicador
        
        # Link costs
        self.link_costs = {}  # arco -> custo
        
        self._initialize_flows_and_duals()
    
    def _initialize_flows_and_duals(self):
        """Initialize all flows and dual variables to zero."""
        for arc in self.graph.edges():
            for origin in self.origins:
                self.flows[(arc, origin)] = 0.0
                self.flows_iter[(arc, origin)] = 0.0
        
        for node in self.graph.nodes():
            for origin in self.origins:
                self.lambda_dual[(node, origin)] = 0.0
        
        for arc in self.graph.edges():
            self.link_costs[arc] = self.graph.edges[arc]['tempo_fluxo_livre']


# =============================================================================
# 2. LINK COLORING (Edge Coloring for Block Decomposition)
# =============================================================================

def color_edges_greedy(graph: nx.DiGraph) -> Dict[Tuple, int]:
    """
    Greedy edge coloring algorithm to partition links into P non-adjacent blocks.
    
    Returns:
        Dict mapping each edge to its color (block number)
    """
    edges_list = list(graph.edges())
    edge_colors = {}
    node_occupied_colors = defaultdict(set)
    
    for edge in edges_list:
        u, v = edge
        used_colors = node_occupied_colors[u] | node_occupied_colors[v]
        
        # Find the smallest color not used by adjacent nodes
        color = 0
        while color in used_colors:
            color += 1
        
        edge_colors[edge] = color
        node_occupied_colors[u].add(color)
        node_occupied_colors[v].add(color)
    
    return edge_colors


def partition_edges_into_blocks(graph: nx.DiGraph) -> Tuple[Dict[int, List], int]:
    """
    Partition edges into independent blocks using edge coloring.
    
    Returns:
        (blocks, num_blocks) where blocks[p] is list of edges in block p
    """
    edge_colors = color_edges_greedy(graph)
    blocks = defaultdict(list)
    
    for edge, color in edge_colors.items():
        blocks[color].append(edge)
    
    num_blocks = max(blocks.keys()) + 1 if blocks else 1
    return dict(blocks), num_blocks


# =============================================================================
# 3. SHORTEST PATH TREES AND FLOW INITIALIZATION
# =============================================================================

def compute_spt(graph: nx.DiGraph, origin: int) -> Dict[int, float]:
    """
    Compute shortest path tree distances from origin using Dijkstra.
    
    Returns:
        Dict mapping destination nodes to shortest distances
    """
    try:
        _, distances = nx.dijkstra_predecessor_and_distance(
            graph, source=origin, weight='custo'
        )
        return distances
    except nx.NodeNotFound:
        return {}


def initialize_flows_frank_wolfe(
    state: NetworkState,
    viagens: Dict[Tuple[int, int], float],
    num_iter_fw: int = 1
) -> None:
    """
    Initialize flows using shortest path method to ensure flow feasibility.
    All demand is routed along shortest paths to guarantee conservation.
    """
    print(f"  Initializing flows by routing all demand along shortest paths...")
    
    # Zerar fluxos iniciais
    for arc in state.graph.edges():
        for origin in state.origins:
            state.flows[(arc, origin)] = 0.0
    
    # Para cada OD, rotear demanda no caminho mais curto
    for origin in state.origins:
        try:
            pred, distances = nx.dijkstra_predecessor_and_distance(
                state.graph, source=origin, weight='tempo_fluxo_livre'
            )
            
            # Para cada destino deste origin
            for (o, d), demand in viagens.items():
                if o == origin and d in distances and demand > 0:
                    # Reconstruir caminho mais curto de origin até d
                    path_nodes = [d]
                    current = d
                    while current != origin:
                        if current not in pred or len(pred[current]) == 0:
                            path_nodes = []
                            break
                        current = pred[current][0]
                        path_nodes.append(current)
                    
                    if not path_nodes:
                        continue
                    
                    path_nodes.reverse()
                    # Converter nós em arcos e atribuir demanda
                    for i in range(len(path_nodes) - 1):
                        arc = (path_nodes[i], path_nodes[i+1])
                        if arc in state.graph.edges():
                            state.flows[(arc, origin)] += demand
        except Exception:
            continue
    
    # Atualizar custos com base nos fluxos iniciais
    update_link_costs(state)


# =============================================================================
# 4. LINK COST FUNCTIONS (BPR Function)
# =============================================================================

def compute_bpr_cost(free_flow_time: float, flow: float, capacity: float) -> float:
    """
    Compute link cost using Bureau of Public Roads (BPR) function.
    
    t(x) = t_0 * (1 + alpha * (x/c)^beta)
    """
    if capacity <= 0:
        return free_flow_time
    
    ratio = flow / capacity
    cost = free_flow_time * (1.0 + BPR_ALPHA * (ratio ** BPR_BETA))
    return cost


def compute_bpr_derivative(free_flow_time: float, flow: float, capacity: float) -> float:
    """Derivative of BPR with respect to flow (for optimization)."""
    if capacity <= 0:
        return 0.0
    
    ratio = flow / capacity
    derivative = free_flow_time * BPR_ALPHA * BPR_BETA * (ratio ** (BPR_BETA - 1)) / capacity
    return derivative


def update_link_costs(state: NetworkState) -> None:
    """Update all link costs based on current aggregate flows."""
    for arc in state.graph.edges():
        # Aggregate flow from all origins
        total_flow = sum(
            state.flows.get((arc, origin), 0.0)
            for origin in state.origins
        )
        
        capacity = state.graph.edges[arc]['capacidade']
        free_flow_time = state.graph.edges[arc]['tempo_fluxo_livre']
        
        cost = compute_bpr_cost(free_flow_time, total_flow, capacity)
        state.link_costs[arc] = cost
        # IMPORTANTE: Atualizar também o grafo para usar nas próximas iterações
        state.graph.edges[arc]['custo'] = cost
        state.graph.edges[arc]['fluxo'] = total_flow


# =============================================================================
# 5. FLOW CONSERVATION CONSTRAINTS
# =============================================================================

def compute_flow_conservation_residual(
    state: NetworkState,
    origin: int,
    node: int,
    viagens: Dict[Tuple[int, int], float]
) -> float:
    """
    Compute the flow conservation constraint residual at a node for an origin.
    
    H_n^o = sum(flow_a_in) - sum(flow_a_out) - g_n^o
    
    where g_n^o is the demand (positive at origin, negative at destination, 0 elsewhere)
    """
    # Outgoing flows
    outgoing = sum(
        state.flows.get(((node, v), origin), 0.0)
        for v in state.graph.successors(node)
    )
    
    # Incoming flows
    incoming = sum(
        state.flows.get(((u, node), origin), 0.0)
        for u in state.graph.predecessors(node)
    )
    
    # Demand
    if node == origin:
        demand = sum(viagens.get((origin, d), 0.0) for d in state.destinations)
    else:
        demand = -viagens.get((origin, node), 0.0)
    
    residual = outgoing - incoming - demand
    return residual


# =============================================================================
# 6A. ALL-OR-NOTHING ASSIGNMENT (FOR STABLE UE SOLUTION)
# =============================================================================

def all_or_nothing_assignment(
    state: NetworkState,
    viagens: Dict[Tuple[int, int], float]
) -> Dict[Tuple[Tuple[int, int], int], float]:
    """Compute an all-or-nothing assignment based on current link costs.

    For each OD pair, all demand is loaded on the current shortest path.
    Returns origin-based link flows v_a^o.
    """
    aon_flows: Dict[Tuple[Tuple[int, int], int], float] = {}
    
    for arc in state.graph.edges():
        for origin in state.origins:
            aon_flows[(arc, origin)] = 0.0
    
    for origin in state.origins:
        try:
            pred, distances = nx.dijkstra_predecessor_and_distance(
                state.graph, source=origin, weight='custo'
            )
        except Exception:
            continue
        
        for (o, d), demand in viagens.items():
            if o != origin or demand <= 0:
                continue
            if d not in distances:
                continue
            
            # Reconstruir caminho mais curto de origin até d
            path_nodes = [d]
            current = d
            while current != origin:
                if current not in pred or len(pred[current]) == 0:
                    path_nodes = []
                    break
                current = pred[current][0]
                path_nodes.append(current)
            
            if not path_nodes:
                continue
            
            path_nodes.reverse()
            for i in range(len(path_nodes) - 1):
                arc = (path_nodes[i], path_nodes[i+1])
                if arc in state.graph.edges():
                    aon_flows[(arc, origin)] += demand
    
    return aon_flows


# =============================================================================
# 6. AUGMENTED LAGRANGIAN AND OPTIMIZATION
# =============================================================================

def compute_augmented_lagrangian_link(
    arc: Tuple[int, int],
    origin: int,
    state: NetworkState,
    rho: float,
    viagens: Dict[Tuple[int, int], float]
) -> float:
    """
    Compute augmented Lagrangian for a single link-origin pair.
    
    L_rho = integral_cost + lambda_term + penalty_term
    """
    v_ao = state.flows.get((arc, origin), 0.0)
    u, v = arc
    
    # Integral of cost function (objective)
    capacity = state.graph.edges[arc]['capacidade']
    free_flow_time = state.graph.edges[arc]['tempo_fluxo_livre']
    integral_cost = free_flow_time * v_ao + \
                    BPR_ALPHA * free_flow_time / (BPR_BETA + 1) * \
                    (v_ao ** (BPR_BETA + 1)) / (capacity ** BPR_BETA)
    
    # Compute flow conservation residuals at both ends
    residual_u = compute_flow_conservation_residual(state, origin, u, viagens)
    residual_v = compute_flow_conservation_residual(state, origin, v, viagens)
    
    # Lambda terms (simplified - each endpoint contributes)
    lambda_term = state.lambda_dual.get((u, origin), 0.0) * (-residual_u) + \
                  state.lambda_dual.get((v, origin), 0.0) * (-residual_v)
    
    # Penalty term
    penalty = 0.5 * rho * (residual_u**2 + residual_v**2)
    
    return integral_cost + lambda_term + penalty


def optimize_link_flows_block(
    block_edges: List[Tuple[int, int]],
    origin: int,
    state: NetworkState,
    rho: float,
    viagens: Dict[Tuple[int, int], float],
    lambda_old: Dict,
) -> Dict:
    """
    Optimize flows for all links in a block for a given origin.
    Uses gradient descent / projected gradient method.
    
    Returns:
        Dictionary of updated flows for this block-origin pair
    """
    updated_flows = {}
    learning_rate = 0.05  # Reduzido para convergência mais estável
    num_inner_iter = 30   # Aumentado para melhor otimização local
    
    for arc in block_edges:
        v_ao = state.flows.get((arc, origin), 0.0)
        u, v = arc
        capacity = state.graph.edges[arc]['capacidade']
        free_flow_time = state.graph.edges[arc]['tempo_fluxo_livre']
        
        # Gradient-based optimization
        for inner_iter in range(num_inner_iter):
            # Gradient of objective + penalty
            residual_u = compute_flow_conservation_residual(state, origin, u, viagens)
            residual_v = compute_flow_conservation_residual(state, origin, v, viagens)
            
            # O gradiente da função objetivo do DUE em relação a v_a^o é o próprio custo do arco t_a(v_a)
            # Calculamos o fluxo total atualizando a parte deste origin
            total_flow = sum(state.flows.get((arc, o), 0.0) for o in state.origins)
            grad_obj = compute_bpr_cost(free_flow_time, total_flow, capacity)

            
            grad_dual = -state.lambda_dual.get((u, origin), 0.0) + \
                        state.lambda_dual.get((v, origin), 0.0)
            
            grad_penalty = rho * (residual_u - residual_v)
            
            grad_total = grad_obj + grad_dual + grad_penalty
            
            # Update com step size adaptativo
            v_ao_new = max(0.0, v_ao - learning_rate * grad_total)
            
            # Adaptive learning rate adjustment
            if inner_iter % 10 == 0 and inner_iter > 0:
                if abs(v_ao_new - v_ao) < 1e-9:
                    break
                learning_rate *= 0.9  # Reduzir step size ao longo do tempo
            
            v_ao = v_ao_new
        
        updated_flows[(arc, origin)] = v_ao
    
    return updated_flows


# =============================================================================
# 7. ADMM-JOR ALGORITHM
# =============================================================================

def solve_due_admm_jor(
    graph: nx.DiGraph,
    viagens: Dict[Tuple[int, int], float],
    rho: float = RHO_PENALTY,
    max_iter: int = MAX_ITERACOES,
    tolerance: float = TOLERANCIA_GAP,
    omega: float = 1.2,
    adaptive_omega: bool = True
) -> Tuple[Dict, float, int]:
    """
    Solve DUE-TAP using ADMM-JOR algorithm.
    
    Parameters:
        graph: Network topology
        viagens: OD demand matrix
        rho: Penalty parameter
        max_iter: Maximum iterations
        tolerance: Convergence tolerance (gap)
        omega: Relaxation factor for JOR
        adaptive_omega: Whether to adaptively adjust omega
    
    Returns:
        (flows, gap, iterations) where:
            flows: Dict[(arc, origin)] -> flow
            gap: Final relative gap
            iterations: Number of iterations executed
    """
    print("\n" + "="*70)
    print("ADMM-JOR ALGORITHM FOR DUE-TAP")
    print("="*70)
    
    start_time = time.time()
    
    # Extract origins and destinations
    origins = sorted(set(o for o, d in viagens.keys()))
    destinations = sorted(set(d for o, d in viagens.keys()))
    
    # Initialize state
    state = NetworkState(graph, origins, destinations)

    # Inicializar fluxos com atribuição de caminho mais curto
    print(f"\nInitializing flows for {len(origins)} origins...")
    initialize_flows_frank_wolfe(state, viagens, num_iter_fw=1)

    # Inicializar a estratégia paralela de decomposição de blocos
    blocks_dict, num_blocks = partition_edges_into_blocks(graph)
    print(f"Network partitioned into {num_blocks} blocks via greedy edge coloring (ADMM-JOR).")

    # Laço principal do algoritmo ADMM-JOR (Zhiyuan24.md)
    gaps: List[float] = []

    for iteration in range(max_iter):
        # Atualizar custos de link com base nos fluxos atuais
        update_link_costs(state)

        # Métrica de convergência: sempre avaliada ANTES da atualização de fluxo
        gap = compute_relative_gap(state, viagens, origins, destinations)
        gaps.append(gap)

        total_flow = sum(state.flows.values())
        # Print mais frequentemente para acompanhar a convergência
        if iteration % 5 == 0 or iteration < 15:
            print(f"  Iter {iteration:4d}: Gap = {gap:.6e}, Flow = {total_flow:.2f}")

        # Forçar pelo menos MIN_ITERATIONS antes de permitir parada
        # Isso garante que o algoritmo faz refinamentos suficientes
        MIN_ITERATIONS = 10
        if gap < tolerance and iteration >= MIN_ITERATIONS:
            print(f"\n✓ Converged at iteration {iteration}")
            break

        # Salvar fluxos e lambda para os cálculos (v^k, lambda^k) antigos
        v_old = state.flows.copy()
        lambda_old = state.lambda_dual.copy()
        f_old = compute_objective_total(state)
        
        v_tilde = state.flows.copy()
        lambda_tilde = state.lambda_dual.copy()

        # FASE 1: Predictor (ADMM Primal e Dual)
        # 1.1 Atualizar fluxos primais partindo os links em blocos (Gauss-Seidel entre blocos)
        for p in range(num_blocks):
            block_edges = blocks_dict[p]
            for origin in origins:
                updated_flows = optimize_link_flows_block(
                    block_edges, origin, state, rho, viagens, state.lambda_dual
                )
                for k, new_flow in updated_flows.items():
                    state.flows[k] = new_flow
                    v_tilde[k] = new_flow
        
        # 1.2 Atualizar variável dual preditora
        for origin in origins:
            for node in graph.nodes():
                residual = compute_flow_conservation_residual(state, origin, node, viagens)
                lambda_tilde[(node, origin)] = state.lambda_dual[(node, origin)] + rho * residual

        # FASE 2: Corrector Diferenciado (Relaxação JOR com omega autoajustável)
        w_curr = omega
        if adaptive_omega:
            eta = 0.75
            while True:
                # Testando o parâmetro JOR e avaliando a função objetivo iterativa
                for key in state.flows.keys():
                    state.flows[key] = w_curr * v_tilde[key] + (1.0 - w_curr) * v_old[key]
                update_link_costs(state)
                f_new = compute_objective_total(state)
                
                # Critério de parada: se a função objetivo decair ou o w_curr for muito brando
                if f_new < f_old or w_curr < 0.05:
                    break
                w_curr *= eta
        
        # 2.1 Aplicar as atualizações definitivas JOR sobre o estado real primário e dual
        for key in state.flows.keys():
            state.flows[key] = w_curr * v_tilde[key] + (1.0 - w_curr) * v_old[key]
        for key in state.lambda_dual.keys():
            state.lambda_dual[key] = w_curr * lambda_tilde[key] + (1.0 - w_curr) * lambda_old[key]
    
    elapsed = time.time() - start_time
    
    print(f"\nFinal relative gap: {gap:.6e}")
    print(f"Total iterations: {iteration + 1}")
    print(f"Computation time: {elapsed:.2f} seconds")
    print("="*70 + "\n")
    
    return state.flows, gap, iteration + 1


# =============================================================================
# 8. CONVERGENCE AND VALIDATION
# =============================================================================

def compute_objective_total(state: NetworkState) -> float:
    """Compute the total objective function for User Equilibrium (integral of link costs)."""
    obj = 0.0
    for arc in state.graph.edges():
        flow = sum(state.flows.get((arc, origin), 0.0) for origin in state.origins)
        capacity = state.graph.edges[arc]['capacidade']
        t0 = state.graph.edges[arc]['tempo_fluxo_livre']
        
        if capacity > 0:
            integral = t0 * flow + t0 * BPR_ALPHA * (flow ** (BPR_BETA + 1.0)) / ((BPR_BETA + 1.0) * (capacity ** BPR_BETA))
        else:
            integral = t0 * flow
        obj += integral
    return obj


def compute_relative_gap(
    state: NetworkState,
    viagens: Dict[Tuple[int, int], float],
    origins: List[int],
    destinations: List[int]
) -> float:
    """
    Compute Relative Gap - standard convergence metric for DUE.
    
    Gap = (Z_UE - Z_LB) / Z_UE
    
    where Z_UE is system travel time at UE, Z_LB is lower bound
    """
    # Current system travel time
    z_ue = 0.0
    for arc in state.graph.edges():
        flow = sum(
            state.flows.get((arc, origin), 0.0)
            for origin in origins
        )
        cost = state.link_costs.get(arc, 0.0)
        z_ue += flow * cost
    
    # Lower bound: shortest path costs
    z_lb = 0.0
    for origin in origins:
        distances = compute_spt(state.graph, origin)
        for (o, d), demand in viagens.items():
            if o == origin and d in distances:
                z_lb += demand * distances[d]
    
    if z_ue <= 1e-6:
        return 1.0
    
    gap = (z_ue - z_lb) / z_ue
    return max(0.0, gap)


def validate_flow_conservation(
    state: NetworkState,
    viagens: Dict[Tuple[int, int], float]
) -> bool:
    """Check if flow conservation constraints are satisfied."""
    print("\nValidating flow conservation constraints...")
    
    max_violation = 0.0
    for origin in state.origins:
        for node in state.graph.nodes():
            residual = compute_flow_conservation_residual(
                state, origin, node, viagens
            )
            max_violation = max(max_violation, abs(residual))
    
    if max_violation < TOLERANCIA_DUAL:
        print(f"✓ Flow conservation satisfied (max violation: {max_violation:.2e})")
        return True
    else:
        print(f"✗ Flow conservation NOT satisfied (max violation: {max_violation:.2e})")
        return False


# =============================================================================
# 9. SOLUTION EXTRACTION AND REPORTING
# =============================================================================

def extract_link_flows(state: NetworkState) -> Dict[Tuple[int, int], float]:
    """Extract aggregate flows on each link."""
    link_flows = {}
    for arc in state.graph.edges():
        total_flow = sum(
            state.flows.get((arc, origin), 0.0)
            for origin in state.origins
        )
        link_flows[arc] = total_flow
    return link_flows


def print_solution_summary(
    state: NetworkState,
    link_flows: Dict[Tuple[int, int], float],
    gap: float
):
    """Print summary of solution."""
    print("\n" + "="*70)
    print("SOLUTION SUMMARY")
    print("="*70)
    print(f"Number of origins: {len(state.origins)}")
    print(f"Number of links with positive flow: {sum(1 for f in link_flows.values() if f > 0.001)}")
    print(f"Total network flow: {sum(link_flows.values()):.2f}")
    print(f"Relative gap: {gap:.6e}")
    print("="*70 + "\n")
