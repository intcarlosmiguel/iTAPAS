"""
Example: Running ADMM-JOR algorithm for DUE-TAP

This script demonstrates how to:
1. Load network topology and OD matrix
2. Configure ADMM-JOR parameters
3. Solve the traffic assignment problem
4. Validate and report results
"""

import sys
import os
from pathlib import Path

# Add paths for imports
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from bibpy.utils import carregar_rede, carregar_viagens, calcular_gap_relativo,gerar_viagens_aleatorias
from bibpy.admm_jor import (
    solve_due_admm_jor,
    extract_link_flows,
    print_solution_summary,
    validate_flow_conservation,
    NetworkState,
    compute_relative_gap
)


def main():
    """Main execution function."""
    
    print("\n" + "="*70)
    print("TRAFFIC ASSIGNMENT USING ADMM-JOR ALGORITHM")
    print("="*70)
    
    # =====================================================================
    # 1. LOAD NETWORK DATA
    # =====================================================================
    
    # Adjust these paths to your network files
    network_file = './fortaleza/edges_fortaleza.txt'
    od_file = './fortaleza/OD2.txt'
    
    print("\n[1] Loading network data...")
    
    try:
        graph = carregar_rede(network_file)
        N = graph.number_of_nodes()
        viagens = gerar_viagens_aleatorias(graph, N*0.02, 1)
    except FileNotFoundError as e:
        print(f"✗ Error loading files: {e}")
        print("\nPlease ensure network and OD matrix files exist.")
        return None
    
    print(f"Network: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} links")
    print(f"OD matrix: {len(viagens)} demand pairs")
    
    # =====================================================================
    # 2. CONFIGURE ALGORITHM PARAMETERS
    # =====================================================================
    
    print("\n[2] Configuring ADMM-JOR parameters...")
    
    params = {
        'rho': 1.5,              # Penalty parameter (sensitivity to network)
        'max_iter': 5000,        # Maximum iterations
        'tolerance': 1e-8,       # Convergence tolerance (relative gap) - mais rigoroso
        'omega': 1.2,            # Initial relaxation factor for JOR
        'adaptive_omega': True   # Adaptive relaxation factor adjustment
    }
    
    print(f"  Penalty parameter (rho): {params['rho']}")
    print(f"  Max iterations: {params['max_iter']}")
    print(f"  Tolerance: {params['tolerance']:.2e}")
    print(f"  Initial relaxation factor: {params['omega']}")
    print(f"  Adaptive omega: {params['adaptive_omega']}")
    
    # =====================================================================
    # 3. SOLVE DUE PROBLEM
    # =====================================================================
    
    print("\n[3] Solving traffic assignment problem...")
    
    flows, final_gap, num_iterations = solve_due_admm_jor(
        graph=graph,
        viagens=viagens,
        **params
    )
    
    # =====================================================================
    # 4. EXTRACT AND VALIDATE SOLUTION
    # =====================================================================
    
    print("\n[4] Post-processing solution...")
    
    # Build state with solution flows
    origens = sorted(set(o for o, d in viagens.keys()))
    destinos = sorted(set(d for o, d in viagens.keys()))

    state = NetworkState(graph, origens, destinos)
    state.flows = flows

    # Extract aggregate link flows from solution state
    link_flows = extract_link_flows(state)

    # Update graph with solution flows
    for arc in graph.edges():
        graph.edges[arc]['fluxo'] = link_flows.get(arc, 0.0)
    
    # Validate solution
    validation_ok = validate_flow_conservation(state, viagens)
    
    # =====================================================================
    # 5. REPORT RESULTS
    # =====================================================================
    
    print_solution_summary(state, link_flows, final_gap)
    
    print("\nTop 10 congested links:")
    print("-" * 70)
    print(f"{'Link':<20} {'Flow':<15} {'Cost':<15} {'v/c Ratio':<15}")
    print("-" * 70)
    
    sorted_links = sorted(
        link_flows.items(),
        key=lambda x: x[1],
        reverse=True
    )
    
    for i, (arc, flow) in enumerate(sorted_links):
        u, v = arc
        capacity = graph.edges[arc]['capacidade']
        cost = graph.edges[arc]['custo']
        vc_ratio = flow / capacity if capacity > 0 else 0.0
        
        print(f"{u:3d} -> {v:<3d} {flow:>14.2f} {cost:>14.4f} {vc_ratio:>14.4f}")
    
    print("-" * 70)
    
    # =====================================================================
    # 6. SUMMARY AND SUCCESS INDICATOR
    # =====================================================================
    
    print("\n" + "="*70)
    print("EXECUTION SUMMARY")
    print("="*70)
    
    status = "SUCCESS" if final_gap < params['tolerance'] else "INCOMPLETE"
    print(f"Status: {status}")
    print(f"Final relative gap: {final_gap:.6e} (tolerance: {params['tolerance']:.2e})")
    print(f"Iterations completed: {num_iterations} / {params['max_iter']}")
    print(f"Flow conservation check: {'PASS' if validation_ok else 'FAIL'}")
    
    print("="*70 + "\n")
    
    return {
        'flows': flows,
        'link_flows': link_flows,
        'gap': final_gap,
        'iterations': num_iterations,
        'graph': graph
    }


if __name__ == "__main__":
    result = main()
    
    if result is None:
        sys.exit(1)
    
    print("✓ ADMM-JOR algorithm completed successfully!")
    sys.exit(0)
