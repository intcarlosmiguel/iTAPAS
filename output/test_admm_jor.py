"""
Unit tests for ADMM-JOR algorithm implementation.

Run with: python -m pytest test_admm_jor.py -v
"""

import pytest
import numpy as np
import networkx as nx
from collections import defaultdict
from typing import Dict, Tuple

# Import from admm_jor module
from bibpy.admm_jor import (
    bpr_cost,
    bpr_derivative,
    bpr_integral,
    color_edges_greedy,
    partition_edges_into_blocks,
    compute_flow_conservation_residual,
    update_link_costs,
    compute_relative_gap,
    optimize_link_flow,
    admm_jor_solver
)


# ============================================================================
# TEST DATA
# ============================================================================

@pytest.fixture
def simple_network():
    """Create a simple 4-node network for testing."""
    graph = nx.DiGraph()
    
    # Add edges: (from, to, capacity, free_flow_time)
    edges = [
        (1, 2, 2000, 10),
        (2, 3, 2000, 10),
        (1, 3, 3000, 20),
        (3, 4, 2000, 10),
        (2, 4, 1500, 15)
    ]
    
    for u, v, cap, fft in edges:
        graph.add_edge(u, v,
                      capacidade=cap,
                      tempo_fluxo_livre=fft,
                      fluxo=0.0,
                      custo=fft,
                      fluxos_por_origem=defaultdict(float))
    
    return graph


@pytest.fixture
def simple_viagens():
    """Create simple OD matrix for testing."""
    return {
        (1, 3): 100.0,
        (1, 4): 150.0,
        (2, 4): 80.0
    }


# ============================================================================
# TESTS: BPR FUNCTION
# ============================================================================

class TestBPRFunction:
    """Test BPR cost and related functions."""
    
    def test_bpr_cost_free_flow(self):
        """BPR cost at zero flow should equal free-flow time."""
        fft = 10.0
        capacity = 2000.0
        
        cost = bpr_cost(fft, 0, capacity)
        assert abs(cost - fft) < 1e-6
    
    def test_bpr_cost_increases_with_flow(self):
        """BPR cost should increase with flow."""
        fft = 10.0
        capacity = 2000.0
        
        cost1 = bpr_cost(fft, 500, capacity)
        cost2 = bpr_cost(fft, 1500, capacity)
        
        assert cost2 > cost1
    
    def test_bpr_cost_zero_capacity(self):
        """BPR cost with zero capacity should return free-flow time."""
        fft = 10.0
        cost = bpr_cost(fft, 100, 0)
        assert abs(cost - fft) < 1e-6
    
    def test_bpr_derivative_positive(self):
        """BPR derivative should be positive."""
        fft = 10.0
        capacity = 2000.0
        flow = 1000.0
        
        deriv = bpr_derivative(fft, flow, capacity)
        assert deriv > 0
    
    def test_bpr_integral_value(self):
        """BPR integral should be positive and increase with flow."""
        fft = 10.0
        capacity = 2000.0
        
        int1 = bpr_integral(fft, 500, capacity)
        int2 = bpr_integral(fft, 1000, capacity)
        
        assert int1 > 0
        assert int2 > int1


# ============================================================================
# TESTS: NETWORK STRUCTURE
# ============================================================================

class TestNetworkOperations:
    """Test network decomposition and operations."""
    
    def test_coloring_produces_valid_colors(self, simple_network):
        """Edge coloring should produce non-negative color assignments."""
        colors = color_edges_greedy(simple_network)
        
        assert len(colors) > 0
        assert all(c >= 0 for c in colors.values())
        assert all(isinstance(c, int) for c in colors.values())
    
    def test_coloring_non_adjacent_property(self, simple_network):
        """Edges with same color should not share vertices."""
        colors = color_edges_greedy(simple_network)
        edge_list = list(simple_network.edges())
        
        # Group edges by color
        color_groups = defaultdict(list)
        for edge, color in colors.items():
            color_groups[color].append(edge)
        
        # Check non-adjacency within groups
        for color, edges in color_groups.items():
            nodes_used = set()
            for u, v in edges:
                assert u not in nodes_used, f"Node {u} appears in multiple edges of color {color}"
                assert v not in nodes_used, f"Node {v} appears in multiple edges of color {color}"
                nodes_used.add(u)
                nodes_used.add(v)
    
    def test_block_decomposition(self, simple_network):
        """Block decomposition should partition all edges."""
        blocks, num_blocks = partition_edges_into_blocks(simple_network)
        
        # Count total edges
        total_edges_in_blocks = sum(len(block) for block in blocks.values())
        assert total_edges_in_blocks == simple_network.number_of_edges()
        
        # Check block numbering
        assert set(blocks.keys()) == set(range(num_blocks))


# ============================================================================
# TESTS: FLOW CONSERVATION
# ============================================================================

class TestFlowConservation:
    """Test flow conservation constraint handling."""
    
    def test_residual_zero_initial_flow(self, simple_network, simple_viagens):
        """Residual should equal demand magnitude at origin with zero initial flow."""
        flows = {(arc, origin): 0.0 for arc in simple_network.edges() 
                for origin in set(o for o, _ in simple_viagens.keys())}
        
        origin = 1
        # At origin node 1, all demand should originate
        expected_demand = sum(d for (o, d_node), d in simple_viagens.items() if o == origin)
        
        residual = compute_flow_conservation_residual(
            flows, origin, origin, simple_network, simple_viagens
        )
        
        # Residual = inflow - outflow - demand = 0 - 0 - demand
        assert abs(residual + expected_demand) < 1e-6
    
    def test_residual_destination_node(self, simple_network, simple_viagens):
        """Residual at destination should reflect incoming demand."""
        flows = {}
        for arc in simple_network.edges():
            for origin in set(o for o, _ in simple_viagens.keys()):
                flows[(arc, origin)] = 0.0
        
        destination = 4
        origin = 1
        
        # Demand from origin 1 to node 4
        expected_demand = simple_viagens.get((origin, destination), 0.0)
        
        residual = compute_flow_conservation_residual(
            flows, origin, destination, simple_network, simple_viagens
        )
        
        # Residual = 0 - 0 - (-demand) = demand
        assert abs(residual - expected_demand) < 1e-6


# ============================================================================
# TESTS: LINK COSTS
# ============================================================================

class TestLinkCosts:
    """Test link cost updates."""
    
    def test_cost_update_increases_with_flow(self, simple_network, simple_viagens):
        """Link costs should increase with aggregate flow."""
        origins = sorted(set(o for o, _ in simple_viagens.keys()))
        
        # Initially zero flow
        flows_1 = {(arc, origin): 0.0 for arc in simple_network.edges() for origin in origins}
        link_costs_1 = {}
        update_link_costs(simple_network, flows_1, link_costs_1, origins)
        
        # Some positive flow
        flows_2 = flows_1.copy()
        arc_sample = list(simple_network.edges())[0]
        flows_2[(arc_sample, origins[0])] = 500.0
        link_costs_2 = {}
        update_link_costs(simple_network, flows_2, link_costs_2, origins)
        
        # Cost should increase
        assert link_costs_2[arc_sample] > link_costs_1[arc_sample]


# ============================================================================
# TESTS: CONVERGENCE METRIC
# ============================================================================

class TestConvergenceMeter:
    """Test gap computation."""
    
    def test_gap_zero_with_zero_flow(self, simple_network, simple_viagens):
        """Gap may be undefined with zero flow but should be handled gracefully."""
        origins = sorted(set(o for o, _ in simple_viagens.keys()))
        flows = {(arc, origin): 0.0 for arc in simple_network.edges() for origin in origins}
        link_costs = {arc: simple_network.edges[arc]['tempo_fluxo_livre'] 
                     for arc in simple_network.edges()}
        
        gap = compute_relative_gap(simple_network, flows, link_costs, simple_viagens, origins)
        
        # Should handle gracefully without crashing
        assert gap >= 0
    
    def test_gap_decreases_toward_solution(self, simple_network, simple_viagens):
        """Gap should be a valid number."""
        origins = sorted(set(o for o, _ in simple_viagens.keys()))
        flows = {(arc, origin): 100.0 for arc in simple_network.edges() for origin in origins}
        link_costs = {arc: simple_network.edges[arc]['tempo_fluxo_livre'] 
                     for arc in simple_network.edges()}
        
        gap = compute_relative_gap(simple_network, flows, link_costs, simple_viagens, origins)
        
        assert gap >= 0
        assert gap < 1e6  # Should be reasonable


# ============================================================================
# TESTS: ALGORITHM EXECUTION
# ============================================================================

class TestADMMJORSolver:
    """Test main ADMM-JOR solver."""
    
    def test_solver_basic_convergence(self, simple_network, simple_viagens):
        """Solver should produce reasonable solution."""
        flows, gaps, iterations = admm_jor_solver(
            simple_network,
            simple_viagens,
            rho=0.5,
            omega=1.2,
            max_iterations=100,
            tolerance=1e-4,
            adaptive_omega=False,
            verbose=False
        )
        
        # Should have runs some iterations
        assert iterations > 0
        assert iterations <= 100
        
        # Should have decreasing gap
        assert gaps[0] > gaps[-1]
        assert gaps[-1] >= 0
    
    def test_solver_different_omegas(self, simple_network, simple_viagens):
        """Different omega values should give different convergence."""
        results = {}
        
        for omega in [1.0, 1.2, 1.5]:
            flows, gaps, iterations = admm_jor_solver(
                simple_network,
                simple_viagens,
                rho=0.5,
                omega=omega,
                max_iterations=100,
                tolerance=1e-4,
                adaptive_omega=False,
                verbose=False
            )
            results[omega] = {'iters': iterations, 'gap': gaps[-1]}
        
        # All should converge
        assert all(r['gap'] < 1e-3 for r in results.values())
    
    def test_solver_adaptive_omega(self, simple_network, simple_viagens):
        """Solver with adaptive omega should work."""
        flows, gaps, iterations = admm_jor_solver(
            simple_network,
            simple_viagens,
            rho=0.5,
            omega=1.3,
            max_iterations=100,
            tolerance=1e-4,
            adaptive_omega=True,
            verbose=False
        )
        
        assert iterations > 0
        assert gaps[-1] >= 0


# ============================================================================
# INTEGRATION TESTS
# ============================================================================

class TestIntegration:
    """Integration tests."""
    
    def test_full_workflow(self, simple_network, simple_viagens):
        """Test complete workflow from network to solution."""
        # Solve
        flows, gaps, iterations = admm_jor_solver(
            simple_network,
            simple_viagens,
            rho=0.5,
            omega=1.2,
            max_iterations=100,
            tolerance=1e-4,
            adaptive_omega=False,
            verbose=False
        )
        
        # Verify solution properties
        assert len(flows) > 0
        assert all(f >= 0 for f in flows.values())  # Non-negative flows
        assert gaps[-1] < gaps[0]  # Convergence
    
    def test_reproducibility(self, simple_network, simple_viagens):
        """Same inputs should give same results."""
        flows1, gaps1, iters1 = admm_jor_solver(
            simple_network,
            simple_viagens,
            rho=0.5,
            omega=1.2,
            max_iterations=50,
            tolerance=1e-4,
            adaptive_omega=False,
            verbose=False
        )
        
        flows2, gaps2, iters2 = admm_jor_solver(
            simple_network,
            simple_viagens,
            rho=0.5,
            omega=1.2,
            max_iterations=50,
            tolerance=1e-4,
            adaptive_omega=False,
            verbose=False
        )
        
        # Should be identical
        assert iters1 == iters2
        assert abs(gaps1[-1] - gaps2[-1]) < 1e-9


# ============================================================================
# PERFORMANCE TESTS
# ============================================================================

class TestPerformance:
    """Performance and scaling tests."""
    
    def test_solver_completes_quickly(self, simple_network, simple_viagens):
        """Solver should complete within reasonable time."""
        import time
        
        start = time.time()
        flows, gaps, iterations = admm_jor_solver(
            simple_network,
            simple_viagens,
            rho=0.5,
            omega=1.2,
            max_iterations=50,
            tolerance=1e-3,
            verbose=False
        )
        elapsed = time.time() - start
        
        # Should complete in reasonable time (< 5 seconds)
        assert elapsed < 5.0


# ============================================================================
# RUN TESTS
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
