"""
Conjugate Direction Frank-Wolfe (CFW) for the Deterministic User Equilibrium
Traffic Assignment Problem (DUE-TAP).

Reference:
    Mitradjieva, M. & Lindberg, P. O. (2013). The Stiff Is Moving —
    Conjugate Direction Frank-Wolfe Methods with Applications to Traffic
    Assignment. Transportation Science, 47(2), 280–294.

This module is written in a purely functional style: every routine receives
its data as arguments and returns results without mutating shared state.
Network topology is stored as a ``networkx.DiGraph`` loaded via
``bibpy.utils.carregar_rede``; edge attributes used are:

    - ``tempo_fluxo_livre``  (float) : free-flow travel time  t₀_a
    - ``capacidade``         (float) : link capacity           c_a

Time  complexity per iteration : O(|O| · SPT_cost + |A|)
Space complexity               : O(|A|)  (two extra vectors beyond FW)
"""

import numpy as np
import networkx as nx
from typing import Dict, Tuple, List, Optional
import time

# ============================================================================
# Constants
# ============================================================================

BPR_ALPHA: float = 0.15
BPR_BETA: float = 4.0
DELTA: float = 1e-5          # safety clamp for alpha_k
EPSILON: float = 1e-4        # default relative-gap tolerance
MAX_ITER: int = 1000


# ============================================================================
# 1. BPR cost functions  –  t_a(x_a), t'_a(x_a), ∫₀^{x_a} t_a(s) ds
# ============================================================================

def custo_bpr(t0: float, fluxo: float, capacidade: float) -> float:
    r"""
    Bureau of Public Roads travel-time function.

    .. math::
        t_a(x_a) = t_0 \bigl(1 + \alpha\,(x_a / c_a)^{\beta}\bigr)

    Time  complexity: O(1).
    Space complexity: O(1).
    """
    if capacidade <= 0.0:
        return t0
    return t0 * (1.0 + BPR_ALPHA * (fluxo / capacidade) ** BPR_BETA)


def derivada_bpr(t0: float, fluxo: float, capacidade: float) -> float:
    r"""
    Derivative of the BPR function w.r.t. flow.

    .. math::
        t'_a(x_a) = \frac{t_0\,\alpha\,\beta}{c_a}
                     \left(\frac{x_a}{c_a}\right)^{\beta - 1}

    Time  complexity: O(1).
    Space complexity: O(1).
    """
    if capacidade <= 0.0:
        return 0.0
    return (t0 * BPR_ALPHA * BPR_BETA
            * (fluxo / capacidade) ** (BPR_BETA - 1.0) / capacidade)


def integral_bpr(t0: float, fluxo: float, capacidade: float) -> float:
    r"""
    Beckmann objective component for one link.

    .. math::
        \int_0^{x_a} t_a(s)\,ds
        = t_0\,x_a + \frac{t_0\,\alpha}{\beta+1}
          \frac{x_a^{\,\beta+1}}{c_a^{\,\beta}}

    Time  complexity: O(1).
    Space complexity: O(1).
    """
    if capacidade <= 0.0:
        return t0 * fluxo
    return (t0 * fluxo
            + t0 * BPR_ALPHA / (BPR_BETA + 1.0)
            * fluxo ** (BPR_BETA + 1.0) / capacidade ** BPR_BETA)


# ============================================================================
# 2. Vectorised helpers over the arc set
# ============================================================================

def _extrair_atributos_arcos(grafo: nx.DiGraph) -> Tuple[np.ndarray,
                                                          np.ndarray,
                                                          list]:
    """
    Returns ordered arrays ``(t0, capacidade)`` and the list of arc keys
    matching that order.  Called once; the ordering is reused throughout.

    Time  complexity: O(|A|).
    Space complexity: O(|A|).
    """
    arcos = list(grafo.edges())
    n = len(arcos)
    t0 = np.empty(n)
    cap = np.empty(n)
    for i, (u, v) in enumerate(arcos):
        dados = grafo[u][v]
        t0[i] = dados['tempo_fluxo_livre']
        cap[i] = dados['capacidade']
    return t0, cap, arcos


def calcular_custos(t0: np.ndarray, x: np.ndarray,
                    cap: np.ndarray) -> np.ndarray:
    """
    Vectorised BPR cost for all arcs.

    Time  complexity: O(|A|).
    Space complexity: O(|A|) for the returned array.
    """
    ratio = np.where(cap > 0.0, x / cap, 0.0)
    return t0 * (1.0 + BPR_ALPHA * ratio ** BPR_BETA)


def calcular_derivadas(t0: np.ndarray, x: np.ndarray,
                       cap: np.ndarray) -> np.ndarray:
    """
    Vectorised BPR derivative for all arcs  (diagonal of the Hessian).

    Time  complexity: O(|A|).
    Space complexity: O(|A|).
    """
    seguro = np.where(cap > 0.0, cap, 1.0)
    ratio = x / seguro
    derivs = t0 * BPR_ALPHA * BPR_BETA * ratio ** (BPR_BETA - 1.0) / seguro
    return np.where(cap > 0.0, derivs, 0.0)


def calcular_objetivo(t0: np.ndarray, x: np.ndarray,
                      cap: np.ndarray) -> float:
    r"""
    Beckmann objective  :math:`T(\mathbf{x}) = \sum_a \int_0^{x_a} t_a(s)\,ds`.

    Time  complexity: O(|A|).
    Space complexity: O(|A|) intermediate.
    """
    seguro = np.where(cap > 0.0, cap, 1.0)
    integrais = (t0 * x
                 + t0 * BPR_ALPHA / (BPR_BETA + 1.0)
                 * x ** (BPR_BETA + 1.0) / seguro ** BPR_BETA)
    # For arcs with cap <= 0, the integral is simply t0 * x.
    integrais = np.where(cap > 0.0, integrais, t0 * x)
    return float(np.sum(integrais))


# ============================================================================
# 3. All-or-Nothing (AON) assignment
# ============================================================================

def atribuicao_tudo_ou_nada(grafo: nx.DiGraph,
                            viagens: Dict[Tuple[int, int], float],
                            arcos: list,
                            custos_arco: np.ndarray) -> np.ndarray:
    r"""
    Compute the AON assignment  :math:`\mathbf{y}^{\text{FW}}`.

    For each origin, a shortest-path tree is built under the current costs;
    all demand is loaded onto those shortest paths.

    Parameters
    ----------
    grafo : nx.DiGraph
        Network topology (nodes and edges only — costs are overwritten).
    viagens : dict
        OD demand  ``{(o, d): volume, ...}``.
    arcos : list
        Ordered list of arcs matching the vectorised arrays.
    custos_arco : np.ndarray
        Current link costs, shape ``(|A|,)``.

    Returns
    -------
    y : np.ndarray
        AON flow vector aligned with ``arcos``, shape ``(|A|,)``.

    Time  complexity: O(|O| · (|A| + |N| log |N|)).
    Space complexity: O(|A| + |N|).
    """
    # Write current costs into the graph for Dijkstra
    for i, (u, v) in enumerate(arcos):
        grafo[u][v]['custo'] = float(custos_arco[i])

    # Index: arc -> position  (built once per call, fast dict lookup)
    indice_arco = {a: i for i, a in enumerate(arcos)}
    n = len(arcos)
    y = np.zeros(n)

    origens = sorted({o for o, _ in viagens})

    for origem in origens:
        try:
            pred, dist = nx.dijkstra_predecessor_and_distance(
                grafo, source=origem, weight='custo'
            )
        except nx.NodeNotFound:
            continue

        for (o, d), demanda in viagens.items():
            if o != origem or demanda <= 0.0:
                continue
            if d not in dist:
                continue

            # Reconstruct shortest path
            caminho = [d]
            atual = d
            while atual != origem:
                if atual not in pred or len(pred[atual]) == 0:
                    caminho = []
                    break
                atual = pred[atual][0]
                caminho.append(atual)

            if not caminho:
                continue

            caminho.reverse()
            for j in range(len(caminho) - 1):
                arco = (caminho[j], caminho[j + 1])
                idx = indice_arco.get(arco)
                if idx is not None:
                    y[idx] += demanda

    return y


# ============================================================================
# 4. Line search  (Newton step with clamping)
# ============================================================================

def line_search_newton(custos: np.ndarray, derivadas: np.ndarray,
                       d: np.ndarray) -> float:
    r"""
    Single Newton step for the one-dimensional line search

    .. math::
        \tau_k = -\frac{\sum_a t_a(x_{k,a})\,d_{k,a}}
                       {\sum_a t'_a(x_{k,a})\,d_{k,a}^2},
        \quad \tau_k \leftarrow \max(0,\,\min(1,\,\tau_k)).

    Parameters
    ----------
    custos : np.ndarray
        :math:`t_a(x_{k,a})` for each arc.
    derivadas : np.ndarray
        :math:`t'_a(x_{k,a})` for each arc (diagonal Hessian).
    d : np.ndarray
        Search direction  :math:`\mathbf{d}_k = \mathbf{s}_k - \mathbf{x}_k`.

    Returns
    -------
    tau : float
        Step size clamped to [0, 1].

    Time  complexity: O(|A|).
    Space complexity: O(|A|) intermediate.
    """
    numerador = -float(np.dot(custos, d))
    denominador = float(np.dot(derivadas, d ** 2))
    if denominador <= 0.0:
        return 1.0 if numerador > 0.0 else 0.0
    tau = numerador / denominador
    return float(np.clip(tau, 0.0, 1.0))


# ============================================================================
# 5. Conjugation step  (CFW-specific)
# ============================================================================

def calcular_alfa_cfw(derivadas: np.ndarray,
                      d_fw: np.ndarray,
                      d_residual: np.ndarray,
                      delta: float = DELTA) -> float:
    r"""
    Compute the conjugation weight  :math:`\alpha_k` for CFW.

    .. math::
        N_k = \sum_a t'_a\;\tilde{d}_{k-1,a}\;d^{\text{FW}}_{k,a},\qquad
        D_k = \sum_a t'_a\;\tilde{d}_{k-1,a}\;
              \bigl(d^{\text{FW}}_{k,a} - \tilde{d}_{k-1,a}\bigr).

    Clamping rules:

    * :math:`D_k \neq 0` and :math:`N_k/D_k \in [0, 1-\delta]`:
      :math:`\alpha_k = N_k / D_k`.
    * :math:`D_k \neq 0` and :math:`N_k/D_k > 1-\delta`:
      :math:`\alpha_k = 1 - \delta`.
    * Otherwise (denominator zero or ratio negative):
      :math:`\alpha_k = 0` (fall back to FW).

    Time  complexity: O(|A|).
    Space complexity: O(|A|) intermediate.
    """
    # Hessian-weighted inner products
    hd_res = derivadas * d_residual           # H_k * d̃_{k-1}
    N_k = float(np.dot(hd_res, d_fw))
    D_k = float(np.dot(hd_res, d_fw - d_residual))

    if abs(D_k) < 1e-15:
        return 0.0

    razao = N_k / D_k
    if razao < 0.0:
        return 0.0
    if razao > 1.0 - delta:
        return 1.0 - delta
    return razao


# ============================================================================
# 6. Bounds and gap
# ============================================================================

def calcular_bounds(t0: np.ndarray, cap: np.ndarray,
                    x_k: np.ndarray, x_kp1: np.ndarray,
                    y_fw: np.ndarray,
                    custos_k: np.ndarray,
                    blb_anterior: float
                    ) -> Tuple[float, float, float, float]:
    r"""
    Compute upper bound, lower bound, best lower bound, and relative gap.

    .. math::
        \text{UBD}_k &= T(\mathbf{x}_{k+1}),\\
        \text{LBD}_k &= T(\mathbf{x}_k)
                        + \sum_a t_a(x_{k,a})\,(y^{\text{FW}}_{k,a} - x_{k,a}),\\
        \text{BLB}_k &= \max(\text{BLB}_{k-1},\;\text{LBD}_k),\\
        \text{RE}_k  &= (\text{UBD}_k - \text{BLB}_k)\,/\,|\text{BLB}_k|.

    Time  complexity: O(|A|).
    Space complexity: O(|A|) intermediate.
    """
    T_xk = calcular_objetivo(t0, x_k, cap)
    UBD = calcular_objetivo(t0, x_kp1, cap)
    LBD = T_xk + float(np.dot(custos_k, y_fw - x_k))
    BLB = max(blb_anterior, LBD)

    if abs(BLB) < 1e-15:
        RE = float('inf')
    else:
        RE = (UBD - BLB) / abs(BLB)

    return UBD, LBD, BLB, RE


# ============================================================================
# 7. Main solver: CFW
# ============================================================================

def resolver_cfw(grafo: nx.DiGraph,
                 viagens: Dict[Tuple[int, int], float],
                 max_iter: int = MAX_ITER,
                 epsilon: float = EPSILON,
                 delta: float = DELTA,
                 verbose: bool = True
                 ) -> Dict:
    r"""
    Solve the DUE-TAP via the Conjugate Direction Frank-Wolfe method.

    The algorithm is detailed in §3 of Mitradjieva & Lindberg (2013).  The
    main loop computes, at each iteration *k*:

    1. **AON direction** :math:`\mathbf{d}_k^{\text{FW}} = \mathbf{y}_k^{\text{FW}} - \mathbf{x}_k`,
    2. **Conjugation** to obtain the target point :math:`\mathbf{s}_k`,
    3. **Newton line search** along :math:`\mathbf{x}_k + \tau\,(\mathbf{s}_k - \mathbf{x}_k)`,
    4. **Update** :math:`\mathbf{x}_{k+1}` and evaluate convergence bounds.

    Parameters
    ----------
    grafo : nx.DiGraph
        Network loaded by :func:`bibpy.utils.carregar_rede`.
    viagens : dict
        OD demand matrix ``{(o, d): volume}``.
    max_iter : int
        Maximum number of iterations.
    epsilon : float
        Relative-gap convergence tolerance.
    delta : float
        Safety parameter to prevent :math:`\alpha_k = 1`.
    verbose : bool
        If ``True``, print iteration log.

    Returns
    -------
    resultado : dict
        Keys:

        - ``'fluxo'``       : np.ndarray, equilibrium link flows.
        - ``'custo'``       : np.ndarray, final link costs.
        - ``'arcos'``       : list, ordered arc keys.
        - ``'gap_relativo'``: float, final relative gap.
        - ``'iteracoes'``   : int, iterations executed.
        - ``'objetivo'``    : float, Beckmann objective at solution.
        - ``'historico_gap'``: list[float], gap per iteration.
        - ``'tempo'``       : float, wall-clock time in seconds.

    Time  complexity per iteration: O(|O| · SPT + |A|).
    Space complexity: O(|A|).
    """
    inicio = time.time()

    # ------------------------------------------------------------------
    # Pre-process: extract ordered arc attributes
    # ------------------------------------------------------------------
    t0, cap, arcos = _extrair_atributos_arcos(grafo)
    n_arcos = len(arcos)

    if verbose:
        origens_set = {o for o, _ in viagens}
        print("\n" + "=" * 70)
        print(" CONJUGATE FRANK-WOLFE  (CFW)")
        print("=" * 70)
        print(f"  Nós       : {grafo.number_of_nodes()}")
        print(f"  Arcos     : {n_arcos}")
        print(f"  Origens   : {len(origens_set)}")
        print(f"  Pares OD  : {len(viagens)}")
        print(f"  ε (tol)   : {epsilon:.1e}")
        print(f"  max iter  : {max_iter}")
        print("-" * 70)

    # ------------------------------------------------------------------
    # Initialisation: AON under free-flow costs  → x_0
    # ------------------------------------------------------------------
    custos_ff = calcular_custos(t0, np.zeros(n_arcos), cap)
    print(f"Primeira atribuição de fluxo...")
    x = atribuicao_tudo_ou_nada(grafo, viagens, arcos, custos_ff)

    # State for conjugation
    s_anterior: Optional[np.ndarray] = None   # s_{k-1}
    tau_anterior: float = 1.0                 # tau_{k-1}

    BLB = -np.inf
    historico_gap: List[float] = []

    # ------------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------------
    for k in range(max_iter):

        # (a) Compute current costs and AON direction
        custos_k = calcular_custos(t0, x, cap)
        print(f"Atribuição de fluxo...")    
        y_fw = atribuicao_tudo_ou_nada(grafo, viagens, arcos, custos_k)
        d_fw = y_fw - x

        # (b) Conjugation
        if k == 0 or tau_anterior >= 1.0 or s_anterior is None:
            # No useful residual direction — fall back to pure FW
            s_k = y_fw.copy()
        else:
            # Residual direction: d̃_{k-1} = s_{k-1} - x_k
            d_residual = s_anterior - x
            derivadas_k = calcular_derivadas(t0, x, cap)
            alfa = calcular_alfa_cfw(derivadas_k, d_fw, d_residual, delta)
            s_k = alfa * s_anterior + (1.0 - alfa) * y_fw

        # (c) Line search (Newton step on τ ∈ [0,1])
        d_k = s_k - x
        derivadas_ls = calcular_derivadas(t0, x, cap)
        tau = line_search_newton(custos_k, derivadas_ls, d_k)

        # (d) Update and convergence test
        x_novo = x + tau * d_k

        UBD, LBD, BLB, RE = calcular_bounds(
            t0, cap, x, x_novo, y_fw, custos_k, BLB
        )
        historico_gap.append(RE)

        if verbose and (k < 10 or k % 10 == 0 or RE < epsilon):
            print(f"  Iter {k:5d} | τ = {tau:.6f} | "
                  f"RE = {RE:.6e} | T(x) = {UBD:.6e}")

        # Prepare next iteration
        s_anterior = s_k
        tau_anterior = tau
        x = x_novo

        if RE < epsilon:
            if verbose:
                print(f"\n✓ Convergência atingida na iteração {k}  "
                      f"(RE = {RE:.4e} < ε = {epsilon:.1e})")
            break
        else:
            if verbose:
                print(f"\n⚠ Máximo de iterações atingido ({max_iter}).  "
                    f"RE = {RE:.4e}")

    elapsed = time.time() - inicio

    # ------------------------------------------------------------------
    # Write final flows and costs back into the graph
    # ------------------------------------------------------------------
    custos_final = calcular_custos(t0, x, cap)
    for i, (u, v) in enumerate(arcos):
        grafo[u][v]['fluxo'] = float(x[i])
        grafo[u][v]['custo'] = float(custos_final[i])

    if verbose:
        obj = calcular_objetivo(t0, x, cap)
        print(f"\n  Objetivo final (Beckmann): {obj:.6e}")
        print(f"  Tempo de execução        : {elapsed:.2f} s")
        print("=" * 70 + "\n")

    return {
        'fluxo': x,
        'custo': custos_final,
        'arcos': arcos,
        'gap_relativo': RE,
        'iteracoes': k + 1,
        'objetivo': calcular_objetivo(t0, x, cap),
        'historico_gap': historico_gap,
        'tempo': elapsed,
    }
