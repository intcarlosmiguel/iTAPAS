/**
 * Conjugate Direction Frank-Wolfe (CFW) for the DUE-TAP.
 *
 * Reference:
 *   Mitradjieva, M. & Lindberg, P. O. (2013). The Stiff Is Moving —
 *   Conjugate Direction Frank-Wolfe Methods with Applications to Traffic
 *   Assignment. Transportation Science, 47(2), 280–294.
 *
 * Purely functional style — no global mutable state.
 * Uses igraph for graph storage and Dijkstra, and the BPR helpers
 * already defined in calc.h / define.h.
 *
 * Time  complexity per iteration : O(|O| · SPT_cost + |A|)
 * Space complexity               : O(|A|)   (two extra vectors beyond FW)
 */

#pragma once

#include "calc.h"
#include "define.h"
#include <igraph/igraph.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

/* Safety clamp that prevents alpha_k = 1 */
#define CFW_DELTA 1e-5

/* ================================================================
   1. All-or-Nothing (AON) assignment
   ================================================================
   Builds y^{FW} by routing all demand on shortest paths under the
   supplied edge costs.  Interface mirrors `atualiza_fluxo` from
   leblanc.h.

   Time  complexity: O(|O| · (|A| + |N| log |N|)).
   Space complexity: O(|N|) per origin (inbound vector).
   ================================================================ */
static inline void cfw_atualiza_fluxo(igraph_t *Grafo, struct OD_MATRIX *OD,
                                      igraph_vector_t *fluxo,
                                      igraph_vector_t *pesos) {
  igraph_vector_fill(fluxo, 0.0);
  for (int i = 0; i < OD->size; i++) {
    int fonte = OD->Elementos[i].fonte;

    igraph_vector_int_t inbounds;
    igraph_vector_int_init(&inbounds, 0);
    igraph_get_shortest_paths_dijkstra(Grafo, NULL, NULL, fonte,
                                       igraph_vss_all(), pesos, IGRAPH_OUT,
                                       NULL, &inbounds);

    int n_alvos = igraph_vector_int_size(&OD->Elementos[i].alvos);
    for (int j = 0; j < n_alvos; j++) {
      int alvo = VECTOR(OD->Elementos[i].alvos)[j];
      double volume = VECTOR(OD->Elementos[i].volumes)[j];
      if (VECTOR(inbounds)[alvo] < 0)
        continue;
      while (alvo != fonte) {
        int edge_id = VECTOR(inbounds)[alvo];
        VECTOR(*fluxo)[edge_id] += volume;
        alvo = IGRAPH_FROM(Grafo, edge_id);
      }
    }
    igraph_vector_int_destroy(&inbounds);
  }
}

/* ================================================================
   2. Conjugation weight  alpha_k   (CFW-specific)
   ================================================================
   Computes:
     N_k = sum_a  t'_a * d_residual_a * d_fw_a
     D_k = sum_a  t'_a * d_residual_a * (d_fw_a - d_residual_a)

   Clamping rules:
     D_k != 0  and  N_k/D_k in [0, 1-delta]  →  alpha = N_k/D_k
     D_k != 0  and  N_k/D_k > 1-delta         →  alpha = 1-delta
     otherwise                                  →  alpha = 0

   Time  complexity: O(|A|).
   Space complexity: O(1)  (scalar accumulation).
   ================================================================ */
static inline double cfw_calcular_alfa(struct PARAMETERS *P,
                                       igraph_vector_t *x,
                                       igraph_vector_t *d_fw,
                                       igraph_vector_t *d_residual) {
  double N_k = 0.0, D_k = 0.0;

  for (int a = 0; a < P->L; a++) {
    double t_prime = single_BPR_derivate(VECTOR(*x)[a],
                                         VECTOR(P->cost_time)[a],
                                         VECTOR(P->capacidade)[a]);
    double hd = t_prime * VECTOR(*d_residual)[a];
    N_k += hd * VECTOR(*d_fw)[a];
    D_k += hd * (VECTOR(*d_fw)[a] - VECTOR(*d_residual)[a]);
  }

  if (fabs(D_k) < 1e-15)
    return 0.0;

  double razao = N_k / D_k;
  if (razao < 0.0)
    return 0.0;
  if (razao > 1.0 - CFW_DELTA)
    return 1.0 - CFW_DELTA;
  return razao;
}

/* ================================================================
   3. Line search  (Newton step with clamping to [0, 1])
   ================================================================
   tau_k = - sum_a t_a(x_k) d_k_a  /  sum_a t'_a(x_k) d_k_a^2

   Time  complexity: O(|A|).
   Space complexity: O(1).
   ================================================================ */
static inline double cfw_line_search(struct PARAMETERS *P,
                                     igraph_vector_t *x,
                                     igraph_vector_t *d_k) {
  double num = 0.0, den = 0.0;

  for (int a = 0; a < P->L; a++) {
    double t_a, t_prime;
    single_BPR_fused(VECTOR(*x)[a], VECTOR(P->cost_time)[a],
                     VECTOR(P->capacidade)[a], &t_a, &t_prime);
    double da = VECTOR(*d_k)[a];
    num -= t_a * da;
    den += t_prime * da * da;
  }

  if (den <= 0.0)
    return (num > 0.0) ? 1.0 : 0.0;

  double tau = num / den;
  if (tau < 0.0)
    tau = 0.0;
  if (tau > 1.0)
    tau = 1.0;
  return tau;
}

/* ================================================================
   4. Bounds and relative gap
   ================================================================
   UBD = T(x_{k+1})                       (Beckmann at new point)
   LBD = T(x_k) + sum t_a(x_k)(y_fw - x_k)  (linearisation bound)
   BLB = max(BLB_{k-1}, LBD)
   RE  = (UBD - BLB) / |BLB|

   Time  complexity: O(|A|).
   Space complexity: O(1).
   ================================================================ */
static inline double cfw_beckmann(struct PARAMETERS *P,
                                  igraph_vector_t *x) {
  double obj = 0.0;
  for (int a = 0; a < P->L; a++) {
    double f = VECTOR(*x)[a];
    double t0 = VECTOR(P->cost_time)[a];
    double c = VECTOR(P->capacidade)[a];
    if (c > 0.0) {
      obj += t0 * f + t0 * ALPHA * pow(f, BETA + 1.0) /
                           ((BETA + 1.0) * pow(c, BETA));
    } else {
      obj += t0 * f;
    }
  }
  return obj;
}

static inline double cfw_calcular_bounds(struct PARAMETERS *P,
                                         igraph_vector_t *x_k,
                                         igraph_vector_t *x_kp1,
                                         igraph_vector_t *y_fw,
                                         igraph_vector_t *custos_k,
                                         double blb_anterior,
                                         double *out_blb) {
  double T_xk = cfw_beckmann(P, x_k);
  double UBD = cfw_beckmann(P, x_kp1);

  double lin = 0.0;
  for (int a = 0; a < P->L; a++)
    lin += VECTOR(*custos_k)[a] * (VECTOR(*y_fw)[a] - VECTOR(*x_k)[a]);

  double LBD = T_xk + lin;
  double BLB = (LBD > blb_anterior) ? LBD : blb_anterior;

  *out_blb = BLB;

  if (fabs(BLB) < 1e-15)
    return 1.0 / 0.0; /* inf */

  return (UBD - BLB) / fabs(BLB);
}

/* ================================================================
   5. Main solver:  CFW
   ================================================================
   Interface mirrors leblanc():
     - BPR_PARAMETERS : network parameters (capacidade, cost_time, L, N)
     - OD             : origin-destination demand
     - Grafo          : igraph directed graph
     - solucao        : output flow vector (allocated by caller or here)
     - WARM_START     : if true, solucao already contains a feasible flow
     - saved_inbounds : array of SPT predecessor vectors (one per origin),
                        updated at convergence for future warm starts.
     - saved_sp_costs : flat array of shortest-path costs per OD pair,
                        updated at convergence.

   Time  complexity per iteration : O(|O| · SPT + |A|).
   Space complexity               : O(|A|) extra.
   ================================================================ */
static inline void cfw_solver(struct PARAMETERS *BPR_PARAMETERS,
                              struct OD_MATRIX *OD, igraph_t *Grafo,
                              igraph_vector_t *solucao, bool warm_start,
                              igraph_vector_int_t *saved_inbounds,
                              double *saved_sp_costs) {
  int L = BPR_PARAMETERS->L;

  /* ── Early exit for warm start already converged ──────── */
  if (warm_start) {
    double GAP_approx = relative_gap_approximate(solucao, saved_sp_costs,
                                                  BPR_PARAMETERS, OD);
    if (GAP_approx < EPSILON) {
      printf("CFW warm start already converged (approx GAP=%e). Skipping.\n",
             GAP_approx);
      return;
    }
  }

  /* ── Scratch vectors ────────────────────────────────────── */
  igraph_vector_t tempo;        /* t_a(x_k)          */
  igraph_vector_t y_fw;         /* AON solution       */
  igraph_vector_t d_fw;         /* y_fw - x           */
  igraph_vector_t s_anterior;   /* s_{k-1}            */
  igraph_vector_t s_k;          /* current target     */
  igraph_vector_t d_k;          /* s_k - x            */
  igraph_vector_t d_residual;   /* s_{k-1} - x_k      */
  igraph_vector_t x_novo;       /* x_{k+1}            */

  igraph_vector_init(&tempo, L);
  igraph_vector_init(&y_fw, L);
  igraph_vector_init(&d_fw, L);
  igraph_vector_init(&s_anterior, L);
  igraph_vector_init(&s_k, L);
  igraph_vector_init(&d_k, L);
  igraph_vector_init(&d_residual, L);
  igraph_vector_init(&x_novo, L);

  /* ── Initialisation ─────────────────────────────────────── */
  if (!warm_start) {
    igraph_vector_init(solucao, L);
    cfw_atualiza_fluxo(Grafo, OD, solucao, &BPR_PARAMETERS->cost_time);
  }

  double tau_anterior = 1.0;
  bool s_anterior_valido = false;
  double BLB = -1e300;
  double RE = 1.0;
  int iteracoes = 0;

  /* ── Main loop ──────────────────────────────────────────── */
  for (int k = 0; k < MAXIMO_ITERACOES; k++) {

    /* (a) Current costs + AON direction */
    BPR(&tempo, BPR_PARAMETERS, solucao);
    cfw_atualiza_fluxo(Grafo, OD, &y_fw, &tempo);

    for (int a = 0; a < L; a++)
      VECTOR(d_fw)[a] = VECTOR(y_fw)[a] - VECTOR(*solucao)[a];

    /* (b) Conjugation */
    if (k == 0 || tau_anterior >= 1.0 || !s_anterior_valido) {
      /* Fall back to pure FW */
      igraph_vector_update(&s_k, &y_fw);
    } else {
      /* d_residual = s_{k-1} - x_k */
      for (int a = 0; a < L; a++)
        VECTOR(d_residual)[a] =
            VECTOR(s_anterior)[a] - VECTOR(*solucao)[a];

      double alfa = cfw_calcular_alfa(BPR_PARAMETERS, solucao,
                                      &d_fw, &d_residual);

      for (int a = 0; a < L; a++)
        VECTOR(s_k)[a] = alfa * VECTOR(s_anterior)[a] +
                         (1.0 - alfa) * VECTOR(y_fw)[a];
    }

    /* (c) Line search */
    for (int a = 0; a < L; a++)
      VECTOR(d_k)[a] = VECTOR(s_k)[a] - VECTOR(*solucao)[a];

    double tau = cfw_line_search(BPR_PARAMETERS, solucao, &d_k);

    /* (d) Update */
    for (int a = 0; a < L; a++)
      VECTOR(x_novo)[a] = VECTOR(*solucao)[a] + tau * VECTOR(d_k)[a];

    RE = cfw_calcular_bounds(BPR_PARAMETERS, solucao, &x_novo,
                             &y_fw, &tempo, BLB, &BLB);
    iteracoes = k + 1;

    if (k < 10 || k % 10 == 0 || RE < EPSILON)
      printf("  CFW Iter %5d | tau = %.6f | RE = %.6e\n", k, tau, RE);

    /* Prepare next iteration */
    igraph_vector_update(&s_anterior, &s_k);
    s_anterior_valido = true;
    tau_anterior = tau;
    igraph_vector_update(solucao, &x_novo);

    if (RE < EPSILON) {
      printf("CFW converged in %d iterations (RE = %e)\n", iteracoes, RE);
      break;
    }
  }

  if (RE >= EPSILON)
    printf("CFW: max iterations reached (%d). RE = %e\n", iteracoes, RE);

  /* ── Save SPT for future warm starts ─────────────────── */
  if (saved_inbounds != NULL || saved_sp_costs != NULL) {
    BPR(&tempo, BPR_PARAMETERS, solucao);
    relative_gap(solucao, Grafo, BPR_PARAMETERS, OD,
                 saved_inbounds, saved_sp_costs);
  }

  igraph_vector_destroy(&tempo);
  igraph_vector_destroy(&y_fw);
  igraph_vector_destroy(&d_fw);
  igraph_vector_destroy(&s_anterior);
  igraph_vector_destroy(&s_k);
  igraph_vector_destroy(&d_k);
  igraph_vector_destroy(&d_residual);
  igraph_vector_destroy(&x_novo);
}
