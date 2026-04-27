#pragma once

#include <float.h> // Para DBL_MAX
#include <igraph/igraph.h>
#include <limits.h> // For INT_MAX
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ALPHA 0.15
#define BETA 4
const double decresse = 0.9, min_step = 1E-8, max_step = 1.0, ftol = 0.2,
             X_TOLERANCIA = 0.0001;
const int MAXIMO_ITERACOES = 10000;
const int THREADS = 10;
double EPSILON = 1e-1;
#define TAPAS_FLOW_TOL 1e-10
#define TAPAS_COST_TOL 1e-10

// Constants from iTAPAS Paper (Xie et al.)
#define TAPAS_EPSILON 1e-1 // epsilon: Flow precision for potential links
#define TAPAS_THETA 1e-16  // theta: Cost precision for potential links
#define TAPAS_NU 0.25      // nu: Flow effective factor (min flow on s2)
#define TAPAS_MU 0.5       // mu: Cost effective factor (min cost diff)

// Algorithm B Optimization Constants
#define INNER_ITERATIONS 20 // Inner iterations per bush
#define MIN_LINK_FLOW 1e-14 // Minimum flow to consider arc active

/**
 * MergeNode: Representa um nó com múltiplos arcos de entrada na bush.
 * Armazena informações para balanceamento de fluxo otimizado.
 */
struct MergeNode {
  int node_id;          // ID do nó no grafo original
  int *approaches;      // IDs dos arcos de entrada na bush (grafo original)
  double *approachFlow; // Fluxo em cada arco de entrada
  int numApproaches;    // Número de arcos de entrada
  int capacity;         // Capacidade alocada para approaches[]
  int SPlink; // Índice do arco mais curto em approaches[] (-1 se inválido)
  int LPlink; // Índice do arco mais longo com fluxo (-1 se inválido)
  int divergenceNode; // LCA entre SP e LP (cached, -1 se não calculado)
};

struct PARAMETERS {
  igraph_vector_t capacidade;
  igraph_vector_t cost_time;
  int L;
  int N;
};
struct ElementOD {
  int fonte;
  igraph_vector_int_t alvos;
  igraph_vector_t volumes;
  igraph_vector_t warm_volumes;
};
struct OD_MATRIX {
  struct ElementOD *Elementos;
  int size;
  int n_elements;
};
struct min_max_bush {
  igraph_vector_int_t min_edges;
  igraph_vector_int_t max_edges;
  igraph_vector_t dist_shortest_local;
  igraph_vector_t dist_longest_local;
};

struct BUSH {
  bool *is_ingraph; // Indica se o nó está no grafo
  igraph_t Grafo;   // Grafo da bush
  int n_alvos;
  struct min_max_bush paths;
  igraph_vector_t flow_per_origin;       // Fluxo por alvo
  igraph_vector_int_t topological_order; // Ordem topológica dos nós
  igraph_vector_int_t topological_value;
  igraph_vector_int_t edge_id;

  // Merge Nodes (otimização Algorithm B)
  struct MergeNode *merges; // Array de merge nodes
  int numMerges;            // Quantidade de merge nodes
  int *nodeToMerge; // Mapeamento: nodeToMerge[node_id] = índice em merges[] (-1
                    // se não é merge)
};

struct CONJUNTO_SOLUCAO {
  igraph_vector_t flow;
  igraph_vector_t time;
};

struct SPT {
  igraph_vector_int_t antecessores;
  igraph_vector_t dist; // distâncias do menor caminho (análogo a dist_shortest_local)
};

struct WARM_START {
  bool has_warm_start;
  struct SPT *SPT;
  igraph_t grafo_bush;
  double INCREMENT;
  igraph_vector_t flow;
  igraph_vector_t time;
  struct PAS *conjunto_pas;
  int num_pas;
};

/* ========================================================================
   GREEDY PATH-BASED ALGORITHM STRUCTURES (Xie et al., 2018)
   ======================================================================== */

#define GREEDY_MAX_INNER_ITER 1000   /* MaxI: Limite do loop interno */
#define GREEDY_UPDATE_DELTA_FREQ 100 /* Frequência para atualizar Delta_rs */

/* Estrutura para representar um caminho */
struct CAMINHO {
  igraph_vector_int_t arcos; /* IDs dos arcos em ordem topológica */
  double fluxo;              /* Fluxo no caminho f_h */
  double custo_v;            /* Custo do caminho v_h */
  double derivada_s;         /* Segunda derivada s_h */
  double termo_c;            /* Termo linearizado c_h = v_h - s_h * g_h */
};

/* Conjunto de caminhos para cada par OD */
struct CONJUNTO_CAMINHOS {
  struct CAMINHO *caminhos; /* Vetor dinâmico de caminhos */
  int num_caminhos;         /* Quantidade atual */
  int capacidade;           /* Capacidade alocada */
  double delta_rs;          /* Discrepância max-min de custo */
};

/* Estado global do algoritmo Greedy */
/* Estado global do algoritmo Greedy (V2 - Indexado por par OD) */
struct GREEDY_STATE {
  struct CONJUNTO_CAMINHOS *conjuntos; /* Array de conjuntos, um por par OD */
  int *origem_idx;    /* origem_idx[i] = índice da origem no OD_MATRIX */
  int *destino_idx;   /* destino_idx[i] = índice do destino na lista de alvos */
  double *demandas;   /* demandas[i] = demanda d_rs */
  int total_pares;    /* Número total de pares OD */
  double RG_anterior; /* Gap relativo da iteração anterior */
};

#define ITAPAS_MAX_ITER 5000
#define TAPAS_TOLERANCE 1e-10

struct PAS_Origin {
  int origem;     // ID da Origem
  double flow_s1; // Fluxo no segmento 1 (min flow)
  double flow_s2; // Fluxo no segmento 2 (min flow)
  double shift;   // Deslocamento calculado
};

struct PAS {
  igraph_vector_int_t c1; // Segmento 1 (Árvore) - ordem topológica

  igraph_vector_int_t c2; // Segmento 2 (Atalho) - ordem topológica

  struct PAS_Origin *origins; // Lista de origens que compartilham este PAS
  int num_origins;
  int capacity_origins;

  // Assinaturas para verificação de unicidade (conjuntos ordenados)
  igraph_vector_int_t s1_sorted;
  igraph_vector_int_t s2_sorted;
};

typedef struct {

    igraph_vector_int_t index;
    igraph_vector_int_t low;
    igraph_vector_bool_t onstack;
    igraph_vector_int_t stack;
    igraph_vector_int_t comp;

    const igraph_vector_bool_t *use_edge;

    int index_counter;
    int comp_counter;

} tarjan_ctx_t;

typedef struct {
  /* Edge-level state arrays (indexed by edge ID) */
  bool *congested_edges; /* true if edge exceeds congestion threshold */

  /* Node-level state arrays (indexed by node ID) */
  bool *nodes_in_giant_component; /* true if node belongs to largest SCC */

  /* Cluster identifiers */
  int giant_component_id; /* ID of the largest cluster */

  /* Normalized cluster size metrics (fraction of total edges) */
  double giant_component_fraction;  /* Size of GCC / total edges */
  double second_component_fraction; /* Size of second largest / total edges */

  /* Statistical metrics for finite clusters (excluding GCC) */
  double finite_cluster_mean; /* Arithmetic mean size of non-GCC clusters */
  double susceptibility;      /* Weighted mean (sum s^2 / sum s) - percolation
                                 indicator */
  double probabilty_congested_edges;
} PercolationState;

typedef struct {
  int id; /* Original cluster identifier from igraph */
  double value;
} Sorting;

#define SEED_DEBUG_ID 0
#define PHASE_ANALYSIS 1
#define PERCOLATION_MODE_RECENT 1
#define PERCOLATION_MODE_RATIO 2
#define PERCOLATION_MODE_FLOW 3
  
