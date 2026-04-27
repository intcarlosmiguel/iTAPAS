
#pragma once

#include "calc.h"
#include "define.h"
#include <igraph/igraph.h>
#include <igraph/igraph_vector.h>
#include <time.h>
const int MAX_ITER = 1000;
const int MAX_ITER_NEWTON = 100;

/* ============================================================================
   MERGE NODE HELPER FUNCTIONS
   ============================================================================
 */

/**
 * Inicializa os campos de merge nodes na estrutura BUSH.
 */
static inline void initMergeNodes(struct BUSH *bush, int numNodes) {
  bush->merges = NULL;
  bush->numMerges = 0;
  bush->nodeToMerge = (int *)malloc(numNodes * sizeof(int));
  for (int i = 0; i < numNodes; i++)
    bush->nodeToMerge[i] = -1; // -1 significa "não é merge node"
}

/**
 * Libera memória dos merge nodes.
 */
static inline void freeMergeNodes(struct BUSH *bush) {
  for (int i = 0; i < bush->numMerges; i++) {
    free(bush->merges[i].approaches);
    free(bush->merges[i].approachFlow);
  }
  free(bush->merges);
  free(bush->nodeToMerge);
  bush->merges = NULL;
  bush->numMerges = 0;
  bush->nodeToMerge = NULL;
}

/**
 * Reconstrói os merge nodes após modificar a topologia da bush.
 * Identifica nós com múltiplos arcos de entrada e cria MergeNodes para eles.
 */
static inline void createMergeNodes(struct BUSH *bush, igraph_t *Grafo,
                                    int fonte, int numNodes) {
  // Primeiro, limpa merge nodes existentes (mantém nodeToMerge alocado)
  for (int i = 0; i < bush->numMerges; i++) {
    free(bush->merges[i].approaches);
    free(bush->merges[i].approachFlow);
  }
  free(bush->merges);
  bush->numMerges = 0;
  for (int i = 0; i < numNodes; i++)
    bush->nodeToMerge[i] = -1;

  // Contar arcos de entrada para cada nó na bush
  int *inDegree = (int *)calloc(numNodes, sizeof(int));
  igraph_vector_int_t incident;
  igraph_vector_int_init(&incident, 0);

  for (int node = 0; node < numNodes; node++) {
    if (node == fonte)
      continue;
    igraph_incident(&bush->Grafo, &incident, node, IGRAPH_IN, IGRAPH_NO_LOOPS);
    inDegree[node] = igraph_vector_int_size(&incident);
    igraph_vector_int_clear(&incident);
  }

  // Contar quantos merge nodes precisamos
  int numMerges = 0;
  for (int i = 0; i < numNodes; i++)
    if (inDegree[i] > 1)
      numMerges++;

  if (numMerges == 0) {
    free(inDegree);
    igraph_vector_int_destroy(&incident);
    bush->merges = NULL;
    return;
  }

  // Alocar e preencher merge nodes
  bush->merges =
      (struct MergeNode *)malloc(numMerges * sizeof(struct MergeNode));
  bush->numMerges = numMerges;

  int mergeIdx = 0;
  for (int node = 0; node < numNodes; node++) {
    if (inDegree[node] <= 1)
      continue;

    struct MergeNode *m = &bush->merges[mergeIdx];
    m->node_id = node;
    m->numApproaches = inDegree[node];
    m->capacity = inDegree[node];
    m->approaches = (int *)malloc(m->numApproaches * sizeof(int));
    m->approachFlow = (double *)malloc(m->numApproaches * sizeof(double));
    m->SPlink = -1;
    m->LPlink = -1;
    m->divergenceNode = -1;

    // Preencher approaches com IDs dos arcos de entrada (do grafo original)
    igraph_incident(&bush->Grafo, &incident, node, IGRAPH_IN, IGRAPH_NO_LOOPS);
    for (int j = 0; j < m->numApproaches; j++) {
      int bushEdge = VECTOR(incident)[j];
      int originalEdge = (int)EAN(&bush->Grafo, "id", bushEdge);
      m->approaches[j] = originalEdge;
      m->approachFlow[j] = VECTOR(bush->flow_per_origin)[originalEdge];
    }
    igraph_vector_int_clear(&incident);

    bush->nodeToMerge[node] = mergeIdx;
    mergeIdx++;
  }

  free(inDegree);
  igraph_vector_int_destroy(&incident);
}

/**
 * Calcula e armazena o nó de divergência (LCA) para cada merge node.
 * Este é o ponto onde os caminhos SP e LP se separam.
 */
static inline void findDivergenceNodes(struct BUSH *bush, igraph_t *Grafo,
                                       int fonte, int numNodes) {
  bool *visited = (bool *)calloc(numNodes, sizeof(bool));

  for (int m = 0; m < bush->numMerges; m++) {
    struct MergeNode *merge = &bush->merges[m];

    // Se não há caminhos SP e LP válidos, skip
    if (merge->SPlink < 0 || merge->LPlink < 0 ||
        merge->SPlink == merge->LPlink) {
      merge->divergenceNode = -1;
      continue;
    }

    int nodeId = merge->node_id;

    // Reset visited
    for (int i = 0; i < numNodes; i++)
      visited[i] = false;
    visited[fonte] = true;

    // Trace SP path backwards, marking nodes
    int current = nodeId;
    while (current != fonte) {
      visited[current] = true;
      int edge = VECTOR(bush->paths.min_edges)[current];
      if (edge < 0)
        break;
      current = IGRAPH_FROM(Grafo, edge);
    }

    // Trace LP path backwards, looking for first visited node
    current = nodeId;
    merge->divergenceNode = -1;
    while (current != fonte) {
      int edge = VECTOR(bush->paths.max_edges)[current];
      if (edge < 0)
        break;
      current = IGRAPH_FROM(Grafo, edge);
      if (visited[current] && current != nodeId) {
        merge->divergenceNode = current;
        break;
      }
    }
  }

  free(visited);
}

static inline bool check_consistency(struct BUSH *bush, igraph_t *Grafo) {
  igraph_integer_t from, to, u, v;
  int id;
  for (long e = 0; e < igraph_ecount(&bush->Grafo); e++) {
    igraph_edge(&bush->Grafo, e, &u, &v);
    id = (int)EAN(&bush->Grafo, "id", e);
    igraph_edge(Grafo, id, &from, &to);
    if (from != u || to != v) {
      printf("Inconsistency found in bush edge %d: bush (%ld -> %ld), original "
             "graph (%ld -> %ld)\n",
             id, u, v, from, to);
      return false;
    }
  }
  return true;
}

static inline bool check_topological_order(struct BUSH *bush) {
  igraph_vector_int_t top_check;
  igraph_vector_int_init(&top_check, 0);
  if (igraph_topological_sorting(&bush->Grafo, &top_check, IGRAPH_OUT) !=
      IGRAPH_SUCCESS) {
    printf("Error: Graph has cycles after adding new edges.\n");
    exit(1);
  }
  int anterior = -1, node;
  for (int i = 0; i < igraph_vector_int_size(&top_check); i++) {
    node = VECTOR(top_check)[i];
    if (anterior != -1) {
      if (VECTOR(bush->topological_value)[node] <
          VECTOR(bush->topological_value)[anterior]) {
        printf("%d (%ld) %d (%ld)\n", node,
               VECTOR(bush->topological_value)[node], anterior,
               VECTOR(bush->topological_value)[anterior]);
        igraph_vector_int_destroy(&top_check);
        return false;
      }
    }
    anterior = node;
  }
  igraph_vector_int_destroy(&top_check);
  return true;
}

static inline void
att_topological_order(struct BUSH *bush,
                      igraph_t *Grafo, // O grafo original, para clareza
                      int fonte) {
  int N = igraph_vcount(Grafo);

  // --- Passo 1: Calcular o "in-degree" de cada nó ---
  // O in-degree é o número de arestas que chegam em um nó.
  // Esta é a informação crucial que nos diz quantas "dependências" um nó tem.
  igraph_vector_int_t in_degrees;
  igraph_vector_int_init(&in_degrees, N);
  igraph_degree(&bush->Grafo, &in_degrees, igraph_vss_all(), IGRAPH_IN,
                IGRAPH_NO_LOOPS);

  // --- Passo 2: Inicializar uma fila com os nós de in-degree 0 ---
  // Esses são os nós que não têm dependências. Dado que o grafo é conexo a
  // partir da 'fonte', a 'fonte' DEVE ser um desses nós.
  igraph_dqueue_int_t q;
  igraph_dqueue_int_init(&q, 0);
  // Na sua premissa, apenas 'fonte' terá in-degree 0 no subgrafo relevante.
  // Se houvesse outros, eles também deveriam ser adicionados aqui.
  igraph_dqueue_int_push(&q, fonte);

  // --- Preparação dos vetores de resultado ---
  // Usar 'clear' e 'resize' é mais seguro e explícito que 'null'.
  igraph_vector_int_clear(&bush->topological_order);
  igraph_vector_int_resize(&bush->topological_order, N);
  igraph_vector_int_clear(&bush->topological_value);
  igraph_vector_int_resize(&bush->topological_value, N);
  igraph_vector_int_fill(&bush->topological_value, 0); // Inicializa com zeros

  // A 'fonte' é o primeiro nó da ordem, e seu "nível" ou "valor" é 1.
  VECTOR(bush->topological_value)[fonte] = 1;

  int n = 0;                    // Contador para a posição na ordem topológica
  igraph_vector_int_t vizinhos; // Vetor para arestas incidentes
  igraph_vector_int_init(&vizinhos, 0);

  // --- Passo 3: Processar a fila ---
  // O loop continua enquanto houver nós sem dependências para processar.
  while (!igraph_dqueue_int_empty(&q)) {
    // Retira um nó da fila. Este nó já teve todas as suas dependências
    // resolvidas.
    int node = igraph_dqueue_int_pop(&q);
    // Adiciona o nó à nossa lista de ordem topológica.
    VECTOR(bush->topological_order)[n] = node;
    n++;

    // Agora, vamos "remover" este nó do grafo, informando seus vizinhos
    // que uma de suas dependências foi satisfeita.
    igraph_incident(&bush->Grafo, &vizinhos, node, IGRAPH_OUT, IGRAPH_NO_LOOPS);
    for (int j = 0; j < igraph_vector_int_size(&vizinhos); j++) {
      int edge_idx = VECTOR(vizinhos)[j];
      int neighbor = IGRAPH_TO(&bush->Grafo, edge_idx);
      // Mantive sua lógica original para 'topological_value'. Em um DAG,
      // isso calcula corretamente o comprimento do caminho mais longo da fonte
      // até cada nó, o que é uma operação muito útil.
      VECTOR(bush->topological_value)
      [neighbor] = fmax(VECTOR(bush->topological_value)[neighbor],
                        VECTOR(bush->topological_value)[node] + 1);

      // Decrementa o in-degree do vizinho.
      VECTOR(in_degrees)[neighbor]--;

      // Se o in-degree do vizinho chegou a 0, significa que todas as suas
      // dependências foram resolvidas. Agora ele está pronto para ser
      // processado.
      if (VECTOR(in_degrees)[neighbor] == 0) {
        igraph_dqueue_int_push(&q, neighbor);
      }
    }
    igraph_vector_int_clear(&vizinhos);
  }
  // --- Passo 4: Limpeza ---
  igraph_vector_int_destroy(&in_degrees);
  igraph_vector_int_destroy(&vizinhos);
  igraph_dqueue_int_destroy(&q);
}

static inline void removeUnusedArcs(struct BUSH *bush, igraph_t *Grafo,
                                    int fonte,
                                    struct PARAMETERS *BPR_PARAMETERS,
                                    igraph_vector_t *total_flow) {

  int i, j, edge_id;
  igraph_vector_int_t new_edges;
  igraph_vector_int_init(&new_edges, 0);
  int L1 = igraph_ecount(&bush->Grafo), L2 = 0;
  for (i = 0; i < L1; i++) {
    edge_id = EAN(&bush->Grafo, "id", i);
    if (VECTOR(bush->flow_per_origin)[edge_id] <= 1e-10) {
      bush->is_ingraph[edge_id] = false;
      igraph_vector_int_push_back(&new_edges, i);
      L2++;
    }
  }
  if (L2 == 0) {
    printf("No edges to remove from bush rooted at %d.\n", fonte);
    exit(0); // Se não houver arestas para remover, sai da função
  }
  igraph_es_t edges;
  igraph_es_vector(&edges, &new_edges);

  igraph_delete_edges(&bush->Grafo, edges);
  igraph_vector_int_clear(&new_edges);
  bool change_topological = false;
  igraph_integer_t from, to;
  int L = 0;
  L1 = igraph_ecount(&bush->Grafo);
  igraph_vector_t ids;
  igraph_vector_init(&ids, 0);
  for (i = 0; i < BPR_PARAMETERS->N; i++) {
    if (i == fonte)
      continue; // Pula a fonte
    edge_id = VECTOR(bush->paths.min_edges)[i];
    if (!bush->is_ingraph[edge_id]) {
      bush->is_ingraph[edge_id] = true;

      from = IGRAPH_FROM(Grafo, edge_id);
      to = IGRAPH_TO(Grafo, edge_id);

      igraph_vector_int_push_back(&new_edges, from);
      igraph_vector_int_push_back(&new_edges, to);
      igraph_vector_push_back(&ids, edge_id);

      if (VECTOR(bush->topological_value)[from] >=
          VECTOR(bush->topological_value)[to])
        change_topological = true;
      L++;
    }
  }
  igraph_add_edges(&bush->Grafo, &new_edges, NULL);
  for (int e = 0; e < L; e++) {
    igraph_integer_t u, v;
    igraph_cattribute_EAN_set(&bush->Grafo, "id", e + L1, VECTOR(ids)[e]);
  }
  att_topological_order(bush, Grafo,
                        fonte); // Atualiza a ordenação topológica da bush
  igraph_vector_destroy(&ids);
  igraph_vector_int_destroy(&new_edges);
  igraph_es_destroy(&edges);
}

static inline double findFlowDelta(struct BUSH *bush, igraph_t *Grafo, int lca,
                                   int alvo, struct PARAMETERS *BPR_PARAMETERS,
                                   igraph_vector_t *total_flow) {

  double mu = DBL_MAX;
  int current_node = alvo;
  while (current_node != lca) {
    int edge_id = VECTOR(bush->paths.max_edges)[current_node];
    if (edge_id == -1) {
      mu = 0;
      break;
    }
    mu = fmin(mu, VECTOR(bush->flow_per_origin)[edge_id]);
    current_node = IGRAPH_FROM(Grafo, edge_id);
  }
  if (mu == 0)
    return 0.0; // Se mu for zero, não há fluxo a ser ajustad
  double delta_x = 0.0, min_path_flow, max_path_flow, min_derivate_flow,
         max_derivate_flow, new_delta_x = 0.0;
  int edge_id;
  double denominator;
  for (int i = 0; i < MAX_ITER_NEWTON; i++) {

    min_path_flow = 0;
    min_derivate_flow = 0;
    max_path_flow = 0;
    max_derivate_flow = 0;
    current_node = alvo;
    while (current_node != lca) {
      edge_id = VECTOR(bush->paths.min_edges)[current_node];
      min_path_flow += single_BPR(VECTOR(*total_flow)[edge_id] + delta_x,
                                  VECTOR(BPR_PARAMETERS->cost_time)[edge_id],
                                  VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
      min_derivate_flow +=
          single_BPR_derivate(VECTOR(*total_flow)[edge_id] + delta_x,
                              VECTOR(BPR_PARAMETERS->cost_time)[edge_id],
                              VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
      current_node = IGRAPH_FROM(Grafo, edge_id);
    }
    current_node = alvo;
    while (current_node != lca) {

      edge_id = VECTOR(bush->paths.max_edges)[current_node];
      max_path_flow += single_BPR(VECTOR(*total_flow)[edge_id] - delta_x,
                                  VECTOR(BPR_PARAMETERS->cost_time)[edge_id],
                                  VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
      max_derivate_flow +=
          single_BPR_derivate(VECTOR(*total_flow)[edge_id] - delta_x,
                              VECTOR(BPR_PARAMETERS->cost_time)[edge_id],
                              VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
      current_node = IGRAPH_FROM(Grafo, edge_id);
    }

    denominator = max_derivate_flow + min_derivate_flow;
    if (denominator == 0)
      exit(0); // Evita divisão por zero
    new_delta_x = delta_x + (max_path_flow - min_path_flow) / denominator;
    if (fabs(new_delta_x - delta_x) < 1e-10) {
      delta_x = fmin(fmax(new_delta_x, 0.0), mu);
      break;
    }
    delta_x = new_delta_x;
  }
  return fmin(fmax(delta_x, 0.0),
              mu); // Garante que delta_x não seja negativo e não exceda mu
}

/**
 * Single Newton Step para um MergeNode.
 * Usa o divergenceNode cacheado em vez de recalcular.
 * Retorna true se houve shift de fluxo.
 */
static inline bool newtonShiftMerge(struct MergeNode *merge, struct BUSH *bush,
                                    igraph_t *Grafo,
                                    struct PARAMETERS *BPR_PARAMETERS,
                                    igraph_vector_t *total_flow,
                                    igraph_vector_t *time) {
  // Debug counters (static to persist across calls)
  static int call_count = 0;
  static int fail_splink = 0, fail_lca = 0, fail_mu = 0, fail_cost = 0,
             fail_denom = 0, fail_shift = 0;
  call_count++;

  // Verifica se há SPlink e LPlink válidos
  if (merge->SPlink < 0 || merge->LPlink < 0 ||
      merge->SPlink == merge->LPlink) {
    fail_splink++;
    return false;
  }

  int lca = merge->divergenceNode;
  if (lca < 0) {
    fail_lca++;
    return false;
  }

  int nodeId = merge->node_id;
  int spEdge = merge->approaches[merge->SPlink];
  int lpEdge = merge->approaches[merge->LPlink];

  // Calcular fluxo máximo que pode ser movido (min flow no LP path)
  double mu = DBL_MAX;
  int current = nodeId;
  while (current != lca) {
    int edge = VECTOR(bush->paths.max_edges)[current];
    if (edge < 0) {
      mu = 0;
      break;
    }
    mu = fmin(mu, VECTOR(bush->flow_per_origin)[edge]);
    current = IGRAPH_FROM(Grafo, edge);
  }

  if (mu <= MIN_LINK_FLOW) {
    fail_mu++;
    return false;
  }

  // Calcular custos dos caminhos SP e LP (do merge node até o LCA)
  double cost_SP = 0.0, cost_LP = 0.0;
  double der_SP = 0.0, der_LP = 0.0;
  int edge;

  // Custo do SP path
  current = nodeId;
  while (current != lca) {
    edge = VECTOR(bush->paths.min_edges)[current];
    if (edge < 0)
      break;
    cost_SP += VECTOR(*time)[edge];
    der_SP += single_BPR_derivate(VECTOR(*total_flow)[edge],
                                  VECTOR(BPR_PARAMETERS->cost_time)[edge],
                                  VECTOR(BPR_PARAMETERS->capacidade)[edge]);
    current = IGRAPH_FROM(Grafo, edge);
  }

  // Custo do LP path
  current = nodeId;
  while (current != lca) {
    edge = VECTOR(bush->paths.max_edges)[current];
    if (edge < 0)
      break;
    cost_LP += VECTOR(*time)[edge];
    der_LP += single_BPR_derivate(VECTOR(*total_flow)[edge],
                                  VECTOR(BPR_PARAMETERS->cost_time)[edge],
                                  VECTOR(BPR_PARAMETERS->capacidade)[edge]);
    current = IGRAPH_FROM(Grafo, edge);
  }

  // Se custos são praticamente iguais, já equilibrado
  if (fabs(cost_LP - cost_SP) < 1e-10)
    return false;

  // Single Newton step
  double denominator = der_SP + der_LP;
  if (denominator <= 1e-10)
    return false;

  double shift = (cost_LP - cost_SP) / denominator;

  // Clamp shift
  shift = fmax(0.0, fmin(shift, mu));

  if (shift <= MIN_LINK_FLOW)
    return false;

  // Aplicar shift: adicionar ao SP path, remover do LP path
  current = nodeId;
  while (current != lca) {
    edge = VECTOR(bush->paths.min_edges)[current];
    if (edge < 0)
      break;
    VECTOR(*total_flow)[edge] += shift;
    VECTOR(bush->flow_per_origin)[edge] += shift;
    VECTOR(*time)
    [edge] = single_BPR(VECTOR(*total_flow)[edge],
                        VECTOR(BPR_PARAMETERS->cost_time)[edge],
                        VECTOR(BPR_PARAMETERS->capacidade)[edge]);
    current = IGRAPH_FROM(Grafo, edge);
  }

  current = nodeId;
  while (current != lca) {
    edge = VECTOR(bush->paths.max_edges)[current];
    if (edge < 0)
      break;
    VECTOR(*total_flow)[edge] -= shift;
    if (VECTOR(*total_flow)[edge] < 1e-10)
      VECTOR(*total_flow)[edge] = 0;
    VECTOR(bush->flow_per_origin)[edge] -= shift;
    if (VECTOR(bush->flow_per_origin)[edge] < 1e-10)
      VECTOR(bush->flow_per_origin)[edge] = 0;
    VECTOR(*time)
    [edge] = single_BPR(VECTOR(*total_flow)[edge],
                        VECTOR(BPR_PARAMETERS->cost_time)[edge],
                        VECTOR(BPR_PARAMETERS->capacidade)[edge]);
    current = IGRAPH_FROM(Grafo, edge);
  }

  return true;
}

static inline void init_bush(struct BUSH *bush, igraph_t *Grafo, int id_bush,
                             struct PARAMETERS *BPR_PARAMETERS,
                             struct OD_MATRIX *OD, igraph_vector_t *flow) {
  igraph_vector_int_t inbound;

  igraph_vector_int_init(&inbound, 0);
  igraph_vector_int_init(&bush->topological_order, BPR_PARAMETERS->N);
  igraph_vector_int_init(&bush->topological_value, BPR_PARAMETERS->N);

  int fonte = OD->Elementos[id_bush].fonte, j, k, id;
  igraph_get_shortest_paths_dijkstra(Grafo, NULL, NULL, fonte, igraph_vss_all(),
                                     &BPR_PARAMETERS->cost_time, IGRAPH_OUT,
                                     NULL, &inbound);

  int alvo, edge_id, from, to;
  double tempo;
  bush->is_ingraph =
      (bool *)calloc(BPR_PARAMETERS->L,
                     sizeof(bool)); // Inicializa o vetor de arestas no grafo
  int target, antecessor;
  igraph_integer_t index;

  igraph_vector_int_t edges_vec;
  igraph_vector_int_init(&edges_vec, 0);
  igraph_vector_t ids;
  igraph_vector_init(&ids, 0);

  for (j = 0; j < BPR_PARAMETERS->N; j++) {
    if (j == fonte)
      continue; // Pula a fonte
    alvo = j;
    while (alvo != fonte) {
      edge_id = VECTOR(inbound)[alvo];
      alvo = IGRAPH_FROM(Grafo, edge_id);
      from = IGRAPH_FROM(Grafo, edge_id);
      to = IGRAPH_TO(Grafo, edge_id);
      if (bush->is_ingraph[edge_id])
        continue; // Se a aresta já está no grafo da bush, pula para a próxima
                  // iteração
      igraph_vector_int_push_back(&edges_vec, from);
      igraph_vector_int_push_back(&edges_vec, to);
      bush->is_ingraph[edge_id] =
          true; // Marca a aresta como parte do grafo da bush
      igraph_vector_push_back(&ids, edge_id);
    }
  }
  // igraph_vector_int_print(&edges_vec);
  igraph_empty(&bush->Grafo, BPR_PARAMETERS->N, IGRAPH_DIRECTED);
  igraph_es_t edges;
  igraph_es_vector(&edges, &edges_vec);
  igraph_add_edges(&bush->Grafo, &edges_vec, NULL);
  igraph_cattribute_EAN_setv(&bush->Grafo, "id", &ids);

  /* libera vetor auxiliar de ids */
  igraph_vector_destroy(&ids);

  /* atualiza ordenação topológica da bush */
  for (j = 0; j < bush->n_alvos; j++) {
    antecessor = VECTOR(OD->Elementos[id_bush].alvos)[j];
    while (antecessor != fonte) {
      index = VECTOR(inbound)[antecessor];
      if (index == -1)
        exit(0); // Se não encontrar a aresta, encerra o programa
      VECTOR(*flow)[index] += VECTOR(OD->Elementos[id_bush].volumes)[j];
      VECTOR(bush->flow_per_origin)
      [index] += VECTOR(OD->Elementos[id_bush].volumes)[j];
      antecessor = IGRAPH_FROM(Grafo, index);
    }
  }
  init_path(&bush->paths, BPR_PARAMETERS->N);
  // char filename[100];
  // sprintf(filename, "bush_%d.txt", fonte);
  // FILE* file = fopen(filename, "w");

  att_topological_order(bush, Grafo, fonte);

  // Inicializa merge nodes
  initMergeNodes(bush, BPR_PARAMETERS->N);

  // fclose(file);
  igraph_vector_int_destroy(&inbound);
  igraph_vector_int_destroy(&edges_vec);
  igraph_es_destroy(&edges);
}

static inline void max_Distance(igraph_vector_t *longest, struct BUSH *bush,
                                igraph_t *Grafo, igraph_integer_t from,
                                struct PARAMETERS *BPR_PARAMETERS,
                                igraph_vector_t *total_flow,
                                igraph_vector_t *time) {
  igraph_vector_init(longest, BPR_PARAMETERS->N);
  igraph_vector_fill(longest, -DBL_MAX); // "Infinito negativo"
  VECTOR(*longest)
  [from] = 0.0; // Inicializa a distância do nó de origem como 0.0

  igraph_vector_int_t incident_edges;
  igraph_vector_int_init(&incident_edges, 0);
  int i, index, id, edge_id, j;
  double fluxo, weight;
  igraph_integer_t eid, u, v_node;

  for (i = 0; i < igraph_vector_int_size(&bush->topological_order); ++i) {
    u = VECTOR(bush->topological_order)[i];
    igraph_incident(&bush->Grafo, &incident_edges, u, IGRAPH_OUT,
                    IGRAPH_NO_LOOPS);
    if (VECTOR(*longest)[u] == -DBL_MAX)
      continue;
    for (j = 0; j < igraph_vector_int_size(&incident_edges); ++j) {

      eid = VECTOR(incident_edges)[j];
      id = EAN(&bush->Grafo, "id", eid);
      v_node = IGRAPH_TO(Grafo, id);
      weight = VECTOR(*time)[id]; // Acessa o peso da aresta

      // Relaxamento para o maior caminho
      // Só estender de 'u' se 'u' já tem um caminho mais longo válido (não
      // -DBL_MAX)
      if (VECTOR(*longest)[u] + weight > VECTOR(*longest)[v_node]) {
        VECTOR(*longest)[v_node] = VECTOR(*longest)[u] + weight;
      }
    }
    igraph_vector_int_clear(&incident_edges); // Reutiliza o vetor de arestas
  }

  // Libera memória dos vetores locais
  igraph_vector_int_destroy(&incident_edges);
}

static inline void find_dag_shortest_longest_costs_and_parents(
    struct BUSH *bush, igraph_t *Grafo, igraph_integer_t from,
    struct PARAMETERS *BPR_PARAMETERS, igraph_vector_t *total_flow,
    igraph_vector_t *time) {
  int N = igraph_vcount(Grafo), j;
  igraph_vector_fill(&bush->paths.dist_shortest_local, DBL_MAX);
  igraph_vector_fill(&bush->paths.dist_longest_local,
                     -DBL_MAX); // "Infinito negativo"
  VECTOR(bush->paths.dist_shortest_local)[from] = 0.0;
  VECTOR(bush->paths.dist_longest_local)[from] = 0.0;
  igraph_vector_int_fill(&bush->paths.max_edges, -1);
  igraph_vector_int_fill(&bush->paths.min_edges, -1); // "Infinito negativo"

  igraph_vector_int_t incident_edges;
  igraph_vector_int_init(&incident_edges, 0);

  int i, index, id, edge_id;
  double fluxo, weight;
  igraph_integer_t eid, u, v_node;
  for (i = 0; i < igraph_vector_int_size(&bush->topological_order); ++i) {
    u = VECTOR(bush->topological_order)[i];
    // Se 'u' não é alcançável a partir de 'from' (para menor caminho),
    if (VECTOR(bush->paths.dist_shortest_local)[u] == DBL_MAX)
      continue;

    igraph_incident(&bush->Grafo, &incident_edges, u, IGRAPH_OUT,
                    IGRAPH_NO_LOOPS);

    for (j = 0; j < igraph_vector_int_size(&incident_edges); ++j) {

      eid = VECTOR(incident_edges)[j];
      v_node = IGRAPH_TO(&bush->Grafo, eid);
      id = (int)EAN(&bush->Grafo, "id", eid);
      fluxo = VECTOR(bush->flow_per_origin)[id]; // Acessa o fluxo da aresta
      weight = VECTOR(*time)[id];                // Acessa o peso da aresta
      if (VECTOR(bush->paths.dist_shortest_local)[u] + weight <
          VECTOR(bush->paths.dist_shortest_local)[v_node]) {

        VECTOR(bush->paths.dist_shortest_local)
        [v_node] = VECTOR(bush->paths.dist_shortest_local)[u] + weight;
        VECTOR(bush->paths.min_edges)[v_node] = id;
      }

      // Relaxamento para o maior caminho
      // Só estender de 'u' se 'u' já tem um caminho mais longo válido (não
      // -DBL_MAX)
      if ((VECTOR(bush->paths.dist_longest_local)[u] != -DBL_MAX) &&
          (fluxo > 0)) {
        if (VECTOR(bush->paths.dist_longest_local)[u] + weight >
            VECTOR(bush->paths.dist_longest_local)[v_node]) {

          VECTOR(bush->paths.dist_longest_local)
          [v_node] = VECTOR(bush->paths.dist_longest_local)[u] + weight;
          VECTOR(bush->paths.max_edges)[v_node] = id;
        }
      }
    }
    igraph_vector_int_clear(&incident_edges); // Reutiliza o vetor de arestas
  }
  // Libera memória dos vetores locais
  igraph_vector_int_destroy(&incident_edges);

  // ============================================================================
  // ATUALIZAR SPlink e LPlink NOS MERGE NODES
  // ============================================================================
  for (int m = 0; m < bush->numMerges; m++) {
    struct MergeNode *merge = &bush->merges[m];
    int nodeId = merge->node_id;

    // Atualizar fluxos nos approaches
    for (int a = 0; a < merge->numApproaches; a++) {
      int edgeId = merge->approaches[a];
      merge->approachFlow[a] = VECTOR(bush->flow_per_origin)[edgeId];
    }

    // Encontrar SPlink (arco com menor custo)
    merge->SPlink = -1;
    double minCost = DBL_MAX;
    for (int a = 0; a < merge->numApproaches; a++) {
      int edgeId = merge->approaches[a];
      int tailNode = IGRAPH_FROM(Grafo, edgeId);
      double cost = VECTOR(bush->paths.dist_shortest_local)[tailNode] +
                    VECTOR(*time)[edgeId];
      if (cost < minCost) {
        minCost = cost;
        merge->SPlink = a;
      }
    }

    // Encontrar LPlink (arco com maior custo E fluxo > 0 - LONGEST_USED_PATH
    // optimization)
    merge->LPlink = -1;
    double maxCost = -DBL_MAX;
    for (int a = 0; a < merge->numApproaches; a++) {
      int edgeId = merge->approaches[a];
      // Só considerar arcos com fluxo (LONGEST_USED_PATH optimization)
      if (merge->approachFlow[a] <= MIN_LINK_FLOW)
        continue;

      int tailNode = IGRAPH_FROM(Grafo, edgeId);
      double cost = VECTOR(bush->paths.dist_longest_local)[tailNode] +
                    VECTOR(*time)[edgeId];
      if (cost > maxCost) {
        maxCost = cost;
        merge->LPlink = a;
      }
    }
  }
}

static inline int find_lca(int fonte, int alvo, struct BUSH *bush,
                           igraph_t *Grafo, bool *visited, int visited_size) {

  int edge_id, current_node = alvo;
  int next;
  int edge_min_parrent = VECTOR(bush->paths.min_edges)[alvo];
  int edge_max_parrent = VECTOR(bush->paths.max_edges)[alvo];
  if (edge_min_parrent == -1 || edge_max_parrent == -1)
    return -1;

  int min_parrent = IGRAPH_FROM(Grafo, edge_min_parrent);
  int max_parrent = IGRAPH_FROM(Grafo, edge_max_parrent);

  if (min_parrent == max_parrent)
    return -1;

  // Reset visited buffer (caller provides pre-allocated buffer)
  memset(visited, 0, visited_size * sizeof(bool));
  visited[alvo] = true;
  visited[fonte] = true;
  while (current_node != fonte) {

    edge_id = VECTOR(bush->paths.min_edges)[current_node];
    current_node = IGRAPH_FROM(Grafo, edge_id);
    visited[current_node] = true;
    if (edge_id == -1)
      break;
  }
  current_node = max_parrent;
  bool find = false;
  int next_node;
  int a = 0;
  while (current_node != fonte) {
    edge_id = VECTOR(bush->paths.max_edges)[current_node];
    next_node = IGRAPH_FROM(Grafo, edge_id);
    if (edge_id == -1) {
      return -1;
    }
    if ((visited[current_node]) && (current_node != alvo)) {
      if (current_node != fonte)
        return current_node;
      return -1;
    }
    if (current_node == next_node)
      a++;
    if (a > 10) {
      printf("Loop detected in max path at node %d\n", current_node);
      exit(0);
    }
    current_node = next_node;
  }
  if ((current_node == fonte) && (visited[current_node])) {
    return current_node;
  }
  return -1;
}

static inline void shift_flow(double delta, int alvo, int lca,
                              struct BUSH *bush, const igraph_t *Grafo,
                              igraph_vector_t *total_flow,
                              igraph_vector_t *time,
                              struct PARAMETERS *BPR_PARAMETERS) {
  int current_node;
  int edge_id;

  // Adiciona fluxo ao caminho mínimo
  current_node = alvo;
  while (current_node != lca) {
    edge_id = VECTOR(bush->paths.min_edges)[current_node];
    if (edge_id == -1)
      break;
    VECTOR(*total_flow)[edge_id] += delta;
    VECTOR(bush->flow_per_origin)[edge_id] += delta; // Atualiza o fluxo da bush
    VECTOR(*time)
    [edge_id] = single_BPR(VECTOR(*total_flow)[edge_id],
                           VECTOR(BPR_PARAMETERS->cost_time)[edge_id],
                           VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
    current_node = IGRAPH_FROM(Grafo, edge_id);
  }

  // Remove fluxo do caminho máximo
  current_node = alvo;
  while (current_node != lca) {
    edge_id = VECTOR(bush->paths.max_edges)[current_node];
    if (edge_id == -1)
      break;
    VECTOR(*total_flow)[edge_id] -= delta;
    if (VECTOR(*total_flow)[edge_id] < 1e-10)
      VECTOR(*total_flow)[edge_id] = 0;              // Evita fluxo negativo
    VECTOR(bush->flow_per_origin)[edge_id] -= delta; // Atualiza o fluxo da bush
    if (VECTOR(bush->flow_per_origin)[edge_id] < 1e-10)
      VECTOR(bush->flow_per_origin)[edge_id] = 0; // Evita fluxo negativo
    VECTOR(*time)
    [edge_id] = single_BPR(VECTOR(*total_flow)[edge_id],
                           VECTOR(BPR_PARAMETERS->cost_time)[edge_id],
                           VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
    current_node = IGRAPH_FROM(Grafo, edge_id);
  }
}

/**
 * @brief Verifica se a bush é ótima e, se não, a melhora adicionando arcos.
 * Implementa a Seção 2.5: "Improved bush".
 */
static inline bool improve_bush(struct BUSH *bush, igraph_t *Grafo,
                                struct PARAMETERS *BPR_PARAMETERS, int fonte,
                                igraph_vector_t *total_flow,
                                igraph_vector_t *time) {

  int i, from, to;
  igraph_integer_t index;
  igraph_vector_int_t new_edges;
  igraph_vector_t ids;
  igraph_vector_init(&ids, 0);
  igraph_vector_int_init(&new_edges, 0); // Inicializa o vetor de novas arestas
  igraph_vector_t longest;
  clock_t t1, t2;

  max_Distance(&longest, bush, Grafo, fonte, BPR_PARAMETERS, total_flow, time);

  int L = 0;

  int L1 = igraph_ecount(&bush->Grafo);
  bool att_top = false;
  for (i = 0; i < BPR_PARAMETERS->L; i++) {

    from = IGRAPH_FROM(Grafo, i);
    to = IGRAPH_TO(Grafo, i);

    if (!bush->is_ingraph[i]) {
      if (VECTOR(longest)[from] == -DBL_MAX || VECTOR(longest)[to] == -DBL_MAX)
        continue; // Se não há caminho para a fonte ou destino, pula esta aresta
      if (VECTOR(longest)[to] - VECTOR(longest)[from] > 1e-10) {
        // Adiciona a aresta ao grafo da bush
        bush->is_ingraph[i] =
            true; // Marca a aresta como parte do grafo da bush
        VECTOR(bush->flow_per_origin)
        [i] = 0; // Inicializa o fluxo para esta aresta
        igraph_vector_int_push_back(
            &new_edges, from); // Adiciona a aresta ao vetor de novas arestas
        igraph_vector_int_push_back(&new_edges, to);
        L++;
        // igraph_vector_int_push_back(&bush->edge_id, i);
        igraph_vector_push_back(&ids, i);
        if (VECTOR(bush->topological_value)[from] >=
            VECTOR(bush->topological_value)[to])
          att_top = true;
      }
    }
  }

  if (L == 0) {
    igraph_vector_destroy(&longest);
    igraph_vector_int_destroy(&new_edges);
    igraph_vector_destroy(&ids);
    return false; // Se não há novas arestas, retorna falso
  }

  igraph_es_t edges;
  igraph_es_vector(&edges, &new_edges);
  t1 = clock();
  igraph_add_edges(&bush->Grafo, &new_edges, NULL);
  for (int e = 0; e < L; e++) {
    igraph_integer_t u, v;
    igraph_cattribute_EAN_set(&bush->Grafo, "id", e + L1, VECTOR(ids)[e]);
  }

  igraph_es_destroy(&edges);
  igraph_vector_int_destroy(
      &new_edges); // Libera memória do vetor de novas arestas
  igraph_vector_destroy(
      &longest); // Libera memória do vetor de maiores distâncias
  igraph_vector_destroy(&ids);

  if (att_top)
    att_topological_order(bush, Grafo,
                          fonte); // Atualiza a ordenação topológica da bush
  return true;
}

static inline void ApplyWarmStart(struct BUSH **bushes, struct OD_MATRIX *OD,
                                  igraph_t *Grafo,
                                  struct PARAMETERS *BPR_PARAMETERS,
                                  igraph_vector_t *solucao,
                                  bool PROPORTIONAL_WARM_START) {
  int i, j, target;
  for (i = 0; i < OD->size; i++) {
    if (!PROPORTIONAL_WARM_START) {
      // STANDARD WARM START (All-or-Nothing to Shortest Path of the Bush)
      for (j = 0; j < (*bushes)[i].n_alvos; j++) {
        target = VECTOR(OD->Elementos[i].alvos)[j];
        int edge_id;
        while (target != OD->Elementos[i].fonte) {
          edge_id = VECTOR((*bushes)[i].paths.min_edges)[target];
          VECTOR((*bushes)[i].flow_per_origin)
          [edge_id] += VECTOR(OD->Elementos[i].warm_volumes)[j];
          target = IGRAPH_FROM(Grafo, edge_id);
        }
      }
    } else {
      /* PROPORTIONAL WARM START: Distribute additional demand using a forward
       * pass with splitting ratios based on existing flows. We distribute the
       * new demand 'd' across all currently used incoming links to the node,
       * recursively (from target back to source). At node v, pull d_k from
       * neighbor u_k proportional to flow on edge (u_k, v). */
      double *demand_at_node =
          (double *)calloc(BPR_PARAMETERS->N, sizeof(double));

      // 1. Accumulate demand at all targets (destinations) for this origin
      for (j = 0; j < (*bushes)[i].n_alvos; j++) {
        target = VECTOR(OD->Elementos[i].alvos)[j];
        demand_at_node[target] += VECTOR(OD->Elementos[i].warm_volumes)[j];
      }

      // 2. Iterate nodes in REVERSE topological order (generally Target ->
      // Source) Standard topological order lists parents before children. So
      // reverse lists children before parents.
      int current_idx =
          igraph_vector_int_size(&(*bushes)[i].topological_order) - 1;

      for (int k = current_idx; k >= 0; k--) {
        int u = VECTOR((*bushes)[i].topological_order)[k];
        double demand = demand_at_node[u];

        if (demand <= 1e-15)
          continue;
        if (u == OD->Elementos[i].fonte)
          continue;

        // Distribute 'demand' to incoming edges (v -> u)
        igraph_vector_int_t incident;
        igraph_vector_int_init(&incident, 0);
        igraph_incident(&(*bushes)[i].Grafo, &incident, u, IGRAPH_IN,
                        IGRAPH_NO_LOOPS);

        int n_incident = igraph_vector_int_size(&incident);
        double total_in_flow = 0.0;

        // Calculate total incoming flow
        for (int inc = 0; inc < n_incident; inc++) {
          int edge_id_bush = VECTOR(incident)[inc];
          int original_edge_id =
              (int)EAN(&(*bushes)[i].Grafo, "id", edge_id_bush);
          total_in_flow +=
              VECTOR((*bushes)[i].flow_per_origin)[original_edge_id];
        }

        if (total_in_flow > 1e-10) {
          // Distribute proportionally
          for (int inc = 0; inc < n_incident; inc++) {
            int edge_id_bush = VECTOR(incident)[inc];
            int original_edge_id =
                (int)EAN(&(*bushes)[i].Grafo, "id", edge_id_bush);
            double flow =
                VECTOR((*bushes)[i].flow_per_origin)[original_edge_id];
            double fraction = flow / total_in_flow;
            double add = demand * fraction;

            VECTOR((*bushes)[i].flow_per_origin)
            [original_edge_id] += add; // Update flow on edge
            int v = IGRAPH_FROM(&(*bushes)[i].Grafo, edge_id_bush);
            demand_at_node[v] += add; // Propagate demand back to parent node v
          }
        } else {
          // Fallback: If no existing flow, push to Shortest Path (min_edge)
          int min_e = VECTOR((*bushes)[i].paths.min_edges)[u];
          if (min_e != -1) {
            VECTOR((*bushes)[i].flow_per_origin)[min_e] += demand;
            int v = IGRAPH_FROM(Grafo, min_e);
            demand_at_node[v] += demand;
          }
        }
        igraph_vector_int_destroy(&incident);
      }
      free(demand_at_node);
    }
  }

  // Update total solution flow
  for (int i = 0; i < BPR_PARAMETERS->L; i++)
    for (j = 0; j < OD->size; j++)
      VECTOR(*solucao)[i] += VECTOR((*bushes)[j].flow_per_origin)[i];
}

static inline void Dial(igraph_t *Grafo, struct OD_MATRIX *OD,
                        struct PARAMETERS *BPR_PARAMETERS,
                        igraph_vector_t *solucao, struct BUSH **bushes,
                        bool WARM_START, bool PROPORTIONAL_WARM_START) {

  igraph_vector_init(solucao, BPR_PARAMETERS->L);
  int i, j, k, iter = 1, fonte, target;
  if (WARM_START && bushes != NULL) {
    ApplyWarmStart(bushes, OD, Grafo, BPR_PARAMETERS, solucao,
                   PROPORTIONAL_WARM_START);
  } else {
    *bushes = (struct BUSH *)malloc(OD->size * sizeof(struct BUSH));
    igraph_bool_t has_multiple;
    for (i = 0; i < OD->size; i++) {
      (*bushes)[i].n_alvos = igraph_vector_int_size(&OD->Elementos[i].alvos);
      igraph_vector_init(&(*bushes)[i].flow_per_origin, BPR_PARAMETERS->L);

      init_bush(&(*bushes)[i], Grafo, i, BPR_PARAMETERS, OD, solucao);
      igraph_has_multiple(&(*bushes)[i].Grafo, &has_multiple);
      if (has_multiple) {
        printf("Bush %d has multiple edges!\n", i);
        exit(0); // Encerra o programa para depuração
      }
    }
  }

  igraph_vector_t time;
  igraph_vector_init(&time, BPR_PARAMETERS->L);
  BPR(&time, BPR_PARAMETERS, solucao);
  
  // Print initial solution values
  // FILE *file = fopen("gap.txt", "w");
  double GAP, previous_GAP = DBL_MAX, a, b, c;
  if (WARM_START) {
    // Use APPROXIMATE gap (avoid Dijkstra) for initial check
    //GAP = relative_gap_approximate(solucao, *bushes, BPR_PARAMETERS, OD);

    if (GAP < EPSILON) {

      // printf("Warm start efficient (Approx GAP=%e). Skipping Dial.\n", GAP);
      igraph_vector_destroy(&time);
      return;
    }
    // Lazy mode: if GAP is small (but > EPSILON), skip improve_bush often
    bool lazy_mode = (GAP < 1e-1);
  } else {
    // Cold start: GAP is definitely high, so don't waste time calculating it
    GAP = DBL_MAX;
  }

  clock_t start_time, end_time, t1, t2;

  int valor = count_files_in_dir("./output/dial");
  char filename[10000];
  // sprintf(filename, "./output/dial/gap_progression_%d.txt", valor + 1);
  // FILE* gap_file = fopen(filename, "w");
  bool check;
  for (int iter = 0; iter < MAX_ITER; iter++) {
    start_time = clock();

    a = b = 0;

    for (i = 0; i < OD->size; i++) {
      fonte = OD->Elementos[i].fonte;


      bool topo_changed = improve_bush(&(*bushes)[i], Grafo, BPR_PARAMETERS,
                                       fonte, solucao, &time);
      // Recria merge nodes se a topologia da bush mudou
      if (topo_changed)
        createMergeNodes(&(*bushes)[i], Grafo, fonte, BPR_PARAMETERS->N);
      find_dag_shortest_longest_costs_and_parents(
          &(*bushes)[i], Grafo, fonte, BPR_PARAMETERS, solucao, &time);
      // Calcula divergence node (LCA) para cada merge node após scan SP/LP
      findDivergenceNodes(&(*bushes)[i], Grafo, fonte, BPR_PARAMETERS->N);
    }


    // =========================================================================
    // FLOW EQUILIBRATION - MergeNode topological pass (Algorithm B style)
    // Itera em ordem topológica DESCENDENTE; só processa merge nodes.
    // O LCA (divergenceNode) já está cacheado em cada MergeNode —
    // elimina o custo de find_lca a cada iteração.
    // =========================================================================
    int shift_count = 0;
    for (i = 0; i < OD->size; i++) {
      fonte = OD->Elementos[i].fonte;
      int N_topo = igraph_vector_int_size(&(*bushes)[i].topological_order);

      bool flow_was_shifted = false;
      for (k = N_topo - 1; k >= 0; k--) {
        int node = VECTOR((*bushes)[i].topological_order)[k];
        if (node == fonte)
          continue;

        int m = (*bushes)[i].nodeToMerge[node]; // O(1)
        if (m < 0)
          continue; // não é merge node — nenhum balanceamento necessário

        bool shifted = newtonShiftMerge(&(*bushes)[i].merges[m], &(*bushes)[i],
                                        Grafo, BPR_PARAMETERS, solucao, &time);
        if (shifted) {
          flow_was_shifted = true;
          shift_count++;
        }
      }
      if (flow_was_shifted)
        removeUnusedArcs(&(*bushes)[i], Grafo, fonte, BPR_PARAMETERS, solucao);
    }

    GAP = relative_gap(solucao, Grafo, BPR_PARAMETERS, OD, NULL, NULL);
    end_time = clock();

    if (GAP < EPSILON)
      break;
    previous_GAP = GAP;
  }
  printf("Iterações: %d, GAP: %e\n", iter, GAP);

  igraph_vector_destroy(&time);
  // fclose(gap_file);  
}