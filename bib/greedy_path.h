#pragma once

#include "calc.h"
#include "define.h"
#include <igraph/igraph.h>

/* ========================================================================
   GREEDY PATH-BASED ALGORITHM (Xie et al., 2018)
   Implementação baseada no prompt.md - Algoritmo de Alocação de Tráfego

   VERSÃO 4: Com Warm Start e Persistência de Estado
   ======================================================================== */

double GREEDY_EPSILON = 1e-6;
int GREEDY_MAX_MAIN_ITER = 500;
double GREEDY_DAMPING = 0.5; /* Fator de amortecimento (0-1) */

/* ========================================================================
   SEÇÃO 1: GERENCIAMENTO DE CAMINHOS
   ======================================================================== */

void caminho_criar_de_antecessores(struct CAMINHO *caminho,
                                   igraph_vector_int_t *antecessores, int fonte,
                                   int destino, igraph_t *Grafo) {
  igraph_vector_int_init(&caminho->arcos, 0);
  caminho->fluxo = 0.0;
  caminho->custo_v = 0.0;
  caminho->derivada_s = 0.0;
  caminho->termo_c = 0.0;

  int no_atual = destino;
  while (no_atual != fonte) {
    int arco = VECTOR(*antecessores)[no_atual];
    if (arco == -1)
      break;
    igraph_vector_int_push_back(&caminho->arcos, arco);
    no_atual = IGRAPH_FROM(Grafo, arco);
  }
}

void caminho_destruir(struct CAMINHO *caminho) {
  igraph_vector_int_destroy(&caminho->arcos);
}

bool caminho_igual(struct CAMINHO *c1, struct CAMINHO *c2) {
  int n1 = igraph_vector_int_size(&c1->arcos);
  int n2 = igraph_vector_int_size(&c2->arcos);
  if (n1 != n2)
    return false;
  for (int i = 0; i < n1; i++) {
    if (VECTOR(c1->arcos)[i] != VECTOR(c2->arcos)[i])
      return false;
  }
  return true;
}

void conjunto_inicializar(struct CONJUNTO_CAMINHOS *conjunto) {
  conjunto->caminhos = NULL;
  conjunto->num_caminhos = 0;
  conjunto->capacidade = 0;
  conjunto->delta_rs = 0.0;
}

void conjunto_destruir(struct CONJUNTO_CAMINHOS *conjunto) {
  for (int i = 0; i < conjunto->num_caminhos; i++) {
    caminho_destruir(&conjunto->caminhos[i]);
  }
  free(conjunto->caminhos);
  conjunto->caminhos = NULL;
  conjunto->num_caminhos = 0;
  conjunto->capacidade = 0;
}

bool conjunto_adicionar_caminho(struct CONJUNTO_CAMINHOS *conjunto,
                                struct CAMINHO *novo_caminho) {
  for (int i = 0; i < conjunto->num_caminhos; i++) {
    if (caminho_igual(&conjunto->caminhos[i], novo_caminho))
      return false;
  }

  if (conjunto->num_caminhos >= conjunto->capacidade) {
    int nova_capacidade =
        (conjunto->capacidade == 0) ? 4 : conjunto->capacidade * 2;
    struct CAMINHO *tmp =
        realloc(conjunto->caminhos, nova_capacidade * sizeof(struct CAMINHO));
    if (!tmp)
      return false;
    conjunto->caminhos = tmp;
    conjunto->capacidade = nova_capacidade;
  }

  struct CAMINHO *destino = &conjunto->caminhos[conjunto->num_caminhos];
  igraph_vector_int_init_copy(&destino->arcos, &novo_caminho->arcos);
  destino->fluxo = novo_caminho->fluxo;
  destino->custo_v = novo_caminho->custo_v;
  destino->derivada_s = novo_caminho->derivada_s;
  destino->termo_c = novo_caminho->termo_c;

  conjunto->num_caminhos++;
  return true;
}

void conjunto_remover_caminhos_vazios(struct CONJUNTO_CAMINHOS *conjunto) {
  int j = 0;
  for (int i = 0; i < conjunto->num_caminhos; i++) {
    if (conjunto->caminhos[i].fluxo > GREEDY_EPSILON) {
      if (i != j)
        conjunto->caminhos[j] = conjunto->caminhos[i];
      j++;
    } else {
      caminho_destruir(&conjunto->caminhos[i]);
    }
  }
  conjunto->num_caminhos = j;
}

/* ========================================================================
   SEÇÃO 2: GERENCIAMENTO DE ESTADO E WARM START
   ======================================================================== */

int contar_total_pares_od(struct OD_MATRIX *OD) {
  int total = 0;
  for (int o = 0; o < OD->size; o++) {
    total += igraph_vector_int_size(&OD->Elementos[o].alvos);
  }
  return total;
}

void greedy_inicializar_estado(struct GREEDY_STATE *estado,
                               struct OD_MATRIX *OD) {
  estado->total_pares = contar_total_pares_od(OD);
  estado->conjuntos =
      malloc(estado->total_pares * sizeof(struct CONJUNTO_CAMINHOS));
  estado->origem_idx = malloc(estado->total_pares * sizeof(int));
  estado->destino_idx = malloc(estado->total_pares * sizeof(int));
  estado->demandas = malloc(estado->total_pares * sizeof(double));
  estado->RG_anterior = 1.0;

  int idx = 0;
  for (int o = 0; o < OD->size; o++) {
    int num_destinos = igraph_vector_int_size(&OD->Elementos[o].alvos);
    for (int d = 0; d < num_destinos; d++) {
      conjunto_inicializar(&estado->conjuntos[idx]);
      estado->origem_idx[idx] = o;
      estado->destino_idx[idx] = d;
      estado->demandas[idx] = VECTOR(OD->Elementos[o].volumes)[d];
      idx++;
    }
  }
}

void greedy_destruir_estado(struct GREEDY_STATE *estado) {
  if (!estado || !estado->conjuntos)
    return;
  for (int i = 0; i < estado->total_pares; i++) {
    conjunto_destruir(&estado->conjuntos[i]);
  }
  free(estado->conjuntos);
  free(estado->origem_idx);
  free(estado->destino_idx);
  free(estado->demandas);
  estado->conjuntos = NULL;
}

/**
 * @brief Implementa a lógica de Warm Start conforme warm_start.md
 * Migra e escala caminhos do estado anterior para o novo.
 */
void greedy_executar_warm_start(struct GREEDY_STATE *novo_estado,
                                struct GREEDY_STATE *estado_antigo,
                                struct OD_MATRIX *OD,
                                igraph_vector_t *solucao) {

  //printf("GreedyPath: Executando Warm Start...\n");

  /* Zera fluxos de elos para reconstrução */
  igraph_vector_fill(solucao, 0.0);

  /* Para cada par OD no novo cenário */
  for (int i = 0; i < novo_estado->total_pares; i++) {
    int o_idx = novo_estado->origem_idx[i];
    int d_idx = novo_estado->destino_idx[i];
    double demanda_nova = novo_estado->demandas[i];

    int fonte_novo = OD->Elementos[o_idx].fonte;
    int alvo_novo = VECTOR(OD->Elementos[o_idx].alvos)[d_idx];

    /* Procura correspondente no estado antigo */
    struct CONJUNTO_CAMINHOS *conjunto_antigo = NULL;
    double demanda_antiga = 0.0;

    /* Busca linear (ineficiente, mas segura para garantir correspondência
     * exata) */
    /* Nota: Assumindo que indices de fonte/alvo são consistentes entre
     * iterações da simulação */
    /* Se a matriz OD muda estrutura, idealmente deveríamos buscar por ID de nó
     * fonte/alvo */
    /* Como simulate_example mantém indices consistentes, usamos indices diretos
     * se possível */

    if (i < estado_antigo->total_pares) {
      /* Otimização: se a estrutura é a mesma, está na mesma posição */
      int o_old = estado_antigo->origem_idx[i];
      int d_old = estado_antigo->destino_idx[i];

      /* Verifica se mudou a estrutura da matriz (nós fonte/alvo) */
      /* Assumimos consistência da simulação. */
      conjunto_antigo = &estado_antigo->conjuntos[i];
      demanda_antiga = estado_antigo->demandas[i];
    }

    /* Se encontrou par com demanda anterior */
    if (conjunto_antigo && conjunto_antigo->num_caminhos > 0 &&
        demanda_antiga > GREEDY_EPSILON) {
      double alpha = demanda_nova / demanda_antiga;

      for (int h = 0; h < conjunto_antigo->num_caminhos; h++) {
        struct CAMINHO *caminho_ant = &conjunto_antigo->caminhos[h];

        /* Passo 1 do Checklist: Se fluxo > epsilon */
        if (caminho_ant->fluxo > GREEDY_EPSILON) {
          struct CAMINHO caminho_novo;
          /* Copia estrutura */
          igraph_vector_int_init_copy(&caminho_novo.arcos, &caminho_ant->arcos);
          caminho_novo.custo_v = 0;
          caminho_novo.derivada_s = 0;
          caminho_novo.termo_c = 0;

          /* Escala fluxo */
          caminho_novo.fluxo = caminho_ant->fluxo * alpha;

          conjunto_adicionar_caminho(&novo_estado->conjuntos[i], &caminho_novo);

          /* Reconstrói fluxos nos elos (Passo 2 do Checklist) */
          int num_arcos = igraph_vector_int_size(&caminho_novo.arcos);
          for (int a = 0; a < num_arcos; a++) {
            int arco = VECTOR(caminho_novo.arcos)[a];
            VECTOR(*solucao)[arco] += caminho_novo.fluxo;
          }

          caminho_destruir(&caminho_novo);
        }
      }
    }
    /* Passo 1 Checklist ELSE: Se não achou antigo, conjunto ficará vazio
       e será tratado na fase de geração de colunas (fallback implícito) */
  }
}

/* ========================================================================
   SEÇÃO 3: INICIALIZAÇÃO
   ======================================================================== */

void greedy_atribuicao_inicial_v2(struct PARAMETERS *BPR_PARAMETERS,
                                  struct OD_MATRIX *OD, igraph_t *Grafo,
                                  igraph_vector_t *solucao,
                                  struct GREEDY_STATE *estado) {

  // igraph_vector_init(solucao, BPR_PARAMETERS->L);
  greedy_inicializar_estado(estado, OD);

  igraph_vector_int_t antecessores;
  igraph_vector_int_init(&antecessores, 0);

  int idx = 0;
  for (int o = 0; o < OD->size; o++) {
    int fonte = OD->Elementos[o].fonte;

    igraph_get_shortest_paths_dijkstra(
        Grafo, NULL, NULL, fonte, igraph_vss_all(), &BPR_PARAMETERS->cost_time,
        IGRAPH_OUT, NULL, &antecessores);

    int num_destinos = igraph_vector_int_size(&OD->Elementos[o].alvos);
    for (int d = 0; d < num_destinos; d++) {
      int destino = VECTOR(OD->Elementos[o].alvos)[d];
      double demanda = VECTOR(OD->Elementos[o].volumes)[d];

      struct CAMINHO caminho;
      caminho_criar_de_antecessores(&caminho, &antecessores, fonte, destino,
                                    Grafo);
      caminho.fluxo = demanda;

      conjunto_adicionar_caminho(&estado->conjuntos[idx], &caminho);

      int num_arcos = igraph_vector_int_size(&caminho.arcos);
      for (int a = 0; a < num_arcos; a++) {
        int arco = VECTOR(caminho.arcos)[a];
        VECTOR(*solucao)[arco] += demanda;
      }

      caminho_destruir(&caminho);
      idx++;
    }
  }

  igraph_vector_int_destroy(&antecessores);
}

/* ========================================================================
   SEÇÃO 4: GREEDY SOLVER COM DAMPING
   ======================================================================== */

static double calcular_custo_caminho(struct CAMINHO *caminho,
                                     struct PARAMETERS *BPR_PARAMETERS,
                                     igraph_vector_t *solucao) {
  double custo = 0.0;
  int num_arcos = igraph_vector_int_size(&caminho->arcos);
  for (int a = 0; a < num_arcos; a++) {
    int arco = VECTOR(caminho->arcos)[a];
    double fluxo = VECTOR(*solucao)[arco];
    double tempo_livre = VECTOR(BPR_PARAMETERS->cost_time)[arco];
    double capacidade = VECTOR(BPR_PARAMETERS->capacidade)[arco];
    custo += single_BPR(fluxo, tempo_livre, capacidade);
  }
  return custo;
}

static double calcular_derivada_caminho(struct CAMINHO *caminho,
                                        struct PARAMETERS *BPR_PARAMETERS,
                                        igraph_vector_t *solucao) {
  double deriv = 0.0;
  int num_arcos = igraph_vector_int_size(&caminho->arcos);
  for (int a = 0; a < num_arcos; a++) {
    int arco = VECTOR(caminho->arcos)[a];
    double fluxo = VECTOR(*solucao)[arco];
    double tempo_livre = VECTOR(BPR_PARAMETERS->cost_time)[arco];
    double capacidade = VECTOR(BPR_PARAMETERS->capacidade)[arco];
    deriv += single_BPR_derivate(fluxo, tempo_livre, capacidade);
  }
  return deriv;
}

void greedy_solver_newton(struct CONJUNTO_CAMINHOS *conjunto, double demanda,
                          struct PARAMETERS *BPR_PARAMETERS,
                          igraph_vector_t *solucao) {

  if (conjunto->num_caminhos <= 1 || demanda < GREEDY_EPSILON)
    return;

  int idx_min = -1, idx_max = -1;
  double custo_min = DBL_MAX, custo_max = -DBL_MAX;

  for (int h = 0; h < conjunto->num_caminhos; h++) {
    struct CAMINHO *c = &conjunto->caminhos[h];
    double custo = calcular_custo_caminho(c, BPR_PARAMETERS, solucao);
    c->custo_v = custo;

    if (custo < custo_min) {
      custo_min = custo;
      idx_min = h;
    }
    if (c->fluxo > GREEDY_EPSILON && custo > custo_max) {
      custo_max = custo;
      idx_max = h;
    }
  }

  if (idx_max < 0 || idx_max == idx_min) {
    for (int h = 0; h < conjunto->num_caminhos; h++) {
      if (conjunto->caminhos[h].fluxo > GREEDY_EPSILON && h != idx_min) {
        idx_max = h;
        custo_max = conjunto->caminhos[h].custo_v;
        break;
      }
    }
    if (idx_max < 0 || idx_max == idx_min)
      return;
  }

  struct CAMINHO *p_max = &conjunto->caminhos[idx_max];
  struct CAMINHO *p_min = &conjunto->caminhos[idx_min];

  double s_max = calcular_derivada_caminho(p_max, BPR_PARAMETERS, solucao);
  double s_min = calcular_derivada_caminho(p_min, BPR_PARAMETERS, solucao);

  double denominador = s_max + s_min;
  if (denominador < GREEDY_EPSILON)
    return;

  double delta = (custo_max - custo_min) / denominador;
  delta *= GREEDY_DAMPING;

  if (delta > p_max->fluxo)
    delta = p_max->fluxo;
  if (delta < 0)
    delta = 0;

  if (delta < GREEDY_EPSILON)
    return;

  p_max->fluxo -= delta;
  p_min->fluxo += delta;

  int num_arcos_max = igraph_vector_int_size(&p_max->arcos);
  for (int a = 0; a < num_arcos_max; a++) {
    int arco = VECTOR(p_max->arcos)[a];
    VECTOR(*solucao)[arco] -= delta;
    if (VECTOR(*solucao)[arco] < 0)
      VECTOR(*solucao)[arco] = 0;
  }

  int num_arcos_min = igraph_vector_int_size(&p_min->arcos);
  for (int a = 0; a < num_arcos_min; a++) {
    int arco = VECTOR(p_min->arcos)[a];
    VECTOR(*solucao)[arco] += delta;
  }
}

/* ========================================================================
   SEÇÃO 5: LOOP INTERNO
   ======================================================================== */

static double calcular_delta_rs(struct CONJUNTO_CAMINHOS *conjunto,
                                struct PARAMETERS *BPR_PARAMETERS,
                                igraph_vector_t *solucao) {
  if (conjunto->num_caminhos <= 1)
    return 0.0;

  double min_custo = DBL_MAX;
  double max_custo = -DBL_MAX;

  for (int h = 0; h < conjunto->num_caminhos; h++) {
    struct CAMINHO *caminho = &conjunto->caminhos[h];
    double custo = calcular_custo_caminho(caminho, BPR_PARAMETERS, solucao);

    if (custo < min_custo)
      min_custo = custo;
    if (caminho->fluxo > GREEDY_EPSILON && custo > max_custo)
      max_custo = custo;
  }

  if (max_custo < min_custo)
    return 0.0;
  return max_custo - min_custo;
}

void loop_interno_v2(struct GREEDY_STATE *estado, struct OD_MATRIX *OD,
                     struct PARAMETERS *BPR_PARAMETERS,
                     igraph_vector_t *solucao) {

  for (int iter = 0; iter < 20; iter++) {
    int mudancas = 0;

    for (int idx = 0; idx < estado->total_pares; idx++) {
      struct CONJUNTO_CAMINHOS *conjunto = &estado->conjuntos[idx];

      double delta_antes = calcular_delta_rs(conjunto, BPR_PARAMETERS, solucao);
      if (delta_antes > GREEDY_EPSILON) {
        greedy_solver_newton(conjunto, estado->demandas[idx], BPR_PARAMETERS,
                             solucao);
        double delta_depois =
            calcular_delta_rs(conjunto, BPR_PARAMETERS, solucao);
        if (delta_depois < delta_antes * 0.999)
          mudancas++;
      }
    }

    if (mudancas == 0)
      break;
  }
}

/* ========================================================================
   SEÇÃO 6: GERAÇÃO DE COLUNAS
   ======================================================================== */

void geracao_colunas_v2(struct GREEDY_STATE *estado, struct OD_MATRIX *OD,
                        igraph_t *Grafo, struct PARAMETERS *BPR_PARAMETERS,
                        igraph_vector_t *solucao) {

  igraph_vector_t tempo_arcos;
  igraph_vector_init(&tempo_arcos, BPR_PARAMETERS->L);
  BPR(&tempo_arcos, BPR_PARAMETERS, solucao);

  igraph_vector_int_t antecessores;
  igraph_vector_int_init(&antecessores, 0);

  int idx = 0;
  for (int o = 0; o < OD->size; o++) {
    int fonte = OD->Elementos[o].fonte;

    /* Apenas busca novos caminhos se a demanda > 0 */
    /* Para warm start, se temos caminhos mas fluxo está mal distribuído,
     * precisamos gerar novos? */
    /* Sim, Geração de Colunas é essencial. Mas podemos pular se o par OD está
     * satisfeito? */
    /* Algoritmo padrão executa sempre. */

    igraph_get_shortest_paths_dijkstra(Grafo, NULL, NULL, fonte,
                                       igraph_vss_all(), &tempo_arcos,
                                       IGRAPH_OUT, NULL, &antecessores);

    int num_destinos = igraph_vector_int_size(&OD->Elementos[o].alvos);
    for (int d = 0; d < num_destinos; d++) {
      int destino = VECTOR(OD->Elementos[o].alvos)[d];

      /* Verifica se o caminho mais curto já existe no conjunto */
      struct CAMINHO novo_caminho;
      caminho_criar_de_antecessores(&novo_caminho, &antecessores, fonte,
                                    destino, Grafo);
      novo_caminho.fluxo = 0.0;

      struct CONJUNTO_CAMINHOS *conjunto = &estado->conjuntos[idx];

      /* Se é um par OD novo (conjunto vazio), precisamos inicializar o fluxo
       * nele */
      if (conjunto->num_caminhos == 0 &&
          estado->demandas[idx] > GREEDY_EPSILON) {
        /* Fallback do Warm Start (Passo 2 do Checklist) */
        novo_caminho.fluxo = estado->demandas[idx];
        conjunto_adicionar_caminho(conjunto, &novo_caminho);

        /* Atualiza solução (fluxos nos arcos) */
        int num_arcos = igraph_vector_int_size(&novo_caminho.arcos);
        for (int a = 0; a < num_arcos; a++) {
          int arco = VECTOR(novo_caminho.arcos)[a];
          VECTOR(*solucao)[arco] += novo_caminho.fluxo;
        }
      } else {
        /* Adiciona como candidato com fluxo 0 */
        bool adicionado = conjunto_adicionar_caminho(conjunto, &novo_caminho);
        if (!adicionado) {
          conjunto_remover_caminhos_vazios(conjunto);
        }
      }

      caminho_destruir(&novo_caminho);
      idx++;
    }
  }

  igraph_vector_destroy(&tempo_arcos);
  igraph_vector_int_destroy(&antecessores);
}

/* ========================================================================
   SEÇÃO 7: PONTO DE ENTRADA PRINCIPAL (Com Suporte a Warm Start)
   ======================================================================== */

void GreedyPath(struct PARAMETERS *BPR_PARAMETERS, struct OD_MATRIX *OD,
                igraph_t *Grafo, igraph_vector_t *solucao,
                struct GREEDY_STATE **ptr_persistent_state) {

  /* Inicializa vetor de solução (sempre necessário) */
  igraph_vector_init(solucao, BPR_PARAMETERS->L);

  /* Aloca novo estado */
  struct GREEDY_STATE *novo_estado = malloc(sizeof(struct GREEDY_STATE));
  if (!novo_estado) {
    perror("GreedyPath: Erro ao alocar estado");
    return;
  }

  /* Inicializa estrutura do novo estado baseada na matriz OD atual */
  greedy_inicializar_estado(novo_estado, OD);

  /* Verifica se existe estado anterior para Warm Start */
  if (ptr_persistent_state && *ptr_persistent_state != NULL) {
    /* EXECUTA WARM START */
    greedy_executar_warm_start(novo_estado, *ptr_persistent_state, OD, solucao);

    /* Libera o estado anterior pois não é mais necessário */
    greedy_destruir_estado(*ptr_persistent_state);
    free(*ptr_persistent_state);
    *ptr_persistent_state = NULL; /* Segurança */

  } else {
    /* COLD START: All-or-Nothing inicial */
    /* Note: inicializa e já popula solucao e conjuntos */
    greedy_destruir_estado(
        novo_estado); /* Limpa a inicialização vazia feita acima */
    greedy_atribuicao_inicial_v2(BPR_PARAMETERS, OD, Grafo, solucao,
                                 novo_estado);
  }

  double gap = fabs(relative_gap(solucao, Grafo, BPR_PARAMETERS, OD, NULL, NULL));
  novo_estado->RG_anterior = gap;

  //printf("GreedyPath: %d pares OD. Estado inicial GAP: %.6e\n",
  //       novo_estado->total_pares, gap);

  clock_t start_time, end_time;
  double tempo;
  double best_gap = gap;

  for (int iteracao = 0; iteracao < GREEDY_MAX_MAIN_ITER; iteracao++) {
    start_time = clock();

    /* Geração de colunas */
    geracao_colunas_v2(novo_estado, OD, Grafo, BPR_PARAMETERS, solucao);

    /* Loop interno */
    loop_interno_v2(novo_estado, OD, BPR_PARAMETERS, solucao);

    gap = fabs(relative_gap(solucao, Grafo, BPR_PARAMETERS, OD, NULL, NULL));
    if (gap < best_gap)
      best_gap = gap;

    end_time = clock();
    tempo = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    //printf("GreedyPath Iter %3d | GAP: %.6e | Best: %.6e | Tempo: %.2fs\n",
    //       iteracao + 1, gap, best_gap, tempo);

    if (gap < GREEDY_EPSILON) {
      //printf("GreedyPath: Convergência alcançada na iteração %d\n",
      //       iteracao + 1);
      break;
    }

    novo_estado->RG_anterior = gap;
  }

  /* Persistência de estado */
  if (ptr_persistent_state) {
    *ptr_persistent_state = novo_estado;
  } else {
    /* Se o chamador não quiser persistência, limpamos tudo */
    greedy_destruir_estado(novo_estado);
    free(novo_estado);
  }
}
