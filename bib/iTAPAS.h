#pragma once
#include "calc.h"
#include "define.h"
#include "mtwister.h"
#include <igraph/igraph.h>

void pas_preparar_assinatura(struct PAS *pas) {
  igraph_vector_int_init_copy(&pas->s1_sorted, &pas->c1);
  igraph_vector_int_sort(&pas->s1_sorted);

  igraph_vector_int_init_copy(&pas->s2_sorted, &pas->c2);
  igraph_vector_int_sort(&pas->s2_sorted);
}

static inline void initialize_warm_start(struct WARM_START *WS, int size) {
  WS->has_warm_start = false;
  WS->INCREMENT = 0.0;
  WS->num_pas = 0;
  WS->conjunto_pas = NULL;

  WS->SPT = (struct SPT *)malloc(size * sizeof(struct SPT));
  if (!WS->SPT) {
    perror("Initialization error: Failed to allocate memory for SPT");
    exit(EXIT_FAILURE);
  }

  // É fundamental inicializar os vetores dentro das estruturas SPT
  for (int i = 0; i < size; i++) {
    igraph_vector_int_init(&WS->SPT[i].antecessores, 0);
    igraph_vector_init(&WS->SPT[i].dist, 0);
  }

  // Inicialização do grafo bush como vazio
  igraph_empty(&WS->grafo_bush, 0, IGRAPH_DIRECTED);

  // Inicialização dos vetores globais de estado
  igraph_vector_init(&WS->flow, 0);
  igraph_vector_init(&WS->time, 0);
}
/**
 * @brief Calcula e armazena a forma canônica de um PAS.
 * Modifica o PAS in-place, preenchendo canon_s1 e canon_s2.
 */

/* ========================================================================
   FUNÇÕES AUXILIARES PARA MANIPULAÇÃO DE PAS
   ======================================================================== */

void pas_add_origin(struct PAS *pas, int origem) {
  if (pas->num_origins >= pas->capacity_origins) {
    int new_cap = pas->capacity_origins == 0 ? 4 : pas->capacity_origins * 2;
    struct PAS_Origin *tmp = (struct PAS_Origin *)realloc(
        pas->origins, new_cap * sizeof(struct PAS_Origin));
    if (!tmp) {
      exit(EXIT_FAILURE);
    }
    pas->origins = tmp;
    pas->capacity_origins = new_cap;
  }
  pas->origins[pas->num_origins].origem = origem;
  pas->origins[pas->num_origins].flow_s1 = 0.0;
  pas->origins[pas->num_origins].flow_s2 = 0.0;
  pas->origins[pas->num_origins].shift = 0.0;
  pas->num_origins++;
}

void pas_print(const struct PAS *pas, igraph_t *Grafo) {
  if (!pas)
    return;

  printf("\nPAS (Shared by %d origins):", pas->num_origins);
  for (int i = 0; i < pas->num_origins; i++) {
    printf(" %d", pas->origins[i].origem);
  }
  printf("\n");

  int tamanho_c1 = igraph_vector_int_size(&pas->c1);
  printf("  c1 (size=%d):", tamanho_c1);
  for (int i = tamanho_c1 - 1; i >= 0; i--) {
    int arco = VECTOR(pas->c1)[i];
    int from = IGRAPH_FROM(Grafo, arco) + 1;
    int to = IGRAPH_TO(Grafo, arco) + 1;
    printf(" (%d, %d)", from, to);
  }

  printf("\n  c2 (size=%ld):", igraph_vector_int_size(&pas->c2));
  for (int i = igraph_vector_int_size(&pas->c2) - 1; i >= 0; i--) {
    int arco = VECTOR(pas->c2)[i];
    int from = IGRAPH_FROM(Grafo, arco) + 1;
    int to = IGRAPH_TO(Grafo, arco) + 1;
    printf(" (%d, %d)", from, to);
  }
  printf("\n\n");
}

void pas_free(struct PAS *pas) {
  if (!pas)
    return;
  igraph_vector_int_destroy(&pas->c1);
  igraph_vector_int_destroy(&pas->c2);
  igraph_vector_int_destroy(&pas->s1_sorted);
  igraph_vector_int_destroy(&pas->s2_sorted);
  if (pas->origins) {
    free(pas->origins);
    pas->origins = NULL;
  }
  pas->num_origins = 0;
  pas->capacity_origins = 0;
}

void pas_remove_at(struct PAS **vec, int *tamanho, int indice) {
  if (!vec || !(*vec) || !tamanho || indice < 0 || indice >= *tamanho) {
    return;
  }

  // 1. Libera os recursos do elemento que será removido.
  pas_free(&(*vec)[indice]);

  // 2. Desloca as estruturas restantes para cobrir o buraco.
  //    Este laço faz cópias superficiais, movendo os ponteiros.
  for (int i = indice; i < (*tamanho - 1); ++i) {
    (*vec)[i] = (*vec)[i + 1];
  }

  // 3. Diminui o tamanho lógico.
  int novo_tamanho = *tamanho - 1;

  if (novo_tamanho > 0) {
    // 4. Tenta encolher o array. Isso remove a cópia duplicada no final.
    struct PAS *tmp =
        (struct PAS *)realloc(*vec, novo_tamanho * sizeof(struct PAS));
    // Apenas atualiza o ponteiro se a realocação for bem-sucedida.
    if (tmp != NULL) {
      *vec = tmp;
    } else {
      // Se realloc falhar, o ponteiro original ainda é válido.
      // O programa pode continuar, mas haverá um vazamento de memória
      // da estrutura extra no final do array. É melhor que um crash.
      fprintf(stderr, "Aviso: realloc falhou em pas_remove_at.\n");
    }
  } else {
    // O array está vazio, libera completamente.
    free(*vec);
    *vec = NULL;
  }

  // 5. Atualiza o contador de tamanho.
  *tamanho = novo_tamanho;
}

/* ========================================================================
   CÁLCULO DE CUSTOS E DERIVADAS
   ======================================================================== */

static void
calcular_custos_caminho(igraph_t *Grafo, struct PARAMETERS *BPR_PARAMETERS,
                        igraph_vector_int_t *caminho, igraph_vector_t solucao,
                        const char *origem_attr, double *tempo_total,
                        double *derivada_total, double *capacidade_minima) {
  *tempo_total = 0.0;
  *derivada_total = 0.0;
  *capacidade_minima = 1e10;

  int num_arcos = igraph_vector_int_size(caminho);
  for (int i = 0; i < num_arcos; i++) {
    int arco_id = VECTOR(*caminho)[i];
    double fluxo = VECTOR(solucao)[arco_id];
    double tempo_livre = VECTOR(BPR_PARAMETERS->cost_time)[arco_id];
    double capacidade = VECTOR(BPR_PARAMETERS->capacidade)[arco_id];

    *tempo_total += single_BPR(fluxo, tempo_livre, capacidade);
    *derivada_total += single_BPR_derivate(fluxo, tempo_livre, capacidade);
    *capacidade_minima =
        fmin(*capacidade_minima, EAN(Grafo, origem_attr, arco_id));
  }
}

/* ========================================================================
   DESLOCAMENTO DE FLUXO
   ======================================================================== */

double deslocar_fluxo_no_pas(igraph_t *Grafo, struct PARAMETERS *BPR_PARAMETERS,
                             struct PAS *pas, igraph_vector_t solucao) {
  // 1. Calcula custos globais e derivadas (baseado em qualquer origem, pois
  // custo é link-based) Mas para capacidade, precisamos checar se todas as
  // origens tem fluxo. Usaremos a primeira origem para cálculo de custos (custo
  // é agnóstico da origem)

  if (pas->num_origins == 0)
    return 0.0;

  char origem_attr_dummy[20];
  sprintf(origem_attr_dummy, "demanda_%d", pas->origins[0].origem);

  double tempo_c1, tempo_c2, derivada_c1, derivada_c2;
  double dummy_cap;

  // Calcula custos (sem se preocupar com capacidade específica de origem ainda)
  calcular_custos_caminho(Grafo, BPR_PARAMETERS, &pas->c1, solucao,
                          origem_attr_dummy, &tempo_c1, &derivada_c1,
                          &dummy_cap);

  calcular_custos_caminho(Grafo, BPR_PARAMETERS, &pas->c2, solucao,
                          origem_attr_dummy, &tempo_c2, &derivada_c2,
                          &dummy_cap);

  // Verificação de custo
  if (tempo_c2 <= tempo_c1 + TAPAS_TOLERANCE) {
    return 0.0;
  }

  // 2. Calcula Shift Total Potencial (Newton step)
  double denominador = derivada_c1 + derivada_c2;
  double total_shift_desired = 0.0;

  if (denominador < 1e-12) {
    total_shift_desired = 1e10;
  } else {
    total_shift_desired = (tempo_c2 - tempo_c1) / denominador;
  }

  // 3. Proporcionalidade: Calcular contribuição de cada origem
  double total_min_flow_on_s2 = 0.0; // Capacidade total de redução no s2

  char origem_attr[20];

  // Loop para calcular fluxos disponíveis por origem
  for (int i = 0; i < pas->num_origins; i++) {
    int org = pas->origins[i].origem;
    sprintf(origem_attr, "demanda_%d", org);

    // Calcula min flow desta origem em s2 (caminho a reduzir)
    double min_f_s2 = 1e10;
    int size_c2 = igraph_vector_int_size(&pas->c2);
    for (int k = 0; k < size_c2; k++) {
      int arco = VECTOR(pas->c2)[k];
      double f = EAN(Grafo, origem_attr, arco);
      min_f_s2 = fmin(min_f_s2, f);
    }

    // Salva no struct para usar depois
    pas->origins[i].flow_s2 = min_f_s2;
    total_min_flow_on_s2 += min_f_s2;
  }

  // 4. Limita o shift total pela capacidade total disponível
  double actual_total_shift = fmin(total_shift_desired, total_min_flow_on_s2);

  if (actual_total_shift < TAPAS_TOLERANCE) {
    return 0.0;
  }

  // 5. Aplica Shifting Proporcional
  // Se total_min_flow_on_s2 é muito pequeno, evitamos divisão por zero
  if (total_min_flow_on_s2 < TAPAS_TOLERANCE)
    return 0.0;

  double total_applied = 0.0;

  for (int i = 0; i < pas->num_origins; i++) {
    // Proporção: quanto esta origem contribui para o fluxo total dispensável
    double origin_share = pas->origins[i].flow_s2 / total_min_flow_on_s2;
    double origin_shift = actual_total_shift * origin_share;

    // Safety check
    origin_shift =
        fmin(origin_shift, pas->origins[i].flow_s2); // Nunca tirar mais que tem

    if (origin_shift < TAPAS_FLOW_TOL)
      continue;

    pas->origins[i].shift = origin_shift;
    total_applied += origin_shift;

    sprintf(origem_attr, "demanda_%d", pas->origins[i].origem);

    // Adiciona em s1
    int size_c1 = igraph_vector_int_size(&pas->c1);
    for (int k = 0; k < size_c1; k++) {
      int arco = VECTOR(pas->c1)[k];
      // Atualiza global
      VECTOR(solucao)[arco] += origin_shift;
      // Atualiza atributo da origem
      double custo_atual = EAN(Grafo, origem_attr, arco);
      igraph_cattribute_EAN_set(Grafo, origem_attr, arco,
                                custo_atual + origin_shift);
    }

    // Remove de s2
    int size_c2 = igraph_vector_int_size(&pas->c2);
    for (int k = 0; k < size_c2; k++) {
      int arco = VECTOR(pas->c2)[k];
      // Atualiza global
      VECTOR(solucao)[arco] -= origin_shift;
      // Atualiza atributo da origem
      double custo_atual = EAN(Grafo, origem_attr, arco);
      double novo_custo = custo_atual - origin_shift;
      if (novo_custo < TAPAS_COST_TOL)
        novo_custo = 0.0;
      igraph_cattribute_EAN_set(Grafo, origem_attr, arco, novo_custo);
    }
  }

  return total_applied;
}

void deslocamento_global_pas(struct PARAMETERS *BPR_PARAMETERS,
                             struct OD_MATRIX *OD, struct WARM_START *WS,
                             igraph_vector_t *solucao) {

  if (WS->num_pas <= 0) {
    // printf("No PAS to process.\n");
    return;
  }
  int n;
  for (n = 0; n < 200; n++) {
    // pas_print(&conjunto_pas[n],Grafo);
    for (int i = 0; i < WS->num_pas; i++) {
      struct PAS *pas_atual = &WS->conjunto_pas[i];

      // Pruning simplificado: Se não tem origens, remove? (Na prática não deve
      // acontecer)
      if (pas_atual->num_origins == 0) {
        pas_remove_at(&WS->conjunto_pas, &WS->num_pas, i);
        i--;
        continue;
      }

      double dx = deslocar_fluxo_no_pas(&WS->grafo_bush, BPR_PARAMETERS,
                                        pas_atual, *solucao);

      // Se dx for muito pequeno repetidamente, poderíamos pensar em remover,
      // mas como são múltiplas origens, melhor manter se pelo menos uma tiver
      // fluxo. O pruning ideal verificaria se total_min_flow_on_s2 < TOL.
    }
  }
}

/* ========================================================================
   IDENTIFICAÇÃO DE PAS
   ======================================================================== */

void identificar_pas_fluxo_maximo(igraph_t *Grafo, int arco_inicial, int fonte,
                                  igraph_vector_int_t *antecessores,
                                  struct PAS *pas) {
  char origem_attr[20];
  sprintf(origem_attr, "demanda_%d", fonte);

  // pas->origem = fonte; REMOVED
  igraph_vector_int_init(&pas->c1, 0);
  igraph_vector_int_init(&pas->c2, 0);
  // Inicialize sorted como vazios por segurança, serão preenchidos depois
  igraph_vector_int_init(&pas->s1_sorted, 0);
  igraph_vector_int_init(&pas->s2_sorted, 0);

  pas->origins = NULL;
  pas->num_origins = 0;
  pas->capacity_origins = 0;

  // Configurações iniciais
  int cabeca = IGRAPH_TO(Grafo, arco_inicial);
  int cauda_arco_deseq = IGRAPH_FROM(Grafo, arco_inicial);
  int num_nos = igraph_vcount(Grafo);

  // Arrays auxiliares para busca O(1)
  bool *nos_em_s1 = (bool *)calloc(num_nos, sizeof(bool));
  if (!nos_em_s1)
    exit(EXIT_FAILURE); // Em HPC, trate o erro adequadamente

  // --- 1. Construção do Segmento 1 (Árvore SPT) ---
  // Rastreia da cabeça para trás até a fonte ou fim da árvore
  int no_atual = cabeca;
  nos_em_s1[no_atual] = true;

  while (no_atual != fonte) {
    int arco = VECTOR(*antecessores)[no_atual];
    if (arco == -1)
      break; // Chegou na raiz ou desconexo

    igraph_vector_int_push_back(&pas->c1, arco);

    no_atual = IGRAPH_FROM(Grafo, arco);
    if (nos_em_s1[no_atual])
      break; // Evita loop na própria árvore
    nos_em_s1[no_atual] = true;
  }

  // --- 2. Construção do Segmento 2 (Backtracking por Fluxo Máximo) ---
  // Do arco desequilibrado para trás até encontrar ALGUÉM de s1

  igraph_vector_int_push_back(&pas->c2, arco_inicial);

  no_atual = cauda_arco_deseq;
  bool encontrou_divergencia = false;
  int no_divergencia = -1;

  bool *visitado_s2 = (bool *)calloc(num_nos, sizeof(bool));
  visitado_s2[no_atual] = true;

  // Verifica se o nó inicial já é a divergência
  if (nos_em_s1[no_atual]) {
    encontrou_divergencia = true;
    no_divergencia = no_atual;
  }

  igraph_vector_int_t arcos_entrantes;
  igraph_vector_int_init(&arcos_entrantes, 0);

  while (!encontrou_divergencia) {
    igraph_incident(Grafo, &arcos_entrantes, no_atual, IGRAPH_IN,
                    IGRAPH_NO_LOOPS);
    int n_incidentes = igraph_vector_int_size(&arcos_entrantes);

    if (n_incidentes == 0)
      break; // Beco sem saída

    int melhor_arco = -1;
    double max_f = -1.0;
    int melhor_pred = -1;

    // Encontra predecessor com maior fluxo para esta origem
    for (int i = 0; i < n_incidentes; i++) {
      int arco = VECTOR(arcos_entrantes)[i];
      double f = EAN(Grafo, origem_attr, arco);
      if (f > max_f) {
        max_f = f;
        melhor_arco = arco;
        melhor_pred = IGRAPH_FROM(Grafo, arco);
      }
    }
    // Se não há fluxo ou arco válido, aborta
    if (melhor_arco == -1 || max_f <= TAPAS_TOLERANCE)
      break;

    // --- TRATAMENTO DE CICLO (CONSISTÊNCIA COM PYTHON) ---
    if (visitado_s2[melhor_pred]) {
      // Ciclo detectado em s2. O Python retorna None.
      // Aqui limpamos os vetores para indicar PAS inválido.
      igraph_vector_int_clear(&pas->c1);
      igraph_vector_int_clear(&pas->c2);
      break; // Sai do loop e retornará vetores vazios
    }

    igraph_vector_int_push_back(&pas->c2, melhor_arco);
    visitado_s2[melhor_pred] = true;
    no_atual = melhor_pred;

    if (nos_em_s1[no_atual]) {
      encontrou_divergencia = true;
      no_divergencia = no_atual;
    }
  }

  igraph_vector_int_destroy(&arcos_entrantes);
  free(visitado_s2);
  free(nos_em_s1);

  // Se falhou por ciclo ou beco sem saída, garante retorno vazio
  if (!encontrou_divergencia || igraph_vector_int_size(&pas->c2) == 0) {
    igraph_vector_int_clear(&pas->c1);
    igraph_vector_int_clear(&pas->c2);
    return;
  }

  // --- 3. Truncar s1 (Cabeça -> Divergência) ---
  // Reconstrói s1 para conter apenas o trecho da divergência até a cabeça
  // (A implementação anterior continha o caminho todo até a raiz)

  igraph_vector_int_clear(&pas->c1);
  int temp_node = cabeca;

  // Segurança contra loops infinitos na reconstrução
  int limit = 0;
  int max_steps = igraph_vcount(Grafo);

  while (temp_node != no_divergencia && limit++ < max_steps) {
    int arco = VECTOR(*antecessores)[temp_node];
    if (arco == -1)
      break;
    igraph_vector_int_push_back(&pas->c1, arco);
    temp_node = IGRAPH_FROM(Grafo, arco);
  }
}

/* ========================================================================
   PROCESSAMENTO POR ORIGEM
   ======================================================================== */

void processar_origem(struct PARAMETERS *BPR_PARAMETERS, struct OD_MATRIX *OD,
                      struct WARM_START *WS, int indice_origem,
                      igraph_vector_t *solucao, int iteracao, double gap) {

  int fonte = OD->Elementos[indice_origem].fonte;
  char origem_attr[20];
  sprintf(origem_attr, "demanda_%d", fonte);

  igraph_vector_t tempo_arcos, custo_nos;
  igraph_vector_int_t arcos_removidos, antecessores;

  igraph_vector_init(&tempo_arcos, BPR_PARAMETERS->L);
  igraph_vector_init(&custo_nos, BPR_PARAMETERS->N);
  igraph_vector_int_init(&arcos_removidos, 0);
  igraph_vector_int_init(&antecessores, 0);

  // 1. Calcula Caminhos Mínimos (SPT)
  BPR(&tempo_arcos, BPR_PARAMETERS, solucao);
  // clock_t start_time = clock();
  igraph_get_shortest_paths_dijkstra(&WS->grafo_bush, NULL, NULL, fonte,
                                     igraph_vss_all(), &tempo_arcos, IGRAPH_OUT,
                                     NULL, &antecessores);
  // clock_t end_time = clock();
  // double elapsed_secs = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  // printf("    Tempos: %.2e ", elapsed_secs);

  igraph_vector_int_update(&WS->SPT[indice_origem].antecessores, &antecessores);

  // 2. Calcula Potenciais dos Nós
  for (int no = 0; no < BPR_PARAMETERS->N; no++) {
    VECTOR(custo_nos)[no] = 0.0;
    int no_atual = no;
    // Reconstrói custo somando arcos (mais preciso que usar dist do Dijkstra se
    // houver turn penalties futuros)
    while (no_atual != fonte) {
      int arco = VECTOR(antecessores)[no_atual];
      if (arco == -1)
        break;
      VECTOR(custo_nos)[no] += VECTOR(tempo_arcos)[arco];
      no_atual = IGRAPH_FROM(&WS->grafo_bush, arco);
    }
  }

  // Salva distâncias para uso no GAP aproximado do warm start
  igraph_vector_update(&WS->SPT[indice_origem].dist, &custo_nos);

  // 3. Identifica Arcos Fora da Árvore (Desequilibrados)
  // start_time = clock();
  for (int arco = 0; arco < BPR_PARAMETERS->L; arco++) {
    double fluxo = EAN(&WS->grafo_bush, origem_attr, arco);
    if (fluxo <= TAPAS_TOLERANCE) // Check vs epsilon
      continue;

    int no_origem = IGRAPH_FROM(&WS->grafo_bush, arco);
    int no_destino = IGRAPH_TO(&WS->grafo_bush, arco);

    // Custo Reduzido: c_uv + pi_u - pi_v
    double custo_reduzido = VECTOR(tempo_arcos)[arco] +
                            VECTOR(custo_nos)[no_origem] -
                            VECTOR(custo_nos)[no_destino];

    // Limiar adaptativo (ref: iTAPAS GeneratePAS - Xie)
    // Iterações iniciais: foco nos desequilíbrios grandes
    // Conforme converge: limiar fica mais fino
    double pre;
    if (iteracao <= 1)
      pre = 1e-1;
    else if (iteracao == 2)
      pre = 1e-3;
    else
      pre = (gap > TAPAS_TOLERANCE) ? gap / 100.0 : TAPAS_THETA;

    if (custo_reduzido > pre && fluxo > TAPAS_FLOW_TOL) {
      // if(fonte+1 == 14013)printf("Arco (%d -> %d) é desequilibrado.
      // Fluxo=%.6f, Custo Reduzido=%.6f\n", no_origem+1, no_destino+1, fluxo,
      // custo_reduzido);
      igraph_vector_int_push_back(&arcos_removidos, arco);
    }
  }
  // end_time = clock();
  // elapsed_secs = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  // printf("%.2e ", elapsed_secs);

  // 4. Processa Candidatos a PAS
  int num_candidatos = igraph_vector_int_size(&arcos_removidos);
  // printf("Origem %d: %d arcos desequilibrados encontrados.\n", fonte+1,
  // num_candidatos); start_time = clock();
  for (int i = 0; i < num_candidatos; i++) {
    struct PAS novo_pas;
    int arco = VECTOR(arcos_removidos)[i];
    // printf("%ld %ld\n", IGRAPH_FROM(&WS->grafo_bush, arco) + 1,
    // IGRAPH_TO(&WS->grafo_bush, arco) + 1);
    identificar_pas_fluxo_maximo(&WS->grafo_bush, arco, fonte, &antecessores,
                                 &novo_pas);

    // Verifica se PAS é válido (retorno não vazio)
    if (igraph_vector_int_size(&novo_pas.c1) == 0 &&
        igraph_vector_int_size(&novo_pas.c2) == 0) {
      // if(IGRAPH_FROM(&WS->grafo_bush, arco) + 1 == 5658)printf("    PAS
      // inválido (ciclo ou beco sem saída). Ignorando. %ld
      // %ld\n",igraph_vector_int_size(&novo_pas.c1),
      // igraph_vector_int_size(&novo_pas.c2));
      pas_free(&novo_pas);
      continue;
    }

    // Prepara assinatura para verificação de unicidade
    pas_preparar_assinatura(&novo_pas);

    struct PAS *pas_target = NULL;
    bool pas_existe = false;

    // Verificação De Unicidade (Topologia)
    for (int j = 0; j < WS->num_pas; j++) {
      struct PAS *existente = &WS->conjunto_pas[j];

      // Compara s1 (Árvore)
      bool s1_igual =
          igraph_vector_int_all_e(&existente->s1_sorted, &novo_pas.s1_sorted);
      if (!s1_igual)
        continue;

      // Compara s2 (Atalho)
      bool s2_igual =
          igraph_vector_int_all_e(&existente->s2_sorted, &novo_pas.s2_sorted);
      if (!s2_igual)
        continue;

      // Se chegou aqui, é idêntico em topologia
      pas_target = existente;
      pas_existe = true;
      break;
    }

    if (pas_existe) {
      // Adiciona a origem ao PAS existente (se ainda não estiver lá, mas nosso
      // código permite duplicatas de origem se não checar) Vamos assumir que
      // processar_origem é chamado uma vez por iteracao por origem, então ok.
      // Precisamos checar se origem já está na lista? `pas_add_origin` não
      // checa. O ideal é checar.
      bool org_already_in = false;
      for (int k = 0; k < pas_target->num_origins; k++) {
        if (pas_target->origins[k].origem == fonte) {
          org_already_in = true;
          break;
        }
      }
      if (!org_already_in) {
        pas_add_origin(pas_target, fonte);
      }

      // Tenta deslocar (agora shift proporcional para TODAS origens nele,
      // inclusive a nova)
      deslocar_fluxo_no_pas(&WS->grafo_bush, BPR_PARAMETERS, pas_target,
                            *solucao);

      pas_free(&novo_pas); // Não precisamos mais da cópia temporária
    } else {
      // Cria novo PAS no conjunto global
      struct PAS *tmp = (struct PAS *)realloc(
          WS->conjunto_pas, (WS->num_pas + 1) * sizeof(struct PAS));
      if (tmp == NULL) {
        perror("Falha crítica de memória ao expandir conjunto PAS");
        pas_free(&novo_pas);
        exit(EXIT_FAILURE);
      }
      WS->conjunto_pas = tmp;

      struct PAS *destino = &WS->conjunto_pas[WS->num_pas];

      // Copia topologia e assinaturas
      igraph_vector_int_init_copy(&destino->c1, &novo_pas.c1);
      igraph_vector_int_init_copy(&destino->c2, &novo_pas.c2);
      igraph_vector_int_init_copy(&destino->s1_sorted, &novo_pas.s1_sorted);
      igraph_vector_int_init_copy(&destino->s2_sorted, &novo_pas.s2_sorted);

      // Inicializa vetor de origens
      destino->origins = NULL;
      destino->num_origins = 0;
      destino->capacity_origins = 0;

      pas_add_origin(destino, fonte);
      WS->num_pas++;

      // Shift inicial
      deslocar_fluxo_no_pas(&WS->grafo_bush, BPR_PARAMETERS, destino, *solucao);

      pas_free(&novo_pas); // Libera vetores da struct temporária (copiamos o
                           // conteúdo)
    }
  }

  // end_time = clock();
  // elapsed_secs = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  // printf(" %d ",*num_pas);
  // printf("%.2e\n", elapsed_secs);
  // printf("    Após processamento da origem %d, total de PAS: %d\n", fonte+1,
  // *num_pas);
  //  Limpeza local
  igraph_vector_destroy(&tempo_arcos);
  igraph_vector_destroy(&custo_nos);
  igraph_vector_int_destroy(&arcos_removidos);
  igraph_vector_int_destroy(&antecessores);
}

/* ========================================================================
   ATRIBUIÇÃO INICIAL (ALL-OR-NOTHING)
   ======================================================================== */

void atribuicao_inicial(struct PARAMETERS *BPR_PARAMETERS, struct OD_MATRIX *OD,
                        igraph_t *Grafo, int indice_origem,
                        igraph_vector_t *solucao) {
  int fonte = OD->Elementos[indice_origem].fonte;

  igraph_vector_int_t antecessores;
  igraph_vector_t demanda_arcos;

  igraph_vector_int_init(&antecessores, 0);
  igraph_vector_init(&demanda_arcos, BPR_PARAMETERS->L);

  /* Calcula caminhos mínimos com tempo livre de fluxo */
  igraph_get_shortest_paths_dijkstra(Grafo, NULL, NULL, fonte, igraph_vss_all(),
                                     &BPR_PARAMETERS->cost_time, IGRAPH_OUT,
                                     NULL, &antecessores);

  /* Atribui demanda aos caminhos mínimos */
  int num_destinos =
      igraph_vector_int_size(&OD->Elementos[indice_origem].alvos);
  for (int i = 0; i < num_destinos; i++) {
    int destino = VECTOR(OD->Elementos[indice_origem].alvos)[i];
    double volume = VECTOR(OD->Elementos[indice_origem].volumes)[i];

    int no_atual = destino;
    while (no_atual != fonte) {
      int arco = VECTOR(antecessores)[no_atual];
      if (arco == -1) {
        fprintf(stderr, "Erro: caminho não encontrado de %d para %d\n", fonte,
                destino);
        exit(EXIT_FAILURE);
      }

      VECTOR(demanda_arcos)[arco] += volume;
      VECTOR(*solucao)[arco] += volume;
      no_atual = IGRAPH_FROM(Grafo, arco);
    }
  }

  char origem_attr[20];
  sprintf(origem_attr, "demanda_%d", fonte);
  igraph_cattribute_EAN_setv(Grafo, origem_attr, &demanda_arcos);

  igraph_vector_int_destroy(&antecessores);
  igraph_vector_destroy(&demanda_arcos);
}

/* ========================================================================
   ALGORITMO PRINCIPAL iTAPAS
   ======================================================================== */

void iTAPAS(struct PARAMETERS *BPR_PARAMETERS, struct OD_MATRIX *OD,
            igraph_t *Grafo, igraph_vector_t *solucao, struct WARM_START *WS) {

  /* Fase 1: Atribuição inicial (all-or-nothing) */
  // clock_t start_time, end_time;
  // start_time = clock();
  int valor = count_files_in_dir("./output/iTAPAS");
  if (!WS->has_warm_start) {
    WS->conjunto_pas = NULL;
    WS->num_pas = 0;
    /* COLD START: Atribuição inicial */
    igraph_vector_init(solucao, BPR_PARAMETERS->L);
    igraph_copy(&WS->grafo_bush, Grafo);
    for (int i = 0; i < OD->size; i++) {
      atribuicao_inicial(BPR_PARAMETERS, OD, &WS->grafo_bush, i, solucao);
    }
  } else {
    /* WARM START: Reconstruir solucao e aplicar incremento proporcional */

    igraph_vector_init(solucao, BPR_PARAMETERS->L);
    igraph_vector_fill(solucao, 0.0);

    // Buffer para leitura/escrita em bulk dos atributos de arco por origem
    igraph_vector_t demanda_origem;
    igraph_vector_init(&demanda_origem, BPR_PARAMETERS->L);

    for (int o = 0; o < OD->size; o++) {
      int fonte = OD->Elementos[o].fonte;
      char origem_attr[20];
      sprintf(origem_attr, "demanda_%d", fonte);

      // Leitura em bulk: 1 lookup em vez de L lookups individuais
      EANV(&WS->grafo_bush, origem_attr, &demanda_origem);

      /* 1. Calcular fluxo total ATUAL saindo desta origem */
      double current_total_flow = 0.0;
      for (int i = 0; i < BPR_PARAMETERS->L; i++) {
        current_total_flow += VECTOR(demanda_origem)[i];
      }

      int num_destinos = igraph_vector_int_size(&OD->Elementos[o].alvos);
      double total_increment = (double)num_destinos * WS->INCREMENT;
      bool can_scale = (current_total_flow > TAPAS_TOLERANCE);

      if (can_scale) {
        // Mantém fluxos antigos intactos na solução
        for (int i = 0; i < BPR_PARAMETERS->L; i++) {
          double old_val = VECTOR(demanda_origem)[i];
          if (old_val > TAPAS_FLOW_TOL) {
            VECTOR(*solucao)[i] += old_val;
          }
        }

        // Atribui demanda nova ao menor caminho anterior (antecessores)
        if (WS->SPT) {
          for (int i = 0; i < num_destinos; i++) {
            int destino = VECTOR(OD->Elementos[o].alvos)[i];
            int no_atual = destino;
            while (no_atual != fonte) {
              if (igraph_vector_int_size(&WS->SPT[o].antecessores) <= no_atual)
                break;
              int arco = VECTOR(WS->SPT[o].antecessores)[no_atual];
              if (arco == -1)
                break;

              VECTOR(*solucao)[arco] += WS->INCREMENT;
              VECTOR(demanda_origem)[arco] += WS->INCREMENT;

              no_atual = IGRAPH_FROM(&WS->grafo_bush, arco);
            }
          }
        }
      } else {
        // Fallback: Se não há fluxo prévio, usa SPT
        if (total_increment > TAPAS_TOLERANCE && WS->SPT) {
          for (int i = 0; i < num_destinos; i++) {
            int destino = VECTOR(OD->Elementos[o].alvos)[i];
            int no_atual = destino;
            while (no_atual != fonte) {
              if (igraph_vector_int_size(&WS->SPT[o].antecessores) <= no_atual)
                break;
              int arco = VECTOR(WS->SPT[o].antecessores)[no_atual];
              if (arco == -1)
                break;

              VECTOR(*solucao)[arco] += WS->INCREMENT;
              VECTOR(demanda_origem)[arco] = WS->INCREMENT;

              no_atual = IGRAPH_FROM(&WS->grafo_bush, arco);
            }
          }
        }
      }
      // Escrita em bulk: 1 operação em vez de múltiplos EAN_set individuais
      igraph_cattribute_EAN_setv(&WS->grafo_bush, origem_attr, &demanda_origem);
    }

    igraph_vector_destroy(&demanda_origem);

    /* 3. GAP aproximado usando distâncias salvas (análogo a dist_shortest_local) */
    igraph_vector_t tempo_arcos_novo;
    igraph_vector_init(&tempo_arcos_novo, BPR_PARAMETERS->L);
    BPR(&tempo_arcos_novo, BPR_PARAMETERS, solucao);

    // TSTT com custos novos
    double total_travel_time = 0.0;
    for (int i = 0; i < BPR_PARAMETERS->L; i++) {
        total_travel_time += VECTOR(*solucao)[i] * VECTOR(tempo_arcos_novo)[i];
    }
    igraph_vector_destroy(&tempo_arcos_novo);

    // SPTT com distâncias ANTIGAS salvas (da última resolução real)
    double total_min_cost = 0.0;
    for (int o = 0; o < OD->size; o++) {
        // Verifica se temos distâncias salvas para esta origem
        if (igraph_vector_size(&WS->SPT[o].dist) < BPR_PARAMETERS->N) {
            total_min_cost = 0.0;
            break;
        }
        int num_destinos = igraph_vector_int_size(&OD->Elementos[o].alvos);
        for (int d = 0; d < num_destinos; d++) {
            int destino = VECTOR(OD->Elementos[o].alvos)[d];
            double volume = VECTOR(OD->Elementos[o].volumes)[d];
            if (volume < TAPAS_TOLERANCE) continue;

            double sp_cost = VECTOR(WS->SPT[o].dist)[destino];
            total_min_cost += sp_cost * volume;
        }
    }

    double gap_aprox = 1.0;
    if (total_travel_time > TAPAS_TOLERANCE && total_min_cost > 0.0) {
        gap_aprox = 1.0 - (total_min_cost / total_travel_time);
    }

    printf("Warm Start GAP aprox: %.6e\n", gap_aprox);
    if (gap_aprox < TAPAS_EPSILON && gap_aprox >= 0.0) {
        WS->has_warm_start = true;
        return;
    }
  }
  // end_time = clock();
  // printf("Tempo gasto na atribuição inicial: %e segundos\n",
  // (double)(end_time - start_time) / CLOCKS_PER_SEC);
  double previous_gap = DBL_MAX, tempo;
  /* Fase 2: Iterações de equilíbrio */
  // printf("\n=== Iterações de Equilíbrio ===\n");
  // double BF_initial = Beckman_function(solucao, &WS->grafo_bush,
  // BPR_PARAMETERS); double BF = BF_initial; printf("Custo inicial (Beckmann):
  // %f\n", BF_initial);
  int iteracao;
  double gap = 0.0;

  for (iteracao = 0; iteracao < ITAPAS_MAX_ITER; iteracao++) {
    /* Processa cada origem */
    // start_time = clock();
    for (int o = 0; o < OD->size; o++) {
      processar_origem(BPR_PARAMETERS, OD, WS, o, solucao, iteracao, gap);
    }
    // end_time = clock();
    // printf("Tempos : %.2e ", (double)(end_time - start_time) /
    // CLOCKS_PER_SEC); printf("Número total de PAS após processamento das
    // origens: %d\n", num_pas);
    /* Deslocamento global de fluxo */
    if (WS->num_pas == 0) {
      printf("Nenhum PAS identificado. Encerrando iterações.\n");
      break;
    }
    // start_time = clock();
    deslocamento_global_pas(BPR_PARAMETERS, OD, WS, solucao);
    // end_time = clock();
    // printf("%.2e ", (double)(end_time - start_time) / CLOCKS_PER_SEC);
    // BF = Beckman_function(solucao, &WS->grafo_bush, BPR_PARAMETERS);
    // printf("%.2e ", BF/BF_initial);
    /* Calcula e verifica critério de convergência */
    // start_time = clock();
    gap = relative_gap(solucao, &WS->grafo_bush, BPR_PARAMETERS, OD, NULL, NULL);
    // end_time = clock();
    // tempo = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    // printf("%.2e\n", tempo);
    printf("Iteração %3d | GAP: %.6e\n", iteracao + 1, gap); 
    if (gap < TAPAS_EPSILON) {
      // printf("\nConvergência alcançada!\n");
      break;
    }
    // for (int i = 0; i < BPR_PARAMETERS->L; i++) printf("Arco (%ld,%ld):
    // %f\n", IGRAPH_FROM(&WS->grafo_bush, i)+1, IGRAPH_TO(&WS->grafo_bush,
    // i)+1, VECTOR(*solucao)[i]);
    previous_gap = gap;
  }
  // printf("Iteração %3d | GAP: %.6e | PAS ativos: %d\n", iteracao , gap,
  // num_pas);

  /* Salva custos BPR finais em WS->time para uso no GAP aproximado do próximo warm start */
  if (igraph_vector_size(&WS->time) != BPR_PARAMETERS->L) {
    igraph_vector_destroy(&WS->time);
    igraph_vector_init(&WS->time, BPR_PARAMETERS->L);
  }
  BPR(&WS->time, BPR_PARAMETERS, solucao);

  /* Libera memória do conjunto de PAS */
  WS->has_warm_start = true;
}