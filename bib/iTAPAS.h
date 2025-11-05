#pragma once
#include "igraph.h"
#include "define.h"
#include "calc.h"
#include "mtwister.h"

int ITAPAS_MAX_ITER = 100;
double TAPAS_EPSILON = 1e-9;

/* Estrutura que representa um Par de Arcos Alternativos (PAS) */
struct PAS {
    igraph_vector_int_t c1;      /* Caminho 1 (arcos) */
    igraph_vector_int_t c2;      /* Caminho 2 (arcos) */
    int origem;                   /* Nó de origem do PAS */
};

/* ========================================================================
   FUNÇÕES AUXILIARES PARA MANIPULAÇÃO DE PAS
   ======================================================================== */

void pas_print(const struct PAS *pas, igraph_t *Grafo) {
    if (!pas) return;
    
    printf("\nPAS(fonte=%d)\n", pas->origem);
    
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
    if (!pas) return;
    igraph_vector_int_destroy(&pas->c1);
    igraph_vector_int_destroy(&pas->c2);
    pas->origem = 0;
}

void pas_remove_at(struct PAS **vec, int *tamanho, int indice) {
    if (!vec || !(*vec) || !tamanho) return;
    if (indice < 0 || indice >= *tamanho) return;
    
    pas_free(&(*vec)[indice]);
    
    if (indice < (*tamanho - 1)) {
        memmove(&(*vec)[indice], &(*vec)[indice + 1], 
                ((*tamanho - 1) - indice) * sizeof(struct PAS));
    }
    
    int novo_tamanho = *tamanho - 1;
    if (novo_tamanho > 0) {
        struct PAS *tmp = (struct PAS*) realloc(*vec, novo_tamanho * sizeof(struct PAS));
        if (tmp) *vec = tmp;
    } else {
        free(*vec);
        *vec = NULL;
    }
    *tamanho = novo_tamanho;
}

/* ========================================================================
   CÁLCULO DE CUSTOS E DERIVADAS
   ======================================================================== */

static void calcular_custos_caminho(igraph_t *Grafo, 
                                     struct PARAMETERS *BPR_PARAMETERS,
                                     igraph_vector_int_t *caminho,
                                     igraph_vector_t solucao,
                                     const char *origem_attr,
                                     double *tempo_total,
                                     double *derivada_total,
                                     double *capacidade_minima) {
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
        *capacidade_minima = fmin(*capacidade_minima, EAN(Grafo, origem_attr, arco_id));
    }
}

/* ========================================================================
   DESLOCAMENTO DE FLUXO
   ======================================================================== */

double deslocar_fluxo_no_pas(igraph_t *Grafo, 
                              struct PARAMETERS *BPR_PARAMETERS,
                              struct PAS *pas, 
                              igraph_vector_t solucao,
                              int fonte) {
    char origem_attr[20];
    sprintf(origem_attr, "demanda_%d", fonte);
    
    double tempo_c1, tempo_c2, derivada_c1, derivada_c2;
    double capacidade_c1, capacidade_c2;
    
    calcular_custos_caminho(Grafo, BPR_PARAMETERS, &pas->c1, solucao,
                            origem_attr, &tempo_c1, &derivada_c1, &capacidade_c1);
    
    calcular_custos_caminho(Grafo, BPR_PARAMETERS, &pas->c2, solucao,
                            origem_attr, &tempo_c2, &derivada_c2, &capacidade_c2);
    
    double denominador = derivada_c1 + derivada_c2;
    double delta_fluxo = (tempo_c2 - tempo_c1) / denominador;
    
    if (fabs(delta_fluxo) <= TAPAS_EPSILON) {
        return 0.0;
    }
    
    /* Determina direção e magnitude do deslocamento */
    igraph_vector_int_t *caminho_origem, *caminho_destino;
    double capacidade_limite;
    
    if (delta_fluxo > 0) {
        delta_fluxo = fmin(delta_fluxo, capacidade_c2);
        caminho_origem = &pas->c1;
        caminho_destino = &pas->c2;
    } else {
        delta_fluxo = fmin(-delta_fluxo, capacidade_c1);
        caminho_origem = &pas->c2;
        caminho_destino = &pas->c1;
    }
    
    /* Atualiza fluxos no caminho de origem (adiciona fluxo) */
    for (int i = 0; i < igraph_vector_int_size(caminho_origem); i++) {
        int arco = VECTOR(*caminho_origem)[i];
        if (arco == -1) break;
        
        VECTOR(solucao)[arco] += delta_fluxo;
        double custo_atual = EAN(Grafo, origem_attr, arco);
        igraph_cattribute_EAN_set(Grafo, origem_attr, arco, custo_atual + delta_fluxo);
    }
    
    /* Atualiza fluxos no caminho de destino (remove fluxo) */
    for (int i = 0; i < igraph_vector_int_size(caminho_destino); i++) {
        int arco = VECTOR(*caminho_destino)[i];
        if (arco == -1) break;
        
        VECTOR(solucao)[arco] -= delta_fluxo;
        double custo_atual = EAN(Grafo, origem_attr, arco);
        igraph_cattribute_EAN_set(Grafo, origem_attr, arco, custo_atual - delta_fluxo);
    }
    
    return delta_fluxo;
}

void deslocamento_global_pas(struct PARAMETERS *BPR_PARAMETERS,
                              struct OD_MATRIX *OD,
                              igraph_t *Grafo,
                              struct PAS *conjunto_pas,
                              int *num_pas,
                              igraph_vector_t *solucao) {
    if (*num_pas <= 0) return;
    
    for (int i = 0; i < *num_pas; i++) {
        struct PAS *pas_atual = &conjunto_pas[i];
        
        /* Verifica capacidade mínima do caminho c2 */
        char origem_attr[20];
        sprintf(origem_attr, "demanda_%d", pas_atual->origem);
        
        double capacidade_min = 1e10;
        int num_arcos_c2 = igraph_vector_int_size(&pas_atual->c2);
        
        for (int j = 0; j < num_arcos_c2; j++) {
            int arco = VECTOR(pas_atual->c2)[j];
            if (arco == -1) break;
            capacidade_min = fmin(capacidade_min, EAN(Grafo, origem_attr, arco));
        }
        
        /* Remove PAS se capacidade é muito pequena */
        if (capacidade_min < TAPAS_EPSILON) {
            pas_free(pas_atual);
            pas_remove_at(&conjunto_pas, num_pas, i);
            i--;
            continue;
        }
        
        deslocar_fluxo_no_pas(Grafo, BPR_PARAMETERS, pas_atual, *solucao, pas_atual->origem);
    }
}

/* ========================================================================
   IDENTIFICAÇÃO DE PAS
   ======================================================================== */

void identificar_pas_fluxo_maximo(igraph_t *Grafo,
                                   int arco_inicial,
                                   int fonte,
                                   igraph_vector_int_t *antecessores,
                                   struct PAS *pas) {
    char origem_attr[20];
    sprintf(origem_attr, "demanda_%d", fonte);
    
    pas->origem = fonte;
    igraph_vector_int_init(&pas->c1, 0);
    igraph_vector_int_init(&pas->c2, 0);
    
    int num_nos = igraph_vector_int_size(antecessores);
    bool *visitado = (bool*) calloc(num_nos, sizeof(bool));
    
    /* Constrói caminho c1: do nó cabeça até a fonte usando antecessores */
    int no_atual = IGRAPH_TO(Grafo, arco_inicial);
    while (no_atual != fonte) {
        int arco = VECTOR(*antecessores)[no_atual];
        if (arco == -1) break;
        
        igraph_vector_int_push_back(&pas->c1, arco);
        visitado[no_atual] = true;
        no_atual = IGRAPH_FROM(Grafo, arco);
    }
    visitado[no_atual] = true;
    
    /* Constrói caminho c2: do nó cauda até encontrar c1, usando arco de maior fluxo */
    no_atual = IGRAPH_FROM(Grafo, arco_inicial);
    igraph_vector_int_push_back(&pas->c2, arco_inicial);
    
    igraph_vector_int_t arcos_entrantes;
    igraph_vector_int_init(&arcos_entrantes, 0);
    
    while (!visitado[no_atual]) {
        igraph_incident(Grafo, &arcos_entrantes, no_atual, IGRAPH_IN, IGRAPH_NO_LOOPS);
        
        int melhor_arco = -1;
        double maior_fluxo = 0.0;
        
        int num_entrantes = igraph_vector_int_size(&arcos_entrantes);
        for (int i = 0; i < num_entrantes; i++) {
            int arco = VECTOR(arcos_entrantes)[i];
            double fluxo = EAN(Grafo, origem_attr, arco);
            
            if (fluxo > maior_fluxo) {
                maior_fluxo = fluxo;
                melhor_arco = arco;
            }
        }
        
        if (melhor_arco == -1) break;
        
        igraph_vector_int_push_back(&pas->c2, melhor_arco);
        no_atual = IGRAPH_FROM(Grafo, melhor_arco);
        
        if (no_atual == fonte) break;
    }
    
    /* Remove elementos de c1 que estão após o ponto de encontro */
    while (no_atual != fonte) {
        int arco = VECTOR(*antecessores)[no_atual];
        if (arco == -1) break;
        
        igraph_vector_int_pop_back(&pas->c1);
        no_atual = IGRAPH_FROM(Grafo, arco);
    }
    
    igraph_vector_int_destroy(&arcos_entrantes);
    free(visitado);
}

/* ========================================================================
   PROCESSAMENTO POR ORIGEM
   ======================================================================== */

void processar_origem(struct PARAMETERS *BPR_PARAMETERS,
                      struct OD_MATRIX *OD,
                      igraph_t *Grafo,
                      int indice_origem,
                      igraph_vector_t *solucao,
                      struct PAS **conjunto_pas,
                      int *num_pas) {
    int fonte = OD->Elementos[indice_origem].fonte;
    char origem_attr[20];
    sprintf(origem_attr, "demanda_%d", fonte);
    
    igraph_vector_t tempo_arcos, custo_nos;
    igraph_vector_int_t arcos_removidos, antecessores;
    
    igraph_vector_init(&tempo_arcos, BPR_PARAMETERS->L);
    igraph_vector_init(&custo_nos, BPR_PARAMETERS->N);
    igraph_vector_int_init(&arcos_removidos, 0);
    igraph_vector_int_init(&antecessores, 0);
    
    /* Calcula tempos de viagem e caminho mais curto */
    BPR(&tempo_arcos, BPR_PARAMETERS, solucao);
    igraph_get_shortest_paths_dijkstra(Grafo, NULL, NULL, fonte, igraph_vss_all(),
                                       &tempo_arcos, IGRAPH_OUT, NULL, &antecessores);
    
    /* Calcula custos mínimos de cada nó até a fonte */
    for (int no = 0; no < BPR_PARAMETERS->N; no++) {
        VECTOR(custo_nos)[no] = 0.0;
        int no_atual = no;
        
        while (no_atual != fonte) {
            int arco = VECTOR(antecessores)[no_atual];
            if (arco == -1) break;
            
            VECTOR(custo_nos)[no] += VECTOR(tempo_arcos)[arco];
            no_atual = IGRAPH_FROM(Grafo, arco);
        }
    }
    
    /* Identifica arcos com custo reduzido positivo e fluxo positivo */
    for (int arco = 0; arco < BPR_PARAMETERS->L; arco++) {
        double fluxo = EAN(Grafo, origem_attr, arco);
        if (fluxo <= TAPAS_EPSILON) continue;
        
        int no_origem = IGRAPH_FROM(Grafo, arco);
        int no_destino = IGRAPH_TO(Grafo, arco);
        
        double custo_reduzido = VECTOR(tempo_arcos)[arco] + 
                                VECTOR(custo_nos)[no_origem] - 
                                VECTOR(custo_nos)[no_destino];
        
        if (custo_reduzido > TAPAS_EPSILON) {
            igraph_vector_int_push_back(&arcos_removidos, arco);
        }
    }
    
    /* Processa cada arco candidato à remoção */
    int num_candidatos = igraph_vector_int_size(&arcos_removidos);
    for (int i = 0; i < num_candidatos; i++) {
        struct PAS novo_pas;
        int arco = VECTOR(arcos_removidos)[i];
        
        identificar_pas_fluxo_maximo(Grafo, arco, fonte, &antecessores, &novo_pas);
        double deslocamento = deslocar_fluxo_no_pas(Grafo, BPR_PARAMETERS, &novo_pas, 
                                                     *solucao, fonte);
        
        if (fabs(deslocamento) < TAPAS_EPSILON) {
            pas_free(&novo_pas);
            continue;
        }
        
        /* Adiciona novo PAS ao conjunto */
        struct PAS *tmp = (struct PAS*) realloc(*conjunto_pas, 
                                                 (*num_pas + 1) * sizeof(struct PAS));
        if (tmp == NULL) {
            pas_free(&novo_pas);
            perror("Falha ao alocar memória para novo PAS");
            exit(EXIT_FAILURE);
        }
        
        *conjunto_pas = tmp;
        (*conjunto_pas)[*num_pas] = novo_pas;
        (*num_pas)++;
    }
    
    igraph_vector_destroy(&tempo_arcos);
    igraph_vector_destroy(&custo_nos);
    igraph_vector_int_destroy(&arcos_removidos);
    igraph_vector_int_destroy(&antecessores);
}

/* ========================================================================
   ATRIBUIÇÃO INICIAL (ALL-OR-NOTHING)
   ======================================================================== */

void atribuicao_inicial(struct PARAMETERS *BPR_PARAMETERS,
                        struct OD_MATRIX *OD,
                        igraph_t *Grafo,
                        int indice_origem,
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
    int num_destinos = igraph_vector_int_size(&OD->Elementos[indice_origem].alvos);
    for (int i = 0; i < num_destinos; i++) {
        int destino = VECTOR(OD->Elementos[indice_origem].alvos)[i];
        double volume = VECTOR(OD->Elementos[indice_origem].volumes)[i];
        
        int no_atual = destino;
        while (no_atual != fonte) {
            int arco = VECTOR(antecessores)[no_atual];
            if (arco == -1) {
                fprintf(stderr, "Erro: caminho não encontrado de %d para %d\n", fonte, destino);
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

void iTAPAS(struct PARAMETERS *BPR_PARAMETERS,
            struct OD_MATRIX *OD,
            igraph_t *Grafo,
            igraph_vector_t *solucao) {
    igraph_t grafo_bush;
    igraph_copy(&grafo_bush, Grafo);
    
    igraph_vector_init(solucao, BPR_PARAMETERS->L);
    
    struct PAS *conjunto_pas = NULL;
    int num_pas = 0;
    
    /* Fase 1: Atribuição inicial (all-or-nothing) */
    printf("=== Atribuição Inicial ===\n");
    for (int i = 0; i < OD->size; i++) {
        atribuicao_inicial(BPR_PARAMETERS, OD, &grafo_bush, i, solucao);
    }
    
    /* Fase 2: Iterações de equilíbrio */
    printf("\n=== Iterações de Equilíbrio ===\n");
    for (int iteracao = 0; iteracao < ITAPAS_MAX_ITER; iteracao++) {
        /* Processa cada origem */
        for (int o = 0; o < OD->size; o++) {
            processar_origem(BPR_PARAMETERS, OD, &grafo_bush, o, solucao, 
                           &conjunto_pas, &num_pas);
        }
        
        /* Deslocamento global de fluxo */
        deslocamento_global_pas(BPR_PARAMETERS, OD, &grafo_bush, 
                               conjunto_pas, &num_pas, solucao);
        
        /* Calcula e verifica critério de convergência */
        double gap = relative_gap(solucao, &grafo_bush, BPR_PARAMETERS, OD);
        printf("Iteração %3d | GAP: %.6e | PAS ativos: %d\n", 
               iteracao + 1, gap, num_pas);
        
        if (gap < 1e-10) {
            printf("\nConvergência alcançada!\n");
            break;
        }
    }
    
    /* Libera memória do conjunto de PAS */
    for (int i = 0; i < num_pas; i++) {
        pas_free(&conjunto_pas[i]);
    }
    free(conjunto_pas);
    
    /* Imprime solução final */
    printf("\n=== Fluxos Finais ===\n");
    for (int i = 0; i < BPR_PARAMETERS->L; i++) {
        if (VECTOR(*solucao)[i] > TAPAS_EPSILON) {
            printf("Arco %4d: %10.2f\n", i, VECTOR(*solucao)[i]);
        }
    }
    
    igraph_destroy(&grafo_bush);
}