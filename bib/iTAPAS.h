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
    igraph_vector_int_t canon_unico; // Segmento menor (lexicograficamente) e ordenado
};


/**
 * @brief Calcula e armazena a forma canônica de um PAS.
 * Modifica o PAS in-place, preenchendo canon_s1 e canon_s2.
 */
void pas_canonizar(struct PAS* pas) {
    // 1. Inicializa o vetor canônico com o conteúdo de c1.
    igraph_vector_int_init_copy(&pas->canon_unico, &pas->c1);
    igraph_vector_int_append(&pas->canon_unico, &pas->c2);
    
    igraph_vector_int_sort(&pas->canon_unico);
    //igraph_vector_int_insert(&pas->canon_unico, L, -1);
}


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
    if (pas->canon_unico.stor_begin != NULL) {
        igraph_vector_int_destroy(&pas->canon_unico);
    }
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
        struct PAS *tmp = (struct PAS*) realloc(*vec, novo_tamanho * sizeof(struct PAS));
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
    sprintf(origem_attr, "demanda_%d", pas->origem);
    
    double tempo_c1, tempo_c2, derivada_c1, derivada_c2;
    double capacidade_c1, capacidade_c2;
    
    calcular_custos_caminho(Grafo, BPR_PARAMETERS, &pas->c1, solucao,
                            origem_attr, &tempo_c1, &derivada_c1, &capacidade_c1);
    
    calcular_custos_caminho(Grafo, BPR_PARAMETERS, &pas->c2, solucao,
                            origem_attr, &tempo_c2, &derivada_c2, &capacidade_c2);
    
    double denominador = derivada_c1 + derivada_c2;

    if(fabs(tempo_c2 - tempo_c1) < TAPAS_EPSILON || fabs(denominador) < TAPAS_EPSILON) {
        return 0.0;
    }
    double delta_fluxo = (tempo_c2 - tempo_c1) / denominador;
    //pas_print(pas, Grafo);
    //printf("%.6f %.6f %.6f %.6f\n", delta_fluxo,denominador, tempo_c2, tempo_c1);
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
    for(int n=0; n < 20; n++){
        //pas_print(&conjunto_pas[n],Grafo);
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
                //pas_free(pas_atual);
                pas_remove_at(&conjunto_pas, num_pas, i);
                i--;
                continue;
            }
            
            double dx = deslocar_fluxo_no_pas(Grafo, BPR_PARAMETERS, pas_atual, *solucao, pas_atual->origem);
            //printf("Deslocamento global no PAS da origem %d: %.6f\n", pas_atual->origem, dx);
        }
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
        if(maior_fluxo < TAPAS_EPSILON){
            //printf("Erro: Nenhum arco com fluxo suficiente encontrado para o nó %d ao construir c2 do PAS.\n", no_atual);
            igraph_vector_int_clear(&pas->c2);
            igraph_vector_int_clear(&pas->c1);
            igraph_vector_int_destroy(&arcos_entrantes);
            free(visitado);
            break; 
        }
        /* Evita adicionar arcos já presentes em c2 (prevenção de loop) */
        if (igraph_vector_int_contains(&pas->c2, melhor_arco)) {
            //printf("Erro: Arco %d já presente em c2; interrompendo construção de c2 para evitar ciclo.\n", melhor_arco);
            igraph_vector_int_clear(&pas->c2);
            igraph_vector_int_clear(&pas->c1);
            igraph_vector_int_destroy(&arcos_entrantes);
            free(visitado);
            break; 
        }
        if (melhor_arco == -1) break;
        
        igraph_vector_int_push_back(&pas->c2, melhor_arco);
        no_atual = IGRAPH_FROM(Grafo, melhor_arco);
        
        if (no_atual == fonte) break;
    }
    /* Remove elementos de c1 que estão após o ponto de encontro */
    if(igraph_vector_int_size(&pas->c2) != 0) {
        
        while (no_atual != fonte) {
            int arco = VECTOR(*antecessores)[no_atual];
            if (arco == -1) break;
            
            igraph_vector_int_pop_back(&pas->c1);
            no_atual = IGRAPH_FROM(Grafo, arco);
        }
        
        igraph_vector_int_destroy(&arcos_entrantes);
        free(visitado);
    }
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
    igraph_get_shortest_paths_dijkstra(Grafo, NULL, NULL, fonte, igraph_vss_all(),&tempo_arcos, IGRAPH_OUT, NULL, &antecessores);
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
        int from = IGRAPH_FROM(Grafo, arco);
        int to = IGRAPH_TO(Grafo, arco);
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
    //printf("Caminhos mínimos calculados para a origem %d.\n", fonte);
    
    /* Processa cada arco candidato à remoção */
    int num_candidatos = igraph_vector_int_size(&arcos_removidos);
    for (int i = 0; i < num_candidatos; i++) {
        //printf("  Processando arco %d/%d para remoção...\n", i + 1, num_candidatos);
        struct PAS novo_pas;
        int arco = VECTOR(arcos_removidos)[i];
        //if(fonte == 34218) printf("Arco candidato à remoção: (%ld, %ld)\n", IGRAPH_FROM(Grafo, arco)+1, IGRAPH_TO(Grafo, arco)+1);
        
        
        identificar_pas_fluxo_maximo(Grafo, arco, fonte, &antecessores, &novo_pas);
        //pas_print(&novo_pas, Grafo);
        if( igraph_vector_int_size(&novo_pas.c1) == 0 && 
            igraph_vector_int_size(&novo_pas.c2) == 0) {
            //if(fonte == 34218) printf("    PAS inválido (duplo caminho vazio). Pulando...\n");
            pas_free(&novo_pas);
            continue;
        }
        double deslocamento = deslocar_fluxo_no_pas(Grafo, BPR_PARAMETERS, &novo_pas, 
                                                     *solucao, fonte);
        //printf("    Deslocamento de fluxo no PAS: %.6f\n", deslocamento);
        if (fabs(deslocamento) < TAPAS_EPSILON) {
            //if(fonte == 34218) printf("    Deslocamento nulo. Pulando...\n");
            pas_free(&novo_pas);
            continue;
        }
        bool pas_duplicado = false;
        pas_canonizar(&novo_pas);
        int k;
        for (int j = 0; j < *num_pas; j++) {
        // A comparação é um teste de tamanho e um único memcmp. Máxima performance.
            if (igraph_vector_int_size(&novo_pas.canon_unico) == igraph_vector_int_size(&(*conjunto_pas)[j].canon_unico) &&
                memcmp(VECTOR(novo_pas.canon_unico), VECTOR((*conjunto_pas)[j].canon_unico), igraph_vector_int_size(&novo_pas.canon_unico) * sizeof(int)) == 0) {
                
                int l1 = igraph_vector_int_size(&novo_pas.c1);
                int l2 = igraph_vector_int_size(&(*conjunto_pas)[j].c1);
                if(l1 != l2) continue; // Tamanhos diferentes, não são iguais.
                /* if(fonte == 34218){
                    igraph_vector_int_print(&novo_pas.canon_unico);
                    igraph_vector_int_print(&(*conjunto_pas)[j].canon_unico);
                } */
                pas_duplicado = true;
                break;
            }
        }

        if (!pas_duplicado) {
            //printf("(%ld, %ld) - %d - %.7f - %.7f - %.7f\n", IGRAPH_FROM(Grafo, arco)+1, IGRAPH_TO(Grafo, arco)+1,fonte+1   , deslocamento,VECTOR(custo_nos)[IGRAPH_FROM(Grafo, arco)], VECTOR(custo_nos)[IGRAPH_TO(Grafo, arco)]);
            /* O PAS é único. Adiciona uma cópia profunda ao conjunto. */
            
            // 1. Aumenta o tamanho do array de estruturas PAS.
            struct PAS *tmp = (struct PAS*) realloc(*conjunto_pas, (*num_pas + 1) * sizeof(struct PAS));
            if (tmp == NULL) {
                perror("Falha ao alocar memória para novo PAS");
                pas_free(&novo_pas); // Limpa a memória do PAS que não pôde ser adicionado.
                exit(EXIT_FAILURE);
            }
            *conjunto_pas = tmp;

            // 2. Obtém um ponteiro para a nova posição no array.
            struct PAS *destino_pas = &(*conjunto_pas)[*num_pas];

            // 3. Realiza a "cópia na força bruta", inicializando cada vetor no destino
            //    e copiando o conteúdo do PAS de origem.
            igraph_vector_int_init_copy(&destino_pas->c1, &novo_pas.c1);
            igraph_vector_int_init_copy(&destino_pas->c2, &novo_pas.c2);
            igraph_vector_int_init_copy(&destino_pas->canon_unico, &novo_pas.canon_unico);
            destino_pas->origem = novo_pas.origem;

            // 4. Incrementa o contador de PAS.
            (*num_pas)++;
            
            // 5. Agora que uma cópia profunda e independente foi criada,
            //    libera a memória do PAS temporário `novo_pas`.
            

        }
        pas_free(&novo_pas);
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

    int valor = count_files_in_dir("./output/iTAPAS");
    char filename[10000];
    sprintf(filename, "./output/iTAPAS/gap_progression_%d.txt", valor + 1);
    FILE* gap_file = fopen(filename, "w");

    for (int i = 0; i < OD->size; i++) {
        atribuicao_inicial(BPR_PARAMETERS, OD, &grafo_bush, i, solucao);
    }
    clock_t start_time, end_time;
    double previous_gap = DBL_MAX ,tempo;
    /* Fase 2: Iterações de equilíbrio */
    //printf("\n=== Iterações de Equilíbrio ===\n");
    for (int iteracao = 0; iteracao < ITAPAS_MAX_ITER; iteracao++) {
        /* Processa cada origem */
        start_time = clock();
        //printf("\n--- Iteração %d ---\n", BPR_PARAMETERS->L);
        for (int o = 0; o < OD->size; o++) {
            //printf("Processando origem %d\n", OD->Elementos[o].fonte);
            processar_origem(BPR_PARAMETERS, OD, &grafo_bush, o, solucao, &conjunto_pas, &num_pas);
        }
        /* Deslocamento global de fluxo */
        deslocamento_global_pas(BPR_PARAMETERS, OD, &grafo_bush, conjunto_pas, &num_pas, solucao);
        
        /* Calcula e verifica critério de convergência */
        double gap = relative_gap(solucao, &grafo_bush, BPR_PARAMETERS, OD);
        printf("Iteração %3d | GAP: %.6e | PAS ativos: %d\n",        iteracao + 1, gap, num_pas);
        end_time = clock();
        tempo = (double)(end_time - start_time) / CLOCKS_PER_SEC;
        fprintf(gap_file, "%d %e %f %f\n", iteracao + 1, gap,fabs(previous_gap - gap)/previous_gap, tempo);
        if (gap < 1e-10) {
            printf("\nConvergência alcançada!\n");
            break;
        }
        //for (int i = 0; i < BPR_PARAMETERS->L; i++) printf("Arco (%ld,%ld): %f\n", IGRAPH_FROM(Grafo, i)+1, IGRAPH_TO(Grafo, i)+1, VECTOR(*solucao)[i]);
        previous_gap = gap;
    }
    
    /* Libera memória do conjunto de PAS */
    for (int i = 0; i < num_pas; i++) {
        pas_free(&conjunto_pas[i]);
    }
    free(conjunto_pas);
    
    /* Imprime solução final */
    /* printf("\n=== Fluxos Finais ===\n");
    for (int i = 0; i < BPR_PARAMETERS->L; i++) {
        if (VECTOR(*solucao)[i] > TAPAS_EPSILON) {
            printf("Arco %4d: %10.2f\n", i, VECTOR(*solucao)[i]);
        }
    } */
    fclose(gap_file);
    
    igraph_destroy(&grafo_bush);
}