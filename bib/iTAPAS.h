#pragma once
#include "igraph.h"
#include "define.h"
#include "calc.h"
#include "network.h"
#include "BPR.h"
#include "mtwister.h"

int ITAPAS_MAX_ITER = 1000;

struct PAS{
    igraph_vector_int_t c1;
    igraph_vector_int_t c2;
    int origem;
};

/* Imprime o conteúdo de uma PAS */
void pas_print(const struct PAS *pas,igraph_t *Grafo){
    if(!pas) return;
    printf("\n");
    int i, n1 = igraph_vector_int_size(&pas->c1), n2 = igraph_vector_int_size(&pas->c2);
    printf("PAS(fonte=%d)\n", pas->origem);
    printf("  c1 (size=%d):", n1);
    int from,to;
    for(i = n1-1; i >= 0; i--){
        to = IGRAPH_TO(Grafo, VECTOR(pas->c1)[i])+1;
        from = IGRAPH_FROM(Grafo, VECTOR(pas->c1)[i])+1;
        printf(" (%d, %d)", from, to);
    }
    printf("\n");
    printf("  c2 (size=%d):", n2);
    for(i = n2-1; i >= 0; i--){
        to = IGRAPH_TO(Grafo, VECTOR(pas->c2)[i])+1;
        from = IGRAPH_FROM(Grafo, VECTOR(pas->c2)[i])+1;
        printf(" (%d, %d)", from, to);
    }
    printf("\n");
    printf("\n");
}

/* Libera os recursos internos de uma PAS (não libera o ponteiro PAS em si) */
void pas_free(struct PAS *pas){
    if(!pas) return;
    /* destruir vetores internos */
    igraph_vector_int_destroy(&pas->c1);
    igraph_vector_int_destroy(&pas->c2);
    /* opcional: limpar campos */
    pas->origem = 0;
}

/* Remove o elemento no índice idx de um vetor dinâmico de PAS.
   vec: ponteiro para o array de PAS (será realocado)
   n: ponteiro para o tamanho atual do array (será atualizado)
*/
void pas_remove_at(struct PAS **vec, int *n, int idx){
    if(!vec || !(*vec) || !n) return;
    if(idx < 0 || idx >= *n) return;
    /* destruir o conteúdo do PAS que será removido */
    pas_free(&(*vec)[idx]);

    /* mover os elementos posteriores para preencher o buraco */
    if(idx < (*n - 1)){
        memmove(&(*vec)[idx], &(*vec)[idx + 1], ((*n - 1) - idx) * sizeof(struct PAS));
    }

    /* reduzir o tamanho do array */
    int new_n = *n - 1;
    if(new_n > 0){
        struct PAS *tmp = (struct PAS*) realloc(*vec, new_n * sizeof(struct PAS));
        if(tmp) {
            *vec = tmp;
        } /* se realloc falhar, mantemos o bloco anterior; o tamanho lógico será atualizado abaixo */
    } else {
        /* novo tamanho zero: liberar o bloco e definir ponteiro NULL */
        free(*vec);
        *vec = NULL;
    }
    *n = new_n;
}



double deslocar_fluxo_no_pas(igraph_t* Grafo, struct PARAMETERS* BPR_PARAMETERS, struct PAS* pas, igraph_vector_t solucao,int fonte){
    char origin[20];
    sprintf(origin, "demanda_%d", fonte);
    int L = igraph_vector_int_size(&pas->c1),alvo,edge_id,i ;
    
    double capacidade, fluxo_atual, delta_fluxo = 0,denominador = 0.0,time_1 = 0,time_2 = 0,dtime_1 = 0,dtime_2 = 0;
    double mu1 = 1e10, mu2 = 1e10;

    for(i = 0; i < L; i++){
        edge_id = VECTOR(pas->c1)[i];
        time_1 += single_BPR(VECTOR(solucao)[edge_id], VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
        dtime_1 += single_BPR_derivate(VECTOR(solucao)[edge_id], VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
        mu1 = fmin(mu1,EAN(Grafo,origin,edge_id));
    }
    L = igraph_vector_int_size(&pas->c2);
    for(i = 0; i < L; i++){
        edge_id = VECTOR(pas->c2)[i];
        time_2 += single_BPR(VECTOR(solucao)[edge_id], VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
        dtime_2 += single_BPR_derivate(VECTOR(solucao)[edge_id], VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
        mu2 = fmin(mu2,EAN(Grafo,origin,edge_id));
    }

    denominador = dtime_1 + dtime_2;
    delta_fluxo = (time_2 - time_1) / denominador;
    //printf("Delta fluxo calculado: %f %f %f\n",delta_fluxo,mu2,mu1);

    if(fabs(delta_fluxo) > 1e-9){
        if(delta_fluxo > 0){
            delta_fluxo = fmin(delta_fluxo, mu2);
            double cost;
            L = igraph_vector_int_size(&pas->c1);
            for(i = 0; i < L; i++){
                edge_id = VECTOR(pas->c1)[i];
                if(edge_id == -1) break;
                VECTOR(solucao)[edge_id] += delta_fluxo;
                cost = EAN(Grafo,origin,edge_id);
                igraph_cattribute_EAN_set(Grafo,origin,edge_id,cost + delta_fluxo);
            }
            L = igraph_vector_int_size(&pas->c2);
            for(i = 0; i < L; i++){
                edge_id = VECTOR(pas->c2)[i];
                if(edge_id == -1) break;
                VECTOR(solucao)[edge_id] -= delta_fluxo;
                cost = EAN(Grafo,origin,edge_id);
                igraph_cattribute_EAN_set(Grafo,origin,edge_id,cost - delta_fluxo);
            }
            return delta_fluxo;

        } else {
            delta_fluxo = -delta_fluxo;
            delta_fluxo = fmin(delta_fluxo, mu1);
            double cost;
            L = igraph_vector_int_size(&pas->c1);
            for(i = 0; i < L; i++){
                edge_id = VECTOR(pas->c1)[i];
                if(edge_id == -1) break;
                VECTOR(solucao)[edge_id] -= delta_fluxo;
                cost = EAN(Grafo,origin,edge_id);
                igraph_cattribute_EAN_set(Grafo,origin,edge_id,cost - delta_fluxo);
            }
            L = igraph_vector_int_size(&pas->c2);
            for(i = 0; i < L; i++){
                edge_id = VECTOR(pas->c2)[i];
                if(edge_id == -1) break;
                VECTOR(solucao)[edge_id] += delta_fluxo;
                cost = EAN(Grafo,origin,edge_id);
                igraph_cattribute_EAN_set(Grafo,origin,edge_id,cost + delta_fluxo);
            }
            return delta_fluxo;
        }
    }
    else return 0.0;
}

void deslocamento_global_pas(struct PARAMETERS* BPR_PARAMETERS,struct OD_MATRIX *OD,igraph_t* Grafo,struct PAS* conjunto_pas,int *n_pas, igraph_vector_t *solucao){
    int i,L,j;
    double dx,mu;
    if(*n_pas > 0){
        for(i = 0; i < *n_pas; i++){
            struct PAS* pas = &conjunto_pas[i];
            L = igraph_vector_int_size(&pas->c2);
            mu = 1e10;
            for(j = 0; j < L; j++){
                int edge_id = VECTOR(pas->c2)[j];
                if(edge_id == -1) break;
                char origin[20];
                sprintf(origin, "demanda_%d", pas->origem);
                mu = fmin(mu,EAN(Grafo,origin,edge_id));
            }
            if(mu < 1e-9){
                pas_free(pas);
                pas_remove_at(&conjunto_pas, n_pas, i);
                i--; // Ajustar índice após remoção
                continue;
            }
            else{
                deslocar_fluxo_no_pas(Grafo, BPR_PARAMETERS, pas, *solucao, pas->origem);
            }
        }

    }
}

void identificar_pas_fluxo_maximo(igraph_t* Grafo,int edge_id,int fonte,igraph_vector_int_t *inbounds,struct PAS* pas){
    char origin[20];
    sprintf(origin, "demanda_%d", fonte);
    int cabeca = IGRAPH_TO(Grafo, edge_id),id,i,id_big;
    double cost = 0.0;
    int N = igraph_vector_int_size(inbounds);
    bool* visited = (bool*) calloc(N, sizeof(bool));
    igraph_vector_int_t in_edges;
    igraph_vector_int_init(&in_edges, 0);

    pas->origem = fonte;
    igraph_vector_int_init(&pas->c1, 0);
    igraph_vector_int_init(&pas->c2, 0);
    while(cabeca != fonte){
        id = VECTOR(*inbounds)[cabeca];
        if(id == -1) break;
        // Processar a aresta id conforme necessário
        cabeca = IGRAPH_FROM(Grafo, id);
        
        igraph_vector_int_push_back(&pas->c1, id);
        visited[cabeca] = true;
    }
    cabeca = IGRAPH_FROM(Grafo, edge_id);
    id_big = -1;
    igraph_vector_int_push_back(&pas->c2, edge_id);
    while(!visited[cabeca]){
        igraph_incident(Grafo, &in_edges, cabeca, IGRAPH_IN,IGRAPH_NO_LOOPS);
        int n_in_edges = igraph_vector_int_size(&in_edges);
        for(int i = 0; i < n_in_edges; i++){
            id = VECTOR(in_edges)[i];
            //printf("%d %f\n",id,EAN(Grafo,origin,id));
            if(EAN(Grafo,origin,id) > cost){
                cost = EAN(Grafo,origin,id);
                id_big = id;
            }
            
        }
        cost = 0.0;
        if(id_big == -1) break;
        igraph_vector_int_push_back(&pas->c2, id_big);
        cabeca = IGRAPH_FROM(Grafo, id_big);
        if(cabeca == fonte) break;
    }
    while(cabeca != fonte){
        id = VECTOR(*inbounds)[cabeca];
        if(id == -1) break;
        // Processar a aresta id conforme necessário
        igraph_vector_int_pop_back(&pas->c1);
        cabeca = IGRAPH_FROM(Grafo, id);
    }
    igraph_vector_int_destroy(&in_edges);
    free(visited);
    //exit(0);
}

void processar_origem(struct PARAMETERS* BPR_PARAMETERS,struct OD_MATRIX *OD,igraph_t* Grafo,int o,igraph_vector_t *solucao, struct PAS* conjunto_pas,int *n_pas){
    int fonte = OD->Elementos[o].fonte,i,from,to;
    double custo;
    char origin[20];
    sprintf(origin, "demanda_%d", fonte);

    igraph_vector_t tempo,cost;
    igraph_vector_int_t removed_edges,max_path,inbounds;

    igraph_vector_int_init(&removed_edges, 0);
    igraph_vector_int_init(&max_path, 0);
    igraph_vector_int_init(&inbounds, 0);
    
    igraph_vector_init(&tempo, BPR_PARAMETERS->L);
    igraph_vector_init(&cost, BPR_PARAMETERS->N);
    BPR(&tempo, BPR_PARAMETERS, solucao); // Atualiza o tempo de cada aresta com base no fluxo atual
    igraph_get_shortest_paths_dijkstra(
        Grafo, NULL,NULL,fonte,igraph_vss_all(),&tempo, IGRAPH_OUT,NULL,&inbounds
    );
    for(i= 0;i<BPR_PARAMETERS->N;i++){
        int alvo = i;
        while(alvo != fonte){
            int edge_id = VECTOR(inbounds)[alvo];
            if(edge_id == -1) break;
            VECTOR(cost)[i] += VECTOR(tempo)[edge_id];
            alvo = IGRAPH_FROM(Grafo, edge_id);
        }
    }
    for(i= 0;i<BPR_PARAMETERS->L;i++){
        from = IGRAPH_FROM(Grafo, i);
        to = IGRAPH_TO(Grafo, i);
        custo = EAN(Grafo,origin,i);
        if(custo > 1e-9){
            custo = VECTOR(tempo)[i] + VECTOR(cost)[from] - VECTOR(cost)[to];
            if(custo > 1e-9){
                //printf("Removed edge: %d from %d to %d with reduced cost %f\n",i,from,to,custo);
                igraph_vector_int_push_back(&removed_edges, i);
            }
        }
        // Armazena ou processa a distância conforme necessário
    }

    int L = igraph_vector_int_size(&removed_edges);
    double dx;
    for(i = 0;i<L;i++){
        struct PAS new_pas;
        int edge_id = VECTOR(removed_edges)[i];
        //printf("Processing removed edge %d from %ld to %ld\n",edge_id,IGRAPH_FROM(Grafo, edge_id)+1,IGRAPH_TO(Grafo, edge_id)+1);
        from = IGRAPH_FROM(Grafo, edge_id);
        to = IGRAPH_TO(Grafo, edge_id);
        identificar_pas_fluxo_maximo(Grafo, edge_id, fonte, &inbounds,&new_pas);
        //pas_print(&new_pas,Grafo);
        dx = deslocar_fluxo_no_pas(Grafo, BPR_PARAMETERS, &new_pas, *solucao, fonte);
        printf("dx = %f\n",dx);
        if(fabs(dx) < 1e-9){
            pas_free(&new_pas);
            continue;
        }
        else{
            conjunto_pas = (struct PAS*) realloc(conjunto_pas,  (*n_pas + 1) * sizeof(struct PAS));
            if(conjunto_pas == NULL) exit(0);
            conjunto_pas[*n_pas] = new_pas;
            (*n_pas)++;
        }
        //exit(0);
        //double delta_fluxo = deslocar_fluxo_no_pas(Grafo, BPR_PARAMETERS, &new_pas.c1, &inbounds, *solucao, fonte);
    }

    //exit(0);
    igraph_vector_destroy(&tempo);
    igraph_vector_int_destroy(&removed_edges);

}

void atribuicao_inicial(struct PARAMETERS* BPR_PARAMETERS,struct OD_MATRIX *OD,igraph_t* Grafo,int o,igraph_vector_t *solucao){

    int i,j,fonte = OD->Elementos[o].fonte;
    igraph_vector_int_t inbound;
    igraph_vector_t demanda;
    igraph_vector_int_init(&inbound, 0);
    igraph_vector_init(&demanda, BPR_PARAMETERS->L);
    igraph_get_shortest_paths_dijkstra(
        Grafo, NULL,NULL,fonte,igraph_vss_all(),&BPR_PARAMETERS->cost_time, IGRAPH_OUT,NULL,&inbound
    );
    int n= igraph_vector_int_size(&OD->Elementos[o].alvos);
    for(j = 0; j < n; j++){
        int alvo = VECTOR(OD->Elementos[o].alvos)[j];
        int antecessor = alvo;
        while(antecessor != fonte){
            int index = VECTOR(inbound)[antecessor];
            if(index == -1) exit(0); // Se não encontrar a aresta, encerra o programa
            VECTOR(demanda)[index] += VECTOR(OD->Elementos[o].volumes)[j];
            VECTOR(*solucao)[index] += VECTOR(OD->Elementos[o].volumes)[j];
            antecessor = IGRAPH_FROM(Grafo, index);
        }
    }
    char origin[20];
    sprintf(origin, "demanda_%d", fonte);
    igraph_cattribute_EAN_setv(Grafo, origin, &demanda);
    igraph_vector_int_destroy(&inbound);
    igraph_vector_destroy(&demanda);    
}

void iTAPAS(struct PARAMETERS* BPR_PARAMETERS,struct OD_MATRIX *OD,igraph_t* Grafo,igraph_vector_t *solucao){
    igraph_t bush_graph;
    igraph_vector_t tempo;
    igraph_vector_init(&tempo, BPR_PARAMETERS->L);
    igraph_vector_init(solucao,BPR_PARAMETERS->L);
    igraph_copy(&bush_graph, Grafo);
    int i, fonte, n_pas = 0;

    struct PAS* conjunto_pas = (struct PAS*) malloc(0 * sizeof(struct PAS));

    for (i = 0; i < OD->size; i++){
        atribuicao_inicial(BPR_PARAMETERS,OD,&bush_graph,i,solucao);
        printf("Initial assignment for origin %d completed.\n", OD->Elementos[i].fonte);
    }



    for(i = 0; i < ITAPAS_MAX_ITER; i++){
        for(int o = 0; o < OD->size; o++){
            printf(" ===== Processing origin %d at iteration %d ===== \n", OD->Elementos[o].fonte, i+1);
            processar_origem(BPR_PARAMETERS,OD,&bush_graph,o,solucao,conjunto_pas,&n_pas);
        }
        exit(0);
        deslocamento_global_pas(BPR_PARAMETERS, OD, &bush_graph, conjunto_pas, &n_pas, solucao);

    }
}