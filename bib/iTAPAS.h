#pragma once
#include "igraph.h"
#include "define.h"
#include "calc.h"
#include "network.h"
#include "BPR.h"
#include "mtwister.h"

int ITAPAS_MAX_ITER = 1000;

void deslocar_fluxo_no_pas(igraph_t* Grafo, struct PARAMETERS* BPR_PARAMETERS, igraph_vector_int_t* PAS,igraph_vector_int_t *inbounds,igraph_vector_t solucao,int fonte){
    char origin[20];
    sprintf(origin, "demanda_%d", fonte);
    int L = igraph_vector_int_size(PAS),alvo,edge_id,i ;
    double capacidade, fluxo_atual, delta_fluxo = 0,denominador = 0.0,time_1 = 0,time_2 = 0,dtime_1 = 0,dtime_2 = 0;
    double mu1 = 1e10, mu2 = 1e10;
    for(i = 0; i < L; i++){
        edge_id = VECTOR(*PAS)[i];
        time_1 += single_BPR(VECTOR(solucao)[edge_id], VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
        dtime_1 += single_BPR_derivate(VECTOR(solucao)[edge_id], VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
        mu1 = fmin(mu1,EAN(Grafo,origin,edge_id));
    }
    alvo = IGRAPH_FROM(Grafo, VECTOR(*PAS)[0]);
    while(alvo != fonte){
        edge_id = VECTOR(*inbounds)[alvo];
        if(edge_id == -1) break;
        time_2 += single_BPR(VECTOR(solucao)[edge_id], VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
        dtime_2 += single_BPR_derivate(VECTOR(solucao)[edge_id], VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
        alvo = IGRAPH_FROM(Grafo, edge_id);
        mu2 = fmin(mu2,EAN(Grafo,origin,edge_id));
    }

    denominador = dtime_1 + dtime_2;
    delta_fluxo = (time_2 - time_1) / denominador;
    if(fabs(delta_fluxo) > 1e-9){
        if(delta_fluxo > 0){
            delta_fluxo = fmin(delta_fluxo, mu1);
        } else {
            delta_fluxo = -delta_fluxo;
            delta_fluxo = fmin(delta_fluxo, mu2);
        }
    }
}

void identificar_pas_fluxo_maximo(igraph_t* Grafo,int edge_id,int fonte,igraph_vector_int_t *inbounds,igraph_vector_int_t* PAS ){
    char origin[20];
    sprintf(origin, "demanda_%d", fonte);
    int cabeca = IGRAPH_TO(Grafo, edge_id),id,i,id_big, L = igraph_vector_int_size(inbounds);
    double cost = 0.0;
    bool* visited = (bool*) calloc(L, sizeof(bool));
    igraph_vector_int_t in_edges;
    igraph_vector_int_init(&in_edges, 0);

    while(cabeca != fonte){
        id = VECTOR(*inbounds)[cabeca];
        if(id == -1) break;
        // Processar a aresta id conforme necessário
        cabeca = IGRAPH_FROM(Grafo, id);
        visited[id] = true;
    }
    cabeca = IGRAPH_FROM(Grafo, edge_id);
    while(true){
        igraph_incident(Grafo, &in_edges, cabeca, IGRAPH_IN,IGRAPH_NO_LOOPS);
        int n_in_edges = igraph_vector_int_size(&in_edges);
        for(int i = 0; i < n_in_edges; i++){
            id = VECTOR(in_edges)[i];
            if(EAN(Grafo,origin,i) > cost){
                cost = EAN(Grafo,origin,i);
                id_big = id;
            }
            
        }
        cost = 0.0;
        if(visited[id_big]) break;
        igraph_vector_int_push_back(PAS, id_big);
        if(cabeca == fonte) break;
        cabeca = IGRAPH_FROM(Grafo, id_big);
    }
    igraph_vector_int_destroy(&in_edges);
    free(visited);
    exit(0);
}

void processar_origem(struct PARAMETERS* BPR_PARAMETERS,struct OD_MATRIX *OD,igraph_t* Grafo,int o,igraph_vector_t *solucao){
    int fonte = OD->Elementos[o].fonte,i,from,to;
    double custo;
    char origin[20];
    sprintf(origin, "demanda_%d", fonte);

    igraph_vector_t tempo,cost;
    igraph_vector_int_t removed_edges,PAS,inbounds;

    igraph_vector_int_init(&removed_edges, 0);
    igraph_vector_int_init(&PAS, 0);
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
                printf("Removed edge: %d from %d to %d with reduced cost %f\n",i,from,to,custo);
                igraph_vector_int_push_back(&removed_edges, i);
            }
        }
        // Armazena ou processa a distância conforme necessário
    }

    int L = igraph_vector_int_size(&removed_edges);

    for(i = 0;i<L;i++){
        int edge_id = VECTOR(removed_edges)[i];
        from = IGRAPH_FROM(Grafo, edge_id);
        to = IGRAPH_TO(Grafo, edge_id);
        igraph_vector_int_clear(&PAS);
        identificar_pas_fluxo_maximo(Grafo, edge_id, fonte, &inbounds, &PAS);
    }

    exit(0);
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
    printf("Initial assignment for origin %d done.\n",fonte);
}

void iTAPAS(struct PARAMETERS* BPR_PARAMETERS,struct OD_MATRIX *OD,igraph_t* Grafo,igraph_vector_t *solucao){
    igraph_t bush_graph;
    igraph_vector_t tempo;
    igraph_vector_init(&tempo, BPR_PARAMETERS->L);
    igraph_vector_init(solucao,BPR_PARAMETERS->L);
    igraph_copy(&bush_graph, Grafo);
    int i, fonte;
    for (i = 0; i < OD->size; i++){

        atribuicao_inicial(BPR_PARAMETERS,OD,&bush_graph,i,solucao);
    }



    for(i = 0; i < ITAPAS_MAX_ITER; i++){
        for(int o = 0; o < OD->size; o++){
            processar_origem(BPR_PARAMETERS,OD,&bush_graph,o,solucao);
        }


    }
}