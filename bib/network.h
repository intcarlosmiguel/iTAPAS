#pragma once

#include <calc.h>
#include "igraph.h"

void print_vector_igraph(igraph_vector_int_t* vetor){
    int N =igraph_vector_int_size(vetor);
    for (int i = 0; i < N; i++){
        if(i != N - 1) printf("%ld,",VECTOR(*vetor)[i]);
        else printf("%ld\n",VECTOR(*vetor)[i]);
    }
    
}

int find_id(int fonte,int alvo,int** edge_list){
    int i = 0;
    while(true){
        if((edge_list[i][0] == fonte)&&(edge_list[i][1] == alvo)) break;
        i++;
    }
    return i;
}

void Dijkstra(igraph_t* Grafo,int fonte,igraph_vector_int_t* alvos,igraph_vector_t* pesos,igraph_vector_int_t* parents){

    igraph_vector_int_list_t vecs, evecs;
    igraph_vector_int_t inbound;
    igraph_vector_int_list_init(&evecs, 0);
    
    igraph_vector_int_init(&inbound, 0);
    igraph_get_shortest_paths_dijkstra(Grafo, NULL,NULL,fonte,igraph_vss_all(),pesos, IGRAPH_OUT,parents,&inbound);
    
    // Print each vector in the lis
    //igraph_get_shortest_paths_bellman_ford(Grafo, &vecs,&evecs,*fonte,alvos_vs,pesos, IGRAPH_IN,parents,&inbound);
    //igraph_vector_int_print(parents);
    igraph_vector_int_destroy(&inbound);
}

void parallel_atualiza_fluxo(igraph_t *Grafo,struct OD_MATRIX* OD,int** edge_list,igraph_vector_t* fluxo, igraph_vector_t *pesos){
    double** fluxos_paralelo = (double**) malloc(THREADS * sizeof(double*));

    int i,j,fonte,c = 0,thread_id;
    for (i = 0; i < THREADS; i++){
        fluxos_paralelo[i] = (double*) calloc( igraph_vector_size(fluxo), sizeof(double));
    }
    igraph_vector_fill(fluxo,0);
    #pragma omp parallel for schedule(dynamic)
    for ( i = 0; i < OD->size; i++){
        thread_id = omp_get_thread_num();
        fonte = OD->Elementos[i].fonte;
        igraph_vector_int_t parents;
        igraph_vector_int_init(&parents, 0);
        Dijkstra(Grafo,fonte,&OD->Elementos[i].alvos,pesos,&parents);

        //printf("Fonte: %d\n",fonte);
        //print_vector_igraph(&parents);
        for (j = 0; j < igraph_vector_int_size(&OD->Elementos[i].alvos); j++) {
            int alvo = VECTOR(OD->Elementos[i].alvos)[j];
            double volume = VECTOR(OD->Elementos[i].volumes)[j];
            if(VECTOR(parents)[alvo] < 0) continue;
            while (alvo != fonte) {
                int antecessor = VECTOR(parents)[alvo];
                int index = find_id(antecessor, alvo, edge_list);
                fluxos_paralelo[thread_id][index] += volume;
                //VECTOR(*fluxo)[index] += volume;
                alvo = antecessor;
            }
            if(VECTOR(OD->Elementos[i].alvos)[j] != fonte) c++;
        }
        igraph_vector_int_destroy(&parents);
    }
    for (i = 0; i < THREADS; i++){
        for (j = 0; j < igraph_vector_size(fluxo); j++){
            VECTOR(*fluxo)[j] += fluxos_paralelo[i][j];
        }
        free(fluxos_paralelo[i]);
    }
    

}
void atualiza_fluxo(igraph_t *Grafo,struct OD_MATRIX* OD,int** edge_list,igraph_vector_t* fluxo, igraph_vector_t *pesos){
    int i,j,fonte,c = 0;
    igraph_vector_fill(fluxo,0);
    for ( i = 0; i < OD->size; i++){
        fonte = OD->Elementos[i].fonte;
        igraph_vector_int_t parents;
        igraph_vector_int_init(&parents, 0);
        Dijkstra(Grafo,fonte,&OD->Elementos[i].alvos,pesos,&parents);

        //printf("Fonte: %d\n",fonte);
        //print_vector_igraph(&parents);
        for (j = 0; j < igraph_vector_int_size(&OD->Elementos[i].alvos); j++) {
            int alvo = VECTOR(OD->Elementos[i].alvos)[j];
            double volume = VECTOR(OD->Elementos[i].volumes)[j];
            if(VECTOR(parents)[alvo] < 0) continue;
            while (alvo != fonte) {
                int antecessor = VECTOR(parents)[alvo];
                int index = find_id(antecessor, alvo, edge_list);
                VECTOR(*fluxo)[index] += volume;
                alvo = antecessor;
            }
            if(VECTOR(OD->Elementos[i].alvos)[j] != fonte) c++;
        }
        igraph_vector_int_destroy(&parents);
    }
}
