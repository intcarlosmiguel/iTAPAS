#pragma once

#include <network.h>
#include <define.h>
#include "igraph.h"

void BPR(igraph_vector_t* tempo,struct PARAMETERS* BPR_PARAMETERS,igraph_vector_t* fluxo){
    double s;
    for (int i = 0; i < BPR_PARAMETERS->L; i++){
        s = pow(VECTOR(*fluxo)[i]/VECTOR(BPR_PARAMETERS->capacidade)[i],BETA);
        VECTOR(*tempo)[i] = (VECTOR(BPR_PARAMETERS->cost_time)[i])*(1 + ALPHA*s);
    }
}

double single_BPR(double fluxo,double free_flow_time,double capacidade){
    double s = pow(fluxo/capacidade,BETA);
    if(isnan(s)) {
        printf("Warning: NaN detected in single_BPR\n");
        exit(0);
        s = 0.0;
    }
    return free_flow_time * (1 + ALPHA * s);
}

double single_BPR_derivate(double fluxo,double free_flow_time,double capacidade){
    double s = pow(fluxo/capacidade,BETA-1);
    if(isnan(s)) {
        printf("Warning: NaN detected in single_BPR_derivate\n");
        exit(0);
        s = 0.0;
    }
    return free_flow_time * ALPHA * BETA * s / capacidade;
}

void BPR_derivate(igraph_vector_t* derivate,struct PARAMETERS* BPR_PARAMETERS,igraph_vector_t* fluxo){
    double s;
    for (int i = 0; i < BPR_PARAMETERS->L; i++){
        s = pow(VECTOR(*fluxo)[i]/VECTOR(BPR_PARAMETERS->capacidade)[i],BETA-1);
        if(isnan(s)) {
            printf("Warning: NaN detected in BPR_derivate\n");
            exit(0);
            s = 0.0;
        }
        VECTOR(*derivate)[i] = (VECTOR(BPR_PARAMETERS->cost_time)[i])*ALPHA*BETA*s/ VECTOR(BPR_PARAMETERS->capacidade)[i];
    }
}

void BPR_gradient(igraph_vector_t* gradiente,struct PARAMETERS* BPR_PARAMETERS,igraph_vector_t* fluxo,double *total_time){
    double s;
    if (total_time != NULL) *total_time = 0;
    for (int i = 0; i < BPR_PARAMETERS->L; i++){
        s = pow(VECTOR(*fluxo)[i]/VECTOR(BPR_PARAMETERS->capacidade)[i],BETA);
        if(isnan(s)) {
            printf("Warning: NaN detected in BPR_derivate\n");
            exit(0);
            s = 0.0;
        }
        if (total_time != NULL) *total_time += (VECTOR(BPR_PARAMETERS->cost_time)[i] * VECTOR(*fluxo)[i]) * (1.0 + 0.03 * s);
        VECTOR(*gradiente)[i] = (VECTOR(BPR_PARAMETERS->cost_time)[i])*(1 + ALPHA*(1+BETA)*s);
    }
}

double frank_wolfe(struct PARAMETERS* BPR_PARAMETERS,igraph_vector_t* fluxo,igraph_vector_t* tempo,igraph_vector_t* direcao,const double* objetivo_incial,igraph_vector_t* y,double* stp){
    double direcao_gradiente = 0, passo = 1.0;
    int i,iteracoes = 0;
    for ( i = 0; i < BPR_PARAMETERS->L; i++) direcao_gradiente += VECTOR(*direcao)[i]*VECTOR(*tempo)[i];
    if(direcao_gradiente > 0){
        printf("Gradiente não ótimo!\n");
        exit(0);
    }
    
    double objetivo = 0;
    double dgtest = ftol * direcao_gradiente;
    
    
    while (true){
        for ( i = 0; i < BPR_PARAMETERS->L; i++) VECTOR(*y)[i] = VECTOR(*fluxo)[i] + VECTOR(*direcao)[i]*passo;
        BPR_gradient(tempo,BPR_PARAMETERS,y,&objetivo);
        iteracoes++;
        //printf("%d %f %f\n",iteracoes,objetivo,*objetivo_incial + passo * dgtest);
        if (objetivo > *objetivo_incial + passo * dgtest) passo *= decresse;
        else break;

        if(passo < min_step) break;
        if(passo > max_step) break;
		if(iteracoes > MAXIMO_ITERACOES) break;
    }
    *stp = passo;
    return objetivo;
    
}

double relative_gap(
    igraph_vector_t *flow,
    igraph_t *Grafo,
    struct PARAMETERS* BPR_PARAMETERS,
    struct OD_MATRIX* OD
) {

    double gap = 0.0;
    double total_cost = 0.0;
    double total_flow = 0.0;
    int i,j,alvo,edge_id,fonte,size;

    igraph_vector_t time;
    igraph_vector_init(&time, BPR_PARAMETERS->L);
    BPR(&time, BPR_PARAMETERS, flow); // Calcula o tempo de cada aresta
    for ( i = 0; i < BPR_PARAMETERS->L; i++) total_flow += VECTOR(*flow)[i]*VECTOR(time)[i];

    
    
    for( i = 0; i<OD->size; i++) {
        fonte = OD->Elementos[i].fonte;
        size = igraph_vector_int_size(&OD->Elementos[i].alvos);

        igraph_vector_int_t inbound;
        igraph_vector_int_init(&inbound, 0);
        igraph_get_shortest_paths_dijkstra(
            Grafo, NULL,NULL,fonte,igraph_vss_all(),&time, IGRAPH_OUT,NULL,&inbound
        );

        for (j = 0; j < size; j++) {
            alvo = VECTOR(OD->Elementos[i].alvos)[j];
            while(alvo != fonte) {
                edge_id = VECTOR(inbound)[alvo];
                total_cost += VECTOR(time)[edge_id] * VECTOR(OD->Elementos[i].volumes)[j];
                alvo = IGRAPH_FROM(Grafo, edge_id);
            }
        }

        igraph_vector_int_destroy(&inbound);
    }
    if (total_flow > 0) gap = 1 - total_cost / total_flow;
    else gap = 1.0; // Evita divisão por zero, assume que o gap é 1 se não houver fluxo
    igraph_vector_destroy(&time);


    return gap;
}

void leblanc(struct PARAMETERS* BPR_PARAMETERS,int** edge_list,struct OD_MATRIX *OD,igraph_t* Grafo,igraph_vector_t *solucao){

    igraph_vector_t tempo;
    igraph_vector_init(&tempo, BPR_PARAMETERS->L);

    int i,iteracoes = 0;
    
    igraph_vector_init(solucao,BPR_PARAMETERS->L);
    igraph_vector_t gradiente;
    igraph_vector_init(&gradiente,BPR_PARAMETERS->L);
    double objetivo = 0,objetivo2 = 0,dx,df,dy;
    atualiza_fluxo(Grafo,OD,edge_list,solucao,&BPR_PARAMETERS->cost_time);
    //else parallel_atualiza_fluxo(Grafo,OD,edge_list,solucao,&BPR_PARAMETERS->cost_time);
    igraph_vector_t novo_fluxo;
    igraph_vector_init(&novo_fluxo,BPR_PARAMETERS->L);

    double stp = 0;
    
    while(true){


        BPR_gradient(&tempo,BPR_PARAMETERS,solucao,&objetivo);
        atualiza_fluxo(Grafo,OD,edge_list,&gradiente,&tempo);
        //else parallel_atualiza_fluxo(Grafo,OD,edge_list,solucao,&BPR_PARAMETERS->cost_time);
        
        for ( i = 0; i < BPR_PARAMETERS->L; i++) VECTOR(gradiente)[i] -= VECTOR(*solucao)[i];

        objetivo2 = frank_wolfe(BPR_PARAMETERS,solucao,&tempo,&gradiente,&objetivo,&novo_fluxo,&stp);

        dx = 0.0;
		for( i = 0; i < BPR_PARAMETERS->L; i++){
			dy = fabs(VECTOR(*solucao)[i]-VECTOR(novo_fluxo)[i]) / (fabs(VECTOR(*solucao)[i]) + 1);
			if(dx < dy) dx = dy;
            VECTOR(*solucao)[i] = VECTOR(novo_fluxo)[i];
		}
        df = (objetivo - objetivo2) / objetivo;
        
        iteracoes++;
        if((iteracoes)%100 == 0)printf("%d %e %e\n",iteracoes,dx,df);
        if(iteracoes > MAXIMO_ITERACOES) break;
        if(dx < X_TOLERANCIA) break;
        if(df < 1E-4) break;
    }
    FILE *file;
    file = fopen("./output/solution.dat","w");
    for ( i = 0; i < BPR_PARAMETERS->L; i++) fprintf(file,"%f ",VECTOR(*solucao)[i]);
    fclose(file);
    igraph_vector_destroy(&novo_fluxo);
    igraph_vector_destroy(&tempo);
    igraph_vector_destroy(&gradiente);
}
