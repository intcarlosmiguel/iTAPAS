#include "define.h"
#include "calc.h"

void atualiza_fluxo(igraph_t *Grafo,struct OD_MATRIX* OD,igraph_vector_t* fluxo, igraph_vector_t *pesos){
    int i,j,fonte;
    int edge_id,alvo;
    double volume;
    igraph_vector_fill(fluxo,0);
    for ( i = 0; i < OD->size; i++){

        fonte = OD->Elementos[i].fonte;
        igraph_vector_int_t inbounds;
        igraph_vector_int_init(&inbounds, 0);
        igraph_get_shortest_paths_dijkstra(Grafo, NULL,NULL,fonte,igraph_vss_all(),pesos, IGRAPH_OUT,NULL,&inbounds);

        for (j = 0; j < igraph_vector_int_size(&OD->Elementos[i].alvos); j++) {

            alvo = VECTOR(OD->Elementos[i].alvos)[j];
            volume = VECTOR(OD->Elementos[i].volumes)[j];

            if(VECTOR(inbounds)[alvo] < 0) continue;

            while (alvo != fonte) {
                edge_id = VECTOR(inbounds)[alvo];
                VECTOR(*fluxo)[edge_id] += volume;
                alvo = IGRAPH_FROM(Grafo, edge_id);
            }

        }
        igraph_vector_int_destroy(&inbounds);
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

void leblanc(struct PARAMETERS* BPR_PARAMETERS,struct OD_MATRIX *OD,igraph_t* Grafo,igraph_vector_t *solucao){

    igraph_vector_t tempo;
    igraph_vector_init(&tempo, BPR_PARAMETERS->L);

    int i,iteracoes = 0;
    
    igraph_vector_init(solucao,BPR_PARAMETERS->L);
    igraph_vector_t gradiente;
    igraph_vector_init(&gradiente,BPR_PARAMETERS->L);
    double objetivo = 0,objetivo2 = 0,dx,df,dy;
    atualiza_fluxo(Grafo,OD,solucao,&BPR_PARAMETERS->cost_time);
    //else parallel_atualiza_fluxo(Grafo,OD,edge_list,solucao,&BPR_PARAMETERS->cost_time);
    igraph_vector_t novo_fluxo;
    igraph_vector_init(&novo_fluxo,BPR_PARAMETERS->L);

    double stp = 0;
    
    while(true){


        BPR_gradient(&tempo,BPR_PARAMETERS,solucao,&objetivo);
        atualiza_fluxo(Grafo,OD,&gradiente,&tempo);
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
