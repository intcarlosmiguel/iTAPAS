#pragma once
#include "igraph.h"
#include "define.h"
#include "calc.h"
#include "network.h"
#include "example.h"
#include "BPR.h"
#include "mtwister.h"
#include "dial.h"
#include "iTAPAS.h"

#include <sys/stat.h>
#include <errno.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir_dir(path) _mkdir(path)
#else
#define mkdir_dir(path) mkdir(path, 0755)
#endif

void check_or_create_dir(const char *path) {
    struct stat st = {0};

    // Checa se a pasta existe
    if (stat(path, &st) == -1) {
        // Pasta não existe, tenta criar
        if (mkdir_dir(path) == 0) {
            printf("Diretório '%s' criado com sucesso.\n", path);
        }
    }
}

void init_time(double* initial_time, struct PARAMETERS BPR_PARAMETERS, struct OD_MATRIX* OD_MATRIX, igraph_t* Grafo){
    igraph_vector_t time;
    igraph_vector_t free_flow;
    igraph_vector_init(&time, BPR_PARAMETERS.L);
    igraph_vector_init(&free_flow, BPR_PARAMETERS.L);
    igraph_vector_update(&free_flow, &BPR_PARAMETERS.capacidade);
    BPR(&time, &BPR_PARAMETERS, &free_flow);
    int index = 0;
    for(int i = 0; i < OD_MATRIX->size; i++){
        int fonte = OD_MATRIX->Elementos[i].fonte;
        igraph_matrix_t res;
        igraph_matrix_init(&res, 0, 0);
        igraph_distances(Grafo,&time,&res,igraph_vss_1(fonte), igraph_vss_all(), IGRAPH_OUT);
        for(int j = 0; j < igraph_vector_int_size(&OD_MATRIX->Elementos[i].alvos); j++){
            int node = VECTOR(OD_MATRIX->Elementos[i].alvos)[j];
            initial_time[index] = MATRIX(res,0,node);
            index++;
        }
        igraph_matrix_destroy(&res);

    }
    igraph_vector_destroy(&time);
    igraph_vector_destroy(&free_flow);
}

void init_simulate(struct PARAMETERS* BPR_PARAMETERS,struct OD_MATRIX *OD,igraph_vector_t *solucao,igraph_vector_int_t * edges){
    igraph_t Grafo;
    igraph_empty(&Grafo, BPR_PARAMETERS->N, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, edges, NULL);
    char algoritmo[1000];
    snprintf(algoritmo, sizeof(algoritmo), "iTAPAS");
    if (strcmp(algoritmo, "Leblanc") == 0) {
        printf("Using Leblanc's algorithm for shortest paths.\n");
        leblanc(BPR_PARAMETERS, NULL, OD, &Grafo, solucao);
    } else if (strcmp(algoritmo, "Dial") == 0) {
        struct BUSH* bushes;
        printf("Using Dial's algorithm for shortest paths.\n");
        Dial(&Grafo, OD, BPR_PARAMETERS, solucao, &bushes, false);
        for (int i = 0; i < OD->size; i++) free_bush(&bushes[i]);
        free(bushes);
    } else if (strcmp(algoritmo, "iTAPAS") == 0) {
        struct BUSH* bushes;
        printf("Using iTAPAS's algorithm for shortest paths.\n");
        iTAPAS(BPR_PARAMETERS, OD, &Grafo, solucao);
    } else {
        printf("Error: Unknown algorithm '%s'\n", algoritmo);
    }
    //struct BUSH* bushes;
    //Dial(&Grafo,OD,BPR_PARAMETERS,solucao,&bushes,false);

    /* for (int i = 0; i < OD->size; i++) {
        for (int j = 0; j < igraph_vector_int_size(&OD->Elementos[i].alvos); j++) {
            VECTOR(OD->Elementos[i].volumes)[j] += 100.0;
            VECTOR(OD->Elementos[i].warm_volumes)[j] = 100.0;
        }
    }

    Dial(&Grafo,OD,BPR_PARAMETERS,solucao,&bushes,true);

    FILE* file_flow;
    file_flow = fopen("./example/flow.txt", "w");
    for(int i = 0; i < BPR_PARAMETERS->L; i++) fprintf(file_flow, "%f %f\n", VECTOR(*solucao)[i],single_BPR(VECTOR(*solucao)[i], VECTOR(BPR_PARAMETERS->cost_time)[i], VECTOR(BPR_PARAMETERS->capacidade)[i]));
    fclose(file_flow); */
    //igraph_vector_print(solucao);
    //optimize(BPR_PARAMETERS,edge_list,OD,&Grafo,solucao);
    igraph_destroy(&Grafo);
}




void simulate_example(const char* arquivoEDGES,const char* arquivoOD){

    struct PARAMETERS BPR_PARAMETERS;
    struct OD_MATRIX OD_MATRIX;


    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    init_parameters(&BPR_PARAMETERS,&edges,arquivoEDGES);

    load_OD_from_file(arquivoOD, &OD_MATRIX);
    print_OD_matrix(&OD_MATRIX);
    igraph_vector_t solucao;
    init_simulate(&BPR_PARAMETERS,&OD_MATRIX,&solucao,&edges);
}

void simulate_percolation(const char* arquivoEDGES,const char* arquivoOD,int K,int seed){
    struct PARAMETERS BPR_PARAMETERS;
    struct OD_MATRIX OD_MATRIX;
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    init_parameters(&BPR_PARAMETERS,&edges,arquivoEDGES);
    load_OD_from_file(arquivoOD, &OD_MATRIX);
    //print_OD_matrix(&OD_MATRIX);

    igraph_vector_t solucao;
    igraph_t Grafo;
    igraph_empty(&Grafo, BPR_PARAMETERS.N, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    double alpha = 0.05;
    double T = 3.;
    double* initial_time = malloc(OD_MATRIX.n_elements * sizeof(double));
    int index = 0;
    igraph_vector_t result;
    igraph_vector_init(&result, 0);
    init_time(initial_time, BPR_PARAMETERS, &OD_MATRIX, &Grafo);

    int congested = 0;  
    int i,j,node;
    double time;
    FILE *flow_per_alpha;
    int iter = 0;
    char filename[50];
    int soma = 0;
    check_or_create_dir("./results");
    char dir[50];
    sprintf(dir, "./results/%d", K);
    check_or_create_dir(dir);
    int fonte;
    sprintf(filename, "./results/%d/removed_edges_%d.txt", K, iter);
    FILE *file_removed_edges;
    file_removed_edges = fopen(filename, "w");
    double eff = 0;
    while(true){
        
        struct BUSH* bushes;
        Dial(&Grafo,&OD_MATRIX,&BPR_PARAMETERS,&solucao,&bushes,false);

        congested = 0;
        index = 0;
        eff = 0;
        if(seed == 0){
            sprintf(filename, "./results/%d/flow_%d.txt", K, iter);
            flow_per_alpha = fopen(filename, "w");
            for(i = 0; i <BPR_PARAMETERS.L ; i++) fprintf(flow_per_alpha, "%f %f\n", VECTOR(solucao)[i],single_BPR(VECTOR(solucao)[i], VECTOR(BPR_PARAMETERS.cost_time)[i], VECTOR(BPR_PARAMETERS.capacidade)[i]));
        }
        for(i = 0; i < OD_MATRIX.size; i++){
            fonte = OD_MATRIX.Elementos[i].fonte;
            for(j = 0; j < bushes[i].n_alvos; j++){
                node = VECTOR(OD_MATRIX.Elementos[i].alvos)[j];
                time = VECTOR(bushes[i].paths.dist_shortest_local)[node];
                if(time/initial_time[index] > T){
                    congested++;
                    int edge_id,target = node;
                    while(target != OD_MATRIX.Elementos[i].fonte){
                        edge_id = VECTOR(bushes[i].paths.min_edges)[target];
                        VECTOR(BPR_PARAMETERS.cost_time)[edge_id] = 1e10;
                        fprintf(file_removed_edges, "%d %d %d %d\n",iter, edge_id,fonte,node);
                        target = IGRAPH_FROM(&Grafo, edge_id);
                    }
                }
                VECTOR(OD_MATRIX.Elementos[i].volumes)[j] += alpha;
                index++;
            }
        }

        for(i = 0; i < BPR_PARAMETERS.L; i++){
            double total_time = single_BPR(VECTOR(solucao)[i], VECTOR(BPR_PARAMETERS.cost_time)[i], VECTOR(BPR_PARAMETERS.capacidade)[i]);
            eff += VECTOR(solucao)[i]*(VECTOR(BPR_PARAMETERS.cost_time)[i]/total_time);
        }
        eff /= BPR_PARAMETERS.L;
        if(seed == 0) printf("%d %d %d %f\n",iter, congested, OD_MATRIX.n_elements, eff);

        igraph_vector_push_back(&result, (double)congested/OD_MATRIX.n_elements);
        for(i = 0; i < OD_MATRIX.size; i++) free_bush(&bushes[i]);
        if(seed == 0) fclose(flow_per_alpha);
        free(bushes);
        iter++;
        if((double)congested/OD_MATRIX.n_elements > 0.98){
            soma += 1;
            if(soma > 10) break;
        }
        if(congested == OD_MATRIX.n_elements) break;
    }
    sprintf(filename, "./results/%d/percolation_%d.txt", K, seed);
    FILE *file = fopen(filename, "w");

    for (int i = 0; i < igraph_vector_size(&result); i++) {
        printf("%f\n", VECTOR(result)[i]);
        fprintf(file, "%.10e\n", VECTOR(result)[i]);
    }
    fclose(file);
    sprintf(filename, "./results/%d/last_od_%d.txt", K, seed);
    FILE *lastr_od = fopen(filename, "w");
    for (i = 0; i < OD_MATRIX.size; i++) {
        for (j = 0; j < igraph_vector_int_size(&OD_MATRIX.Elementos[i].alvos); j++) {
            fprintf(lastr_od, "%f\n", VECTOR(OD_MATRIX.Elementos[i].volumes)[j]);
        }
    }
    fclose(lastr_od);

    igraph_destroy(&Grafo);
    free(initial_time);
    igraph_vector_destroy(&BPR_PARAMETERS.cost_time);
    igraph_vector_destroy(&BPR_PARAMETERS.capacidade);
    igraph_vector_int_destroy(&edges);
    for (int i = 0; i < OD_MATRIX.size; i++) {
        igraph_vector_int_destroy(&OD_MATRIX.Elementos[i].alvos);
        igraph_vector_destroy(&OD_MATRIX.Elementos[i].volumes);
    }
    free(OD_MATRIX.Elementos);
    igraph_vector_destroy(&solucao);
    igraph_vector_destroy(&result);
}

