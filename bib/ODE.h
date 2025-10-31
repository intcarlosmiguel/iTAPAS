#pragma once
#include "igraph.h"
#include "network.h"
#include "define.h"
#include "dial.h"
#include "simulate.h"

void create_random_OD_matrix(struct OD_MATRIX* OD, int N, int min_volume, int max_volume) {
    OD->size = N;
    OD->Elementos = (struct ElementOD*) malloc(N * sizeof(struct ElementOD));
    
    for(int i = 0; i < N; i++) {
        igraph_vector_int_init(&OD->Elementos[i].alvos, 0);
        igraph_vector_init(&OD->Elementos[i].volumes, 0);
        
        // For each origin, create random number of destinations
        int num_destinations = 1 + rand() % (N / 10); // At least 1 destination, up to N/10
        for(int j = 0; j < num_destinations; j++) {
            int target;
            do {
                target = rand() % N;
            } while(target == i); // Avoid self-loops
            
            igraph_vector_int_push_back(&OD->Elementos[i].alvos, target);
            int volume = min_volume + (rand() % (max_volume - min_volume + 1));
            igraph_vector_push_back(&OD->Elementos[i].volumes, volume);
        }
    }
}

void remove_small_volumes(struct OD_MATRIX* OD) {
    for (int i = 0; i < OD->size; i++) {
        int L = igraph_vector_int_size(&OD->Elementos[i].alvos);
        int j = 0;
        while (j < L) {
            if (VECTOR(OD->Elementos[i].volumes)[j] < 1) {
                // Remove element at index j
                igraph_vector_remove(&OD->Elementos[i].volumes, j);
                igraph_vector_int_remove(&OD->Elementos[i].alvos, j);
                L--;
            } else {
                j++;
            }
        }
    }
}


double F2(
    igraph_vector_t *flow,
    igraph_vector_t *flow_hat,
    int L
) {
    // Compute the objective function value
    double objective_value = 0.0;
    for (int i = 0; i <L; i++) {
        objective_value += pow(VECTOR(*flow_hat)[i] - VECTOR(*flow)[i], 2)/2;
    }
    //printf("Objective Function Value: %.20f\n", objective_value);
    return objective_value;
}

double F1(
    struct OD_MATRIX* OD
) {
    // Compute the objective function value
    double objective_value = 0.0, volume = 0.0;
    int L;
    for (int i = 0; i < OD->size; i++) {
        L = igraph_vector_int_size(&OD->Elementos[i].alvos);
        for (int j = 0; j < L; j++) {
            volume = VECTOR(OD->Elementos[i].volumes)[j];
            if (volume <= 0) continue; // Skip if volume is zero or negative
            objective_value += volume * log10(volume);
        }
    }
    //printf("Objective Function Value: %.20f\n", objective_value);
    return objective_value;
}

double ObjectiveFunction(
    igraph_vector_t *time,
    igraph_vector_t *time_hat,
    struct OD_MATRIX* OD,
    int L,
    double lambda
) {
    // Compute the objective function value
    double objective_value = 0.0;
    objective_value += F1(OD);
    objective_value += lambda*F2(time, time_hat, L);
    //printf("Objective Function Value: %.20f\n", objective_value);
    return objective_value;
}

double calc_probability(
    struct BUSH* bush,
    int alvo,
    struct PARAMETERS* BPR_PARAMETERS,
    igraph_vector_t *flow,
    igraph_vector_t *flow_hat
){
    double* prob = (double*) calloc(BPR_PARAMETERS->L, sizeof(double));
    igraph_vector_int_t incident_edges;
    igraph_vector_int_init(&incident_edges, 0);
    igraph_vector_int_t nodes;
    igraph_vector_int_init(&nodes, 0);
    igraph_vector_int_push_back(&nodes, alvo);
    bool * visited_nodes = (bool*) calloc(BPR_PARAMETERS->N, sizeof(bool));
    visited_nodes[alvo] = true;
    int n = 1,node,i,j,edge,from;
    for(i = 0; i < n; i++) {
        node = VECTOR(nodes)[i];
        igraph_incident(&bush->Grafo, &incident_edges, node, IGRAPH_IN,IGRAPH_NO_LOOPS);
        for(j = 0; j < igraph_vector_int_size(&incident_edges); j++) {
            edge = VECTOR(incident_edges)[j];
            from = IGRAPH_FROM(&bush->Grafo, edge);
            if(!visited_nodes[from]) {
                igraph_vector_int_push_back(&nodes, from);
                visited_nodes[from] = true;
                n++;
            }
        }
    }
    double* kappa = (double*) calloc(BPR_PARAMETERS->N, sizeof(double));
    kappa[alvo] = 1.0; // Set kappa for the target node
    double total_flow = 0.0;
    for(i = 0; i < n; i++) {
        total_flow = 0.0;
        node = VECTOR(nodes)[i];
        igraph_incident(&bush->Grafo, &incident_edges, node, IGRAPH_IN, IGRAPH_NO_LOOPS);
        int tamanho = igraph_vector_int_size(&incident_edges);
        igraph_vector_t phi;
        igraph_vector_init(&phi, tamanho);
        for(j = 0; j < tamanho; j++) {
            edge = VECTOR(incident_edges)[j];
            total_flow += VECTOR(bush->flow_per_origin)[VECTOR(bush->edge_id)[edge]];
        }
        for(j = 0; j < tamanho; j++) {
            edge = VECTOR(incident_edges)[j];
            VECTOR(phi)[j] = VECTOR(bush->flow_per_origin)[VECTOR(bush->edge_id)[edge]] / total_flow;
            if(VECTOR(phi)[j] > 0){
                prob[edge] = kappa[node] * VECTOR(phi)[j];
                from = IGRAPH_FROM(&bush->Grafo, edge);
                kappa[from] += prob[edge]; // Update kappa for the source node
            }
        }
        igraph_vector_destroy(&phi);
    }
    double d = 0;
    for(i = 0; i < BPR_PARAMETERS->L; i++) {
        d += prob[i]*(VECTOR(*flow)[i] - VECTOR(*flow_hat)[i]);
    }
    for(i = 0; i < BPR_PARAMETERS->L; i++) {
        if(prob[i] > 0)printf("(%ld,%ld) - %f\n", IGRAPH_FROM(&bush->Grafo, i) + 1, IGRAPH_TO(&bush->Grafo, i) + 1, prob[i]);
    }
    free(prob);
    free(visited_nodes);
    free(kappa);
    igraph_vector_int_destroy(&incident_edges);
    igraph_vector_int_destroy(&nodes);
    printf("d: %.20f\n", d);
    return d;
    

}


void OD_gradient(
    struct OD_MATRIX* OD,
    struct BUSH* bush,
    struct PARAMETERS* BPR_PARAMETERS,
    igraph_vector_t *flow,
    igraph_vector_t *flow_hat,
    double lambda,
    double eta
){

    int size = OD->size;
    double norma = 0.0;
    for(int i = 0; i < size; i++) {
        int L = igraph_vector_int_size(&OD->Elementos[i].alvos);
        if(L == 0) continue;
        for(int j = 0; j < L; j++) {
            int alvo = VECTOR(OD->Elementos[i].alvos)[j];
            double P = calc_probability(&bush[i], alvo, BPR_PARAMETERS, flow, flow_hat);
            double volume = VECTOR(OD->Elementos[i].volumes)[j];
            if(volume <= 0) continue; // Skip if volume is zero or negative
            VECTOR(OD->Elementos[i].volumes)[j] -= eta*( log10(volume) + lambda * P);
            norma += pow(eta*( log10(volume) + lambda * P), 2);
            if (VECTOR(OD->Elementos[i].volumes)[j] < 1) {
                // Remove element at index j
                igraph_vector_remove(&OD->Elementos[i].volumes, j);
                igraph_vector_int_remove(&OD->Elementos[i].alvos, j);
                L--;
            }
            else j++;

        }
    }
    norma = sqrt(norma);
    double obj = ObjectiveFunction(flow, flow_hat, OD, BPR_PARAMETERS->L, lambda);
    double RMSE = F2(flow, flow_hat, BPR_PARAMETERS->L);
    printf("Updated Objective Function Value: %.20f\n", obj);
    printf("Norm of the gradient: %.20f\n", norma);
    printf("RMSE: %.20f\n", RMSE);
    printf("OD volumes updated.\n");
}

void read_flow_hat(const char* filename, igraph_vector_t* flow_hat) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    igraph_vector_init(flow_hat, 0);
    double value;
    while (fscanf(file, "%lf", &value) == 1) {
        igraph_vector_push_back(flow_hat, value);
    }

    fclose(file);
}

void OD_ESTIMATION(
    double lambda,
    double eta,
    const char* file_edges,
    const char* file_flow_hat,
    int iter
){
    struct PARAMETERS BPR_PARAMETERS;
    struct OD_MATRIX OD_MATRIX;
    igraph_vector_t flow_hat;
    igraph_vector_t flow;
    igraph_vector_int_t edges;

    read_flow_hat(file_flow_hat, &flow_hat);
    igraph_vector_int_init(&edges, 0);
    init_parameters(&BPR_PARAMETERS,&edges,file_edges);

    //load_OD_from_file(file_od, &OD_MATRIX);
    print_OD_matrix(&OD_MATRIX);

    igraph_t Grafo;
    igraph_empty(&Grafo, BPR_PARAMETERS.N, IGRAPH_DIRECTED);
    igraph_add_edges(&Grafo, &edges, NULL);
    struct BUSH* bushes;
    for(int i = 0; i < iter; i++) {
        Dial(&Grafo, &OD_MATRIX, &BPR_PARAMETERS, &flow, &bushes,false);
        OD_gradient(&OD_MATRIX, bushes, &BPR_PARAMETERS, &flow, &flow_hat, lambda, eta);
    }
        

    igraph_destroy(&Grafo);
}