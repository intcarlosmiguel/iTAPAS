
#pragma once

#include "igraph.h"
#include "define.h"
#include <time.h>
const int MAX_ITER = 5000;


void att_topological_order(
    struct BUSH *bush,
    igraph_t *Grafo
) {
    igraph_integer_t index;

    igraph_vector_int_destroy(&bush->edge_id);
    igraph_vector_int_init(&bush->edge_id, 0);
    for (long e = 0; e < igraph_ecount(&bush->Grafo); e++) {
        igraph_integer_t from, to;
        igraph_edge(&bush->Grafo, e, &from, &to);
        igraph_get_eid(
            Grafo,&index, from, to, IGRAPH_UNDIRECTED, -1
        );
        if (index == -1) {
            printf("Erro ao encontrar o ID da aresta (%ld,%ld) no grafo original.\n", from, to);
            exit(0); // Encerra o programa se não encontrar a aresta
        }
        igraph_vector_int_push_back(&bush->edge_id, index);
        //printf("Aresta %d (%ld,%ld) - %f %f\n",index, from, to, VECTOR(*flow)[index],bush->flow_per_alvo[index]);
    }
    // Update topological order after adding new edges
    igraph_vector_int_destroy(&bush->topological_order);
    igraph_vector_int_init(&bush->topological_order, 0);
    if (igraph_topological_sorting(&bush->Grafo, &bush->topological_order, IGRAPH_OUT) != IGRAPH_SUCCESS) {
        printf("Error: Graph has cycles after adding new edges.\n");
        exit(1);
    }
    
}

void removeUnusedArcs(
    struct BUSH *bush,
    igraph_t *Grafo,
    int fonte,
    struct PARAMETERS* BPR_PARAMETERS,
    igraph_vector_t *total_flow
) {

    int i, j, edge_id;
    igraph_vector_int_t new_edges; 
    igraph_vector_int_init(&new_edges, 0);

    for (i = 0; i < igraph_vector_int_size(&bush->edge_id); i++) {
        edge_id = VECTOR(bush->edge_id)[i];
        bush->is_ingraph[edge_id] = true;
        if (VECTOR(bush->flow_per_origin)[edge_id] <= 0.0) bush->is_ingraph[edge_id] = false;
        else igraph_vector_int_push_back(&new_edges, edge_id);
    }

    for (i = 0; i < BPR_PARAMETERS->N; i++) {
        if(i == fonte) continue; // Pula a fonte
        edge_id = VECTOR(bush->paths.min_edges)[i];
        if(!bush->is_ingraph[edge_id]) {
            bush->is_ingraph[edge_id] = true;
            igraph_vector_int_push_back(&new_edges, edge_id);
        }
    }
    igraph_es_t edges;
    igraph_es_vector(&edges, &new_edges);
    igraph_destroy(&bush->Grafo); 

    igraph_subgraph_from_edges(
        Grafo, &bush->Grafo, edges, false
    );
    att_topological_order(bush,Grafo); // Atualiza a ordenação topológica da bush
    igraph_vector_int_destroy(&new_edges);
    igraph_es_destroy(&edges);
}


double findFlowDelta(
    struct BUSH *bush,
    igraph_t *Grafo,
    int lca,
    int alvo,
    struct PARAMETERS* BPR_PARAMETERS,
    igraph_vector_t *total_flow
){

    double mu = DBL_MAX;
    int current_node = alvo;
    while(current_node != lca) {
        int edge_id = VECTOR(bush->paths.max_edges)[current_node];
        if (edge_id == -1) { mu = 0; break; }
        mu = fmin(mu, VECTOR(bush->flow_per_origin)[edge_id]);
        current_node = IGRAPH_FROM(Grafo, edge_id);
    }
    if(mu == 0) return 0.0; // Se mu for zero, não há fluxo a ser ajustad
    double delta_x = 0.0,min_path_flow,max_path_flow,min_derivate_flow,max_derivate_flow,new_delta_x = 0.0;
    int edge_id;
    double denominator;
    for(int i = 0; i < MAX_ITER; i++){

        min_path_flow = 0;
        min_derivate_flow = 0;
        max_path_flow = 0;
        max_derivate_flow = 0;
        current_node = alvo;
        while(current_node!=lca){
            edge_id = VECTOR(bush->paths.min_edges)[current_node];
            min_path_flow += single_BPR(VECTOR(*total_flow)[edge_id]+delta_x, VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
            min_derivate_flow += single_BPR_derivate(VECTOR(*total_flow)[edge_id]+delta_x, VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
            current_node = IGRAPH_FROM(Grafo, edge_id);
        }
        current_node = alvo;
        while(current_node!=lca){

            edge_id = VECTOR(bush->paths.max_edges)[current_node];
            max_path_flow += single_BPR(VECTOR(*total_flow)[edge_id]-delta_x, VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
            max_derivate_flow += single_BPR_derivate(VECTOR(*total_flow)[edge_id]-delta_x, VECTOR(BPR_PARAMETERS->cost_time)[edge_id], VECTOR(BPR_PARAMETERS->capacidade)[edge_id]);
            current_node = IGRAPH_FROM(Grafo, edge_id);
        }



        denominator = max_derivate_flow + min_derivate_flow;
        if(denominator == 0) exit(0); // Evita divisão por zero
        new_delta_x = delta_x + (max_path_flow - min_path_flow) / denominator;
        if(fabs(new_delta_x - delta_x) < EPSILON) {
            double ant = delta_x;

            //printf("Convergência alcançada após %d iterações.\n", i+1);
            delta_x = fmin(fmax(new_delta_x,0.0),mu);
            break;
        }
        delta_x = new_delta_x;
    }
    return fmin(fmax(delta_x,0.0),mu); // Garante que delta_x não seja negativo e não exceda mu
}

void init_bush(
    struct BUSH *bush, 
    igraph_t *Grafo, 
    int id_bush,
    struct PARAMETERS* BPR_PARAMETERS,
    struct OD_MATRIX* OD,
    igraph_vector_t *flow,
    bool s
) {
    igraph_vector_int_t inbound;
    
    igraph_vector_int_init(&inbound, 0);
    igraph_vector_int_init(&bush->topological_order,0);

    int fonte = OD->Elementos[id_bush].fonte,j,k,id;
    //printf("Fonte: %d!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", fonte);
    
    //BPR(&time, BPR_PARAMETERS, flow);
    //igraph_vector_print(&BPR_PARAMETERS->cost_time);
    igraph_get_shortest_paths_dijkstra(
        Grafo, NULL,NULL,fonte,igraph_vss_all(),&BPR_PARAMETERS->cost_time, IGRAPH_OUT,NULL,&inbound
    );
    
    int alvo,edge_id,from,to;
    double tempo;
    igraph_vector_int_init(&bush->edge_id, 0);
    bush->is_ingraph = (bool*) calloc(BPR_PARAMETERS->L, sizeof(bool)); // Inicializa o vetor de arestas no grafo
    int target,antecessor;
    igraph_integer_t index;

    igraph_vector_int_t edges_vec;
    igraph_vector_int_init(&edges_vec, 0);

    
    for ( j = 0; j< BPR_PARAMETERS->N; j++){
        if(j == fonte) continue; // Pula a fonte
        alvo = j;
        while(alvo!=fonte){
            edge_id = VECTOR(inbound)[alvo];
            alvo = IGRAPH_FROM(Grafo, edge_id);
            if(bush->is_ingraph[edge_id]) continue; // Se a aresta já está no grafo da bush, pula para a próxima iteração
            igraph_vector_int_push_back(&edges_vec, edge_id);
            bush->is_ingraph[edge_id] = true; // Marca a aresta como parte do grafo da bush
        }
    }
    //igraph_vector_int_print(&edges_vec);
    igraph_empty(&bush->Grafo, BPR_PARAMETERS->N, IGRAPH_DIRECTED);
    igraph_es_t edges;
    igraph_es_vector(&edges,&edges_vec);
    igraph_subgraph_from_edges(
        Grafo, &bush->Grafo, edges, false
    );
    for(j = 0; j < bush->n_alvos; j++) {
        antecessor = VECTOR(OD->Elementos[id_bush].alvos)[j];
        while(antecessor != fonte) {
            index = VECTOR(inbound)[antecessor];
            if(index == -1) exit(0); // Se não encontrar a aresta, encerra o programa
            VECTOR(*flow)[index] += VECTOR(OD->Elementos[id_bush].volumes)[j];
            VECTOR(bush->flow_per_origin)[index] += VECTOR(OD->Elementos[id_bush].volumes)[j];
            antecessor = IGRAPH_FROM(Grafo, index);
        }
    }
    init_path(&bush->paths, BPR_PARAMETERS->N);
    //char filename[100];
    //sprintf(filename, "bush_%d.txt", fonte);
    //FILE* file = fopen(filename, "w");
    for (long e = 0; e < igraph_ecount(&bush->Grafo); e++) {
        igraph_integer_t from, to;
        igraph_edge(&bush->Grafo, e, &from, &to);
        igraph_get_eid(Grafo,&index, from, to, IGRAPH_DIRECTED, -1);
        igraph_vector_int_push_back(&bush->edge_id, index);
    }

    if (igraph_topological_sorting(&bush->Grafo, &bush->topological_order, IGRAPH_OUT) != IGRAPH_SUCCESS) {
        printf("Erro ao ordenar topologicamente o grafo.\n");
        exit(0); // Falha ao ordenar topologicamente, encerra o programa
    }
    
    //fclose(file);
    igraph_vector_int_destroy(&inbound);
    igraph_vector_int_destroy(&edges_vec);
    igraph_es_destroy(&edges);
    
}


void max_Distance(
    igraph_vector_t* longest,
    struct BUSH *bush,
    igraph_t *Grafo,
    igraph_integer_t from,
    struct PARAMETERS *BPR_PARAMETERS,
    igraph_vector_t *total_flow
){
    igraph_vector_init(longest, BPR_PARAMETERS->N);
    igraph_vector_fill(longest, -DBL_MAX); // "Infinito negativo"
    VECTOR(*longest)[from] = 0.0; // Inicializa a distância do nó de origem como 0.0
    
    igraph_vector_int_t incident_edges;
    igraph_vector_int_init(&incident_edges, 0);
    igraph_vector_t time;
    igraph_vector_init(&time, BPR_PARAMETERS->L);
    int i,index,id,edge_id,j;
    double fluxo,weight;
    BPR(&time, BPR_PARAMETERS, total_flow); // Calcula o tempo de cada aresta
    igraph_integer_t eid,u,v_node;

    for (i = 0; i < igraph_vector_int_size(&bush->topological_order); ++i) {
        u = VECTOR(bush->topological_order)[i];
        igraph_incident(&bush->Grafo, &incident_edges, u, IGRAPH_OUT,IGRAPH_NO_LOOPS);
        if (VECTOR(*longest)[u] == -DBL_MAX) continue;
        for (j = 0; j < igraph_vector_int_size(&incident_edges); ++j) {

            eid = VECTOR(incident_edges)[j];
            id = VECTOR(bush->edge_id)[eid]; // Acessa o id da aresta
            v_node = IGRAPH_TO(Grafo, id);
            weight = VECTOR(time)[id]; // Acessa o peso da aresta

            // Relaxamento para o maior caminho
            // Só estender de 'u' se 'u' já tem um caminho mais longo válido (não -DBL_MAX)
            if (VECTOR(*longest)[u] + weight > VECTOR(*longest)[v_node]) {
                VECTOR(*longest)[v_node] = VECTOR(*longest)[u] + weight;
            }
        }
        igraph_vector_int_clear(&incident_edges); // Reutiliza o vetor de arestas
    }
    

    // Libera memória dos vetores locais
    igraph_vector_int_destroy(&incident_edges);
    igraph_vector_destroy(&time);
}

void find_dag_shortest_longest_costs_and_parents(
    struct BUSH *bush,
    igraph_t *Grafo,
    igraph_integer_t from,
    struct PARAMETERS *BPR_PARAMETERS,
    igraph_vector_t *total_flow
) {
    int N = igraph_vcount(Grafo),j;
    igraph_vector_fill(&bush->paths.dist_shortest_local, DBL_MAX);
    igraph_vector_fill(&bush->paths.dist_longest_local, -DBL_MAX); // "Infinito negativo"
    VECTOR(bush->paths.dist_shortest_local)[from] = 0.0;
    VECTOR(bush->paths.dist_longest_local)[from] = 0.0;
    igraph_vector_int_fill(&bush->paths.max_edges, -1);
    igraph_vector_int_fill(&bush->paths.min_edges, -1); // "Infinito negativo"
        
    
    igraph_vector_int_t incident_edges;
    igraph_vector_int_init(&incident_edges, 0);
    igraph_vector_t time;
    igraph_vector_init(&time, BPR_PARAMETERS->L);

    int i,index,id,edge_id;
    double fluxo,weight;
    BPR(&time, BPR_PARAMETERS, total_flow); // Calcula o tempo de cada aresta
    igraph_integer_t eid,u,v_node;

    for (i = 0; i < igraph_vector_int_size(&bush->topological_order); ++i) {
        u = VECTOR(bush->topological_order)[i];
        // Se 'u' não é alcançável a partir de 'from' (para menor caminho),
        if (VECTOR(bush->paths.dist_shortest_local)[u] == DBL_MAX) continue;
        
        igraph_incident(&bush->Grafo, &incident_edges, u, IGRAPH_OUT,IGRAPH_NO_LOOPS);
        
        for (j = 0; j < igraph_vector_int_size(&incident_edges); ++j) {

            eid = VECTOR(incident_edges)[j];
            v_node = IGRAPH_TO(&bush->Grafo, eid);
            id = VECTOR(bush->edge_id)[eid]; // Acessa o id da aresta

            fluxo = VECTOR(bush->flow_per_origin)[id]; // Acessa o fluxo da aresta
            weight = VECTOR(time)[id]; // Acessa o peso da aresta            
            if (VECTOR(bush->paths.dist_shortest_local)[u] + weight < VECTOR(bush->paths.dist_shortest_local)[v_node]) {
                
                VECTOR(bush->paths.dist_shortest_local)[v_node] = VECTOR(bush->paths.dist_shortest_local)[u] + weight;
                VECTOR(bush->paths.min_edges)[v_node] = id;
                
            }
            
            // Relaxamento para o maior caminho
            // Só estender de 'u' se 'u' já tem um caminho mais longo válido (não -DBL_MAX)
            if ((VECTOR(bush->paths.dist_longest_local)[u] != -DBL_MAX) && (fluxo > 0)) {
                if (VECTOR(bush->paths.dist_longest_local)[u] + weight > VECTOR(bush->paths.dist_longest_local)[v_node]) {

                    VECTOR(bush->paths.dist_longest_local)[v_node] = VECTOR(bush->paths.dist_longest_local)[u] + weight;
                    VECTOR(bush->paths.max_edges)[v_node] = id;
                    
                }
            }
        }
        igraph_vector_int_clear(&incident_edges); // Reutiliza o vetor de arestas
    }
    // Libera memória dos vetores locais
    igraph_vector_int_destroy(&incident_edges);
    igraph_vector_destroy(&time);
    

}

int find_lca(
    int fonte,
    int alvo,
    struct BUSH *bush,
    igraph_t *Grafo
) {

    int edge_id,current_node = alvo;
    int next;
    int edge_min_parrent = VECTOR(bush->paths.min_edges)[alvo];
    int edge_max_parrent = VECTOR(bush->paths.max_edges)[alvo];
    if(edge_min_parrent == -1 || edge_max_parrent == -1) return -1; // Se não houver caminho mínimo ou máximo, retorna -1

    int min_parrent = IGRAPH_FROM(Grafo, edge_min_parrent);
    int max_parrent = IGRAPH_FROM(Grafo, edge_max_parrent);


    if(min_parrent == max_parrent) return -1; // Se os pais são iguais, retorna -1

    bool* visited = (bool*) calloc(igraph_vcount(Grafo), sizeof(bool));
    visited[alvo] = true; // Marca o nó alvo como visitado
    visited[fonte] = true; // Marca a fonte como visitada
    while(current_node != fonte) {

        edge_id = VECTOR(bush->paths.min_edges)[current_node];
        current_node = IGRAPH_FROM(Grafo, edge_id);
        visited[current_node] = true; // Marca o nó atual como visitado
        if (edge_id == -1) break; // Se não houver aresta, sai do loop
    }
    current_node = max_parrent;
    bool find = false;
    int next_node;
    while(current_node != fonte) {
        edge_id = VECTOR(bush->paths.max_edges)[current_node];
        next_node = IGRAPH_FROM(Grafo, edge_id);
        if (edge_id == -1){
            free(visited); // Libera a memória alocada para o vetor visited
            return -1; // Se não houver aresta, retorna -1
        }
        //current_node = IGRAPH_FROM(Grafo, edge_id);
        if ((visited[current_node]) && (current_node != alvo)) {
            free(visited); // Libera a memória alocada para o vetor visited
            if(current_node != fonte) return current_node; // Retorna o nó onde os caminhos divergem
            return -1; // Se não encontrar um LCA válido, retorna -1
        }
        current_node = next_node;
    }
    if((current_node == fonte) && (visited[current_node])){
        free(visited); // Libera a memória alocada para o vetor visited
        return current_node;
    }
    free(visited);
    return -1; // Se não encontrar um LCA válido, retorna -1
}

void shift_flow(
    double delta,
    int alvo,
    int lca,
    struct BUSH *bush,
    const igraph_t *Grafo,
    igraph_vector_t *total_flow
) {
    int current_node;
    int edge_id;

    // Adiciona fluxo ao caminho mínimo
    current_node = alvo;
    while (current_node != lca) {
        edge_id = VECTOR(bush->paths.min_edges)[current_node];
        if (edge_id == -1) break;
        VECTOR(*total_flow)[edge_id] += delta;
        VECTOR(bush->flow_per_origin)[edge_id] += delta; // Atualiza o fluxo da bush
        current_node = IGRAPH_FROM(Grafo, edge_id);
    }

    // Remove fluxo do caminho máximo
    current_node = alvo;
    while (current_node != lca) {
        edge_id = VECTOR(bush->paths.max_edges)[current_node];
        if (edge_id == -1) break;
        VECTOR(*total_flow)[edge_id] -= delta;
        if(VECTOR(*total_flow)[edge_id] < 1e-8) VECTOR(*total_flow)[edge_id] = 0; // Evita fluxo negativo
        VECTOR(bush->flow_per_origin)[edge_id] -= delta; // Atualiza o fluxo da bush
        if(VECTOR(bush->flow_per_origin)[edge_id] < 1e-8) VECTOR(bush->flow_per_origin)[edge_id] = 0; // Evita fluxo negativo
        current_node = IGRAPH_FROM(Grafo, edge_id);
    }
}

/**
 * @brief Verifica se a bush é ótima e, se não, a melhora adicionando arcos.
 * Implementa a Seção 2.5: "Improved bush".
 */
bool improve_bush(
    struct BUSH *bush,
    igraph_t *Grafo,
    struct PARAMETERS *BPR_PARAMETERS,
    int fonte,
    igraph_vector_t *total_flow,
    bool s
) {
    int i,from,to;
    igraph_integer_t index;
    igraph_vector_int_t new_edges;
    igraph_vector_int_init(&new_edges, 0); // Inicializa o vetor de novas arestas
    igraph_vector_t longest;
    max_Distance(&longest, bush, Grafo,fonte, BPR_PARAMETERS, total_flow);
    int L1 = igraph_vector_int_size(&bush->edge_id);
    int L = 0;
    for(i = 0; i < BPR_PARAMETERS->L; i++) {
        from = IGRAPH_FROM(Grafo, i);
        to = IGRAPH_TO(Grafo, i);
        //if(s) printf("Adicionando aresta (%d,%d) - %f - %f %d\n", from+1, to+1, VECTOR(longest)[from], VECTOR(longest)[to],BPR_PARAMETERS->L);
        if (!bush->is_ingraph[i]){
            if(VECTOR(longest)[from] == -DBL_MAX || VECTOR(longest)[to] == -DBL_MAX) continue; // Se não há caminho para a fonte ou destino, pula esta aresta
            if(VECTOR(longest)[from] < VECTOR(longest)[to]) {
                // Adiciona a aresta ao grafo da bush
                bush->is_ingraph[i] = true; // Marca a aresta como parte do grafo da bush
                VECTOR(bush->flow_per_origin)[i] = 0; // Inicializa o fluxo para esta aresta
                igraph_vector_int_push_back(&new_edges, i); // Adiciona a aresta ao vetor de novas arestas
                L++;
            }
        }
        else igraph_vector_int_push_back(&new_edges, i);
        L1 = igraph_vector_int_size(&new_edges);
    }
    if (L == 0) {
        igraph_vector_destroy(&longest);
        igraph_vector_int_destroy(&new_edges);
        return false; // Se não há novas arestas, retorna falso
    }
    igraph_es_t edges;
    igraph_es_vector(&edges,&new_edges);
    igraph_destroy(&bush->Grafo); 
    igraph_subgraph_from_edges(
        Grafo, &bush->Grafo, edges, false
    );
    att_topological_order(bush,Grafo); // Atualiza a ordenação topológica da bush
    igraph_vector_int_destroy(&new_edges); // Libera memória do vetor de novas arestas
    igraph_vector_destroy(&longest); // Libera memória do vetor de maiores distâncias

    // Check if edge IDs match the correct edges in the original graph

    for (long e = 0; e < igraph_ecount(&bush->Grafo); e++) {
        igraph_integer_t from, to;
        igraph_edge(&bush->Grafo, e, &from, &to);
        igraph_integer_t graph_edge_id = VECTOR(bush->edge_id)[e];
        
        // Verify if the edge in the original graph matches
        if (from != IGRAPH_FROM(Grafo, graph_edge_id) || 
            to != IGRAPH_TO(Grafo, graph_edge_id)) {
            printf("Edge mismatch found: Bush edge (%ld,%ld) doesn't match graph edge %ld (%ld,%ld)\n",
                   from, to, graph_edge_id, 
                   IGRAPH_FROM(Grafo, graph_edge_id), 
                   IGRAPH_TO(Grafo, graph_edge_id));
            exit(1);
        }
    }
    return true; // Return true if edges were added
}

void Dial(
    igraph_t *Grafo,
    struct OD_MATRIX* OD,
    struct PARAMETERS* BPR_PARAMETERS,
    igraph_vector_t *solucao,
    struct BUSH **bushes,
    bool WARM_START 
) {

    igraph_vector_init(solucao, BPR_PARAMETERS->L);
    int i, j, k, iter = 1,fonte,target;
    if(WARM_START && bushes != NULL){
        for(i = 0; i < OD->size; i++){
            for(j = 0; j < (*bushes)[i].n_alvos; j++){
                target = VECTOR(OD->Elementos[i].alvos)[j];
                int edge_id;
                while(target != OD->Elementos[i].fonte){
                    edge_id = VECTOR((*bushes)[i].paths.min_edges)[target];

                    VECTOR((*bushes)[i].flow_per_origin)[edge_id] += VECTOR(OD->Elementos[i].warm_volumes)[j];
                    
                    target = IGRAPH_FROM(Grafo, edge_id);
                }
            }
        }
        for(i = 0; i < BPR_PARAMETERS->L; i++){
            
            for(j = 0; j < OD->size; j++){
                VECTOR(*solucao)[i] += VECTOR((*bushes)[j].flow_per_origin)[i];
            }
        }
    }
    else{
        *bushes = (struct BUSH*) malloc(OD->size * sizeof(struct BUSH));
        igraph_bool_t has_multiple;
        for (i = 0; i < OD->size; i++){
    
            (*bushes)[i].n_alvos = igraph_vector_int_size(&OD->Elementos[i].alvos);
            igraph_vector_init(&(*bushes)[i].flow_per_origin, BPR_PARAMETERS->L);
    
            init_bush(&(*bushes)[i], Grafo, i, BPR_PARAMETERS, OD, solucao,i == 147);
            igraph_has_multiple(&(*bushes)[i].Grafo, &has_multiple);
            if (has_multiple) {
                printf("Bush %d has multiple edges!\n", i);
                exit(0); // Encerra o programa para depuração
            }
        }

    }
    // Print initial solution values
    bool ALL_BUSHES_ARE_OPTIMAL = false;
    //FILE *file = fopen("gap.txt", "w");
    double GAP, previous_GAP = DBL_MAX;
    clock_t start_time, end_time;

    for(int iter = 0; iter < MAX_ITER; iter++) {
        start_time = clock();
        //printf("🔢🔢🔢🔢🔢🔢🔢🔢 Iteração Global %d 🔢🔢🔢🔢🔢🔢🔢🔢\n", iter+1);
        ALL_BUSHES_ARE_OPTIMAL = true; // Assume que todas são ótimas nesta iteração
        double max_deltax = 0.;
            
        for (i = 0; i < OD->size; i++){
            //printf("🌳 Melhorando Bush %d de %d 🌳\n", i+1, OD->size);
            fonte = OD->Elementos[i].fonte;
            improve_bush(&(*bushes)[i], Grafo,BPR_PARAMETERS,fonte, solucao,i == 0);
            bool flow_was_shifted = false;

            find_dag_shortest_longest_costs_and_parents(&(*bushes)[i], Grafo, fonte, BPR_PARAMETERS, solucao);

            for (j = 0; j < BPR_PARAMETERS->N; j++) {
                if(j == fonte) continue; // Pula a fonte
                int lca = find_lca(fonte, j, &(*bushes)[i],Grafo);
                if (lca == -1) continue;

                double delta_x = findFlowDelta(&(*bushes)[i],Grafo, lca, j, BPR_PARAMETERS, solucao);

                if(delta_x > EPSILON){
                    if(delta_x > max_deltax) max_deltax = delta_x;
                    //printf("  Delta_x calculado: %f\n", delta_x);
                    shift_flow(delta_x, j, lca, &(*bushes)[i], Grafo, solucao);
                    flow_was_shifted = true; // Indica que um fluxo foi deslocado
                    
                }
                if(delta_x == 0) continue;
                if(delta_x < 0) {
                    printf("  Delta_x negativo encontrado: %e, ajustando para 0\n", delta_x);
                    delta_x = 0; // Ajusta delta_x para 0 se for negativo
                    exit(0); // Encerra o programa se delta_x for negativo
                }
            }
            //printf("  Fluxo deslocado: %s\n", flow_was_shifted ? "Sim" : "Não");
            
            if(flow_was_shifted){
                ALL_BUSHES_ARE_OPTIMAL = false;
                removeUnusedArcs(&(*bushes)[i], Grafo,fonte, BPR_PARAMETERS, solucao);
                //printf(" Ligações Atualizadas!\n");
            }
            
        }
        GAP = relative_gap(solucao, Grafo, BPR_PARAMETERS, OD);
        end_time = clock();

        /* if(fabs(previous_GAP - GAP)/previous_GAP < 1e-6) {
            printf(" (Mudança: %e)\n", fabs(previous_GAP - GAP));
            //printf("Convergência alcançada após %d iterações globais.\n", iter+1);
            break;
        } */
        printf("%d %e %e %e\n", iter+1, GAP, fabs(previous_GAP - GAP)/previous_GAP, ((double)(end_time - start_time)) / CLOCKS_PER_SEC);
        if(GAP < 1e-10) break;
        previous_GAP = GAP;
    }
    printf("\n=== Fluxos Finais ===\n");
    for (int i = 0; i < BPR_PARAMETERS->L; i++) {
        printf("Arco (%ld, %ld): %10.2f\n", IGRAPH_FROM(Grafo,i), IGRAPH_TO(Grafo,i), VECTOR(*solucao)[i]);
    }
}