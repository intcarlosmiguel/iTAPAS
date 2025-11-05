#pragma once

#include <define.h>
#include "igraph.h"
#include <ctype.h>

void print_vector_igraph(igraph_vector_int_t* vetor){
    int N =igraph_vector_int_size(vetor);
    for (int i = 0; i < N; i++){
        if(i != N - 1) printf("%ld,",VECTOR(*vetor)[i]);
        else printf("%ld\n",VECTOR(*vetor)[i]);
    }
    
}

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

void load_OD_from_file(const char* filename, struct OD_MATRIX* OD_MATRIX) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file: %s\n", filename);
        return;
    }
    int current_source = -1;

    int source, target;
    double flow;
    OD_MATRIX->size = 0;
    OD_MATRIX->Elementos = malloc(sizeof(struct ElementOD) * 0);
    OD_MATRIX->n_elements = 0;
    while (fscanf(file, "%d %d %lf", &source, &target, &flow) == 3) {
        OD_MATRIX->n_elements++;
        if (current_source != source - 1) {
            OD_MATRIX->Elementos = realloc(OD_MATRIX->Elementos, sizeof(struct ElementOD) * (OD_MATRIX->size + 1));
            OD_MATRIX->Elementos[OD_MATRIX->size].fonte = source - 1;

            igraph_vector_int_init(&OD_MATRIX->Elementos[OD_MATRIX->size].alvos, 0);
            igraph_vector_init(&OD_MATRIX->Elementos[OD_MATRIX->size].volumes, 0);
            igraph_vector_init(&OD_MATRIX->Elementos[OD_MATRIX->size].warm_volumes, 0);

            igraph_vector_int_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size].alvos, target - 1);
            igraph_vector_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size].volumes, flow);
            igraph_vector_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size].warm_volumes, 0.0);
            OD_MATRIX->size++;
            current_source = source - 1;
        }
        else{
            igraph_vector_int_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size - 1].alvos, target - 1);
            igraph_vector_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size - 1].volumes, flow);
            igraph_vector_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size - 1].warm_volumes, 0.0);
        }
    }

    fclose(file);
}

void print_OD_matrix(struct OD_MATRIX* OD_MATRIX) {
    printf("Origin-Destination Matrix:\n");
    for (int i = 0; i < OD_MATRIX->size; i++) {
        printf("Origin %d -> Destinations: ", OD_MATRIX->Elementos[i].fonte+1);
        for (int j = 0; j < igraph_vector_int_size(&OD_MATRIX->Elementos[i].alvos); j++) {
            printf("(D:%ld, V:%f) ", VECTOR(OD_MATRIX->Elementos[i].alvos)[j]+1, VECTOR(OD_MATRIX->Elementos[i].volumes)[j]);
        }
        printf("\n");
    }
}



void free_bush(struct BUSH *bush) {
    igraph_vector_int_destroy(&bush->edge_id);
    igraph_vector_int_destroy(&bush->topological_order);
    igraph_vector_destroy(&bush->flow_per_origin);
    free(bush->is_ingraph);
    igraph_destroy(&bush->Grafo);
    igraph_vector_destroy(&bush->paths.dist_longest_local);
    igraph_vector_destroy(&bush->paths.dist_shortest_local);
    igraph_vector_int_destroy(&bush->paths.max_edges);
    igraph_vector_int_destroy(&bush->paths.min_edges);

}

void print_edges(
    int fonte,
    int alvo,
    igraph_vector_int_t *edges,
    igraph_t *Grafo
){
    int k = alvo,edge_id;
    printf("%d -> ",k);
    while(k != fonte) {
        edge_id = VECTOR(*edges)[k];
        printf("%ld -> ", IGRAPH_FROM(Grafo, edge_id));
        k = IGRAPH_FROM(Grafo, edge_id);
    }
    printf("\n");
}

void print_flow(
    igraph_vector_t *flow,
    igraph_t *Grafo,
    char* filename
) {
    if (igraph_ecount(Grafo) != igraph_vector_size(flow)) {
        printf("Error: Number of edges (%ld) does not match flow vector size (%ld)\n", 
               igraph_ecount(Grafo), igraph_vector_size(flow));
        exit(1);
    }
    if (filename == NULL) {
        printf("Flow values:\n");
        for (long i = 0; i < igraph_ecount(Grafo); i++) {
            igraph_integer_t from, to;
            igraph_edge(Grafo, i, &from, &to);
            printf("%ld %ld %f\n", from, to, VECTOR(*flow)[i]);
        }
        printf("\n");
    }
    else{
        FILE *file = fopen(filename, "w");
        for (long i = 0; i < igraph_ecount(Grafo); i++) {
            igraph_integer_t from, to;
            igraph_edge(Grafo, i, &from, &to);
            fprintf(file,"%ld %ld %f\n", from, to, VECTOR(*flow)[i]);
        }
        fclose(file);

    }
}

void init_path(struct min_max_bush *path, int N){
    igraph_vector_int_init(&path->min_edges, N);
    igraph_vector_int_init(&path->max_edges, N);
    
    igraph_vector_init(&path->dist_shortest_local, N);
    igraph_vector_fill(&path->dist_shortest_local, DBL_MAX);
    igraph_vector_init(&path->dist_longest_local, N);
    igraph_vector_fill(&path->dist_longest_local, -DBL_MAX); // "Infinito negativo"
    
}

void erase_path(struct min_max_bush *path){
    igraph_vector_int_destroy(&path->min_edges);
    igraph_vector_int_destroy(&path->max_edges);
}

void print_vetor(void* array,int N,int check){
    if(check == sizeof(int)){
        int* intArray = (int*)array;
        for (int i = 0; i < N; i++){
            if(i!=N-1) printf("%d ",intArray[i]);
            else printf("%d\n",intArray[i]);
        }
    }
    if(check == sizeof(double)){
        double* doubleArray = (double*)array;
        for (int i = 0; i < N; i++){
            if(i!=N-1) printf("%.2f ",doubleArray[i]);
            else printf("%.2f\n",doubleArray[i]);
        }
    }
}

int contarLinhasNoArquivo(const char *nomeArquivo) {
    FILE *arquivo = fopen(nomeArquivo, "r");
    if (!arquivo) {
        perror("Erro ao abrir o arquivo");
        return -1; // Retorna -1 em caso de erro
    }

    char buffer[1024]; // Um buffer para armazenar cada linha
    int linhasUteis = 0;

    // Lê o arquivo linha por linha
    while (fgets(buffer, sizeof(buffer), arquivo) != NULL) {
        int ehUtil = 0; // Flag para marcar se a linha é útil
        
        // Itera sobre a linha lida para ver se tem algo além de espaços
        for (int i = 0; buffer[i] != '\0'; i++) {
            // isspace() verifica se o caractere é espaço, tab, newline, etc.
            if (!isspace((unsigned char)buffer[i])) {
                ehUtil = 1; // Encontrou um caractere não-espaço
                break;      // Já sabemos que a linha é útil, podemos parar de verificar
            }
        }

        if (ehUtil) {
            linhasUteis++;
        }
    }

    fclose(arquivo);
    return linhasUteis;
}

double** lerArquivo(const char *nomeArquivo, int nColunas,int* size) {
    int N = contarLinhasNoArquivo(nomeArquivo);
    double** data = (double**) malloc(N*sizeof(double*));
    
    FILE *arquivo;
    char buffer[1024];

    arquivo = fopen(nomeArquivo, "r");
    if (!arquivo) {
        perror("Erro ao abrir o arquivo");
        exit(0);
    }

    // Lê o arquivo linha por linha
    int linha = 0,i;
    while (fgets(buffer, 1024, arquivo) != NULL) {
        // Processa cada linha do arquivo aqui. Neste exemplo, vamos assumir que os dados são números inteiros.
        char *token = strtok(buffer, " "); // Supondo que os dados sejam separados por espaços. Ajuste o delimitador conforme necessário.
        data[linha] =(double*) malloc(nColunas*sizeof(double));
        for (i = 0; i < nColunas && token != NULL; i++) {
            // Converte o token (string) para o tipo de dado desejado, neste caso, int.
            sscanf(token, "%lf", &data[linha][i]);
            
            // Processa o dado da coluna aqui. Neste exemplo, estamos apenas imprimindo.
            //printf("Dado da coluna %d: %f\n", i + 1, data[linha][i]);

            // Avança para o próximo token (próxima coluna)
            token = strtok(NULL, " ");
        }
        linha++;
    }

    // Fecha o arquivo ao terminar de processar.
    fclose(arquivo);
    *size = N;
    return data;
}
void init_parameters(struct PARAMETERS* BPR_PARAMETERS, igraph_vector_int_t* edges, const char* nomeDoArquivo) {

    int i,e1,e2;
    double capacidade,tempo;

    int N = contarLinhasNoArquivo(nomeDoArquivo);

    igraph_vector_init(&BPR_PARAMETERS->capacidade, 0);
    igraph_vector_init(&BPR_PARAMETERS->cost_time, 0);
    

    BPR_PARAMETERS->L = N;
    BPR_PARAMETERS->N = 0;
    FILE *arquivo;
    arquivo = fopen(nomeDoArquivo, "r");
    if (!arquivo) {
        perror("Erro ao abrir o arquivo");
        exit(0);
    }
    for (i = 0; i < N; i++){
        if(fscanf(arquivo, "%d %d %lf %lf\n", &e1, &e2, &capacidade, &tempo));
        igraph_vector_int_push_back(edges,e1-1);
        igraph_vector_int_push_back(edges,e2-1);
        igraph_vector_push_back(&BPR_PARAMETERS->capacidade,capacidade);
        igraph_vector_push_back(&BPR_PARAMETERS->cost_time,tempo);
        if(BPR_PARAMETERS->N < e1) BPR_PARAMETERS->N = e1;
        if(BPR_PARAMETERS->N < e2) BPR_PARAMETERS->N = e2;


    }
    fclose(arquivo);

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
