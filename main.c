#include "bib/simulate.h"
#include "bib/ODE.h"
#include "igraph.h"

int main(int argc, char *argv[]) {
    igraph_set_attribute_table(&igraph_cattribute_table);
    switch (argc){
    case 3:{
        char *arquivoEDGES = argv[1];
        char *arquivoOD = argv[2];
        simulate_example(arquivoEDGES, arquivoOD);

        break;
    }
    case 4:{

    
        char *arquivoEDGES = argv[1];
        int K = atoi(argv[2]);
        int seed = atoi(argv[3]);
        
        char arquivoOD[100];
        sprintf(arquivoOD, "./fortaleza/OD_0/%d/od_%d.txt",K, seed);
        simulate_percolation(arquivoEDGES, arquivoOD,K,seed);

        break;
    }
    case 6:
        char *file_edges = argv[1];
        char *file_flow_hat = argv[2];
        double lambda = atof(argv[3]);
        double eta = atof(argv[4]);
        int iter = atoi(argv[5]);
        OD_ESTIMATION(lambda, eta, file_edges, file_flow_hat, iter);
        break;
    
    default:
        printf("Error: Insufficient parameters\n");
        printf("Usage: %s <file1> <file2> for simulate_example\n", argv[0]);
        printf("Or more parameters for OD_ESTIMATION\n");
        break;
    }
}