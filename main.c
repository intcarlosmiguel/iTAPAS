#include "bib/simulate.h"
//#include "bib/ODE.h"
#include <igraph/igraph.h>
#include <omp.h>

int main(int argc, char *argv[]) {
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_setup();
    char *arquivoEDGES = argv[1];
    char *algoritmo = argv[2];
    int N_ODS = atoi(argv[3]);
    int VOLUMES = atoi(argv[4]);
    int i;

    //#pragma omp parallel for num_threads(8)
    for (i = 0; i < 1; i++) {
        char arquivoOD_local[10000];
        if(VOLUMES == 0){
            snprintf(arquivoOD_local, sizeof(arquivoOD_local), "./example/od.txt");
        }else{
            snprintf(arquivoOD_local, sizeof(arquivoOD_local), "./od_outputs/OD_%d_%d/OD_%d.txt", N_ODS, VOLUMES, i);
        }
        printf("Simulating with OD file: %s\n", arquivoOD_local);
        //snprintf(arquivoOD_local, sizeof(arquivoOD_local), "./example/od.txt");
        simulate_simple(arquivoEDGES, arquivoOD_local, algoritmo);
        printf("Simulação %d de %d concluída.\n", i + 1, N_ODS);
    }
    //sprintf(arquivoOD, "./od_outputs/od_%d_%d.txt", N_ODS, VOLUMES);
    //simulate_example(arquivoEDGES, arquivoOD, algoritmo);
}