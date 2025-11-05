#include "bib/simulate.h"
#include "bib/ODE.h"
#include "igraph.h"

int main(int argc, char *argv[]) {
    igraph_set_attribute_table(&igraph_cattribute_table);
    char *arquivoEDGES = argv[1];
    char *arquivoOD = argv[2];
    char *algoritmo = argv[3];
    simulate_example(arquivoEDGES, arquivoOD, algoritmo);
}