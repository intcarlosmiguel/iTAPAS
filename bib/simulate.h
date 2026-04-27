#pragma once
#include "calc.h"
#include "define.h"
#include "cfw.h"
#include "dial.h"
#include "greedy_path.h"
#include "iTAPAS.h"
#include "leblanc.h"
#include "mtwister.h"
#include "percolation.h"
#include <igraph/igraph.h>

#include <errno.h>
#include <string.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir_dir(path) _mkdir(path)
#else
#define mkdir_dir(path) mkdir(path, 0755)
#endif

static inline void check_or_create_dir(const char *path) {
  struct stat st = {0};

  // Checa se a pasta existe
  if (stat(path, &st) == -1) {
    // Pasta não existe, tenta criar
    if (mkdir_dir(path) == 0) {
      printf("Diretório '%s' criado com sucesso.\n", path);
    }
  }
}

static inline void init_time(double *initial_time,
                             struct PARAMETERS BPR_PARAMETERS,
                             struct OD_MATRIX *OD_MATRIX, igraph_t *Grafo) {
  igraph_vector_t time;
  igraph_vector_t free_flow;
  igraph_vector_init(&time, BPR_PARAMETERS.L);
  igraph_vector_init(&free_flow, BPR_PARAMETERS.L);
  igraph_vector_update(&free_flow, &BPR_PARAMETERS.capacidade);
  BPR(&time, &BPR_PARAMETERS, &free_flow);
  int index = 0;
  for (int i = 0; i < OD_MATRIX->size; i++) {
    int fonte = OD_MATRIX->Elementos[i].fonte;
    igraph_matrix_t res;
    igraph_matrix_init(&res, 0, 0);
    igraph_distances(Grafo, &time, &res, igraph_vss_1(fonte), igraph_vss_all(),
                     IGRAPH_OUT);
    for (int j = 0; j < igraph_vector_int_size(&OD_MATRIX->Elementos[i].alvos);
         j++) {
      int node = VECTOR(OD_MATRIX->Elementos[i].alvos)[j];
      initial_time[index] = MATRIX(res, 0, node);
      index++;
    }
    igraph_matrix_destroy(&res);
  }
  igraph_vector_destroy(&time);
  igraph_vector_destroy(&free_flow);
}

static inline void init_simulate(struct PARAMETERS *BPR_PARAMETERS,
                                 struct OD_MATRIX *OD, igraph_vector_t *solucao,
                                 igraph_t *Grafo,
                                 const char *algoritmo) {
  if (strcmp(algoritmo, "Leblanc") == 0) {
    // printf("Using Leblanc's algorithm for shortest paths.\n");
    leblanc(BPR_PARAMETERS, OD, Grafo, solucao, false, NULL, NULL);
  } else if (strcmp(algoritmo, "CFW") == 0) {
    printf("Using CFW algorithm for shortest paths.\n");
    cfw_solver(BPR_PARAMETERS, OD, Grafo, solucao, false, NULL, NULL);
    
  } else if (strcmp(algoritmo, "Dial") == 0) {

    struct BUSH *bushes;
    printf("Using Dial's algorithm for shortest paths.\n");
    Dial(Grafo, OD, BPR_PARAMETERS, solucao, &bushes, false, false);
    for (int i = 0; i < OD->size; i++)
      free_bush(&bushes[i]);
    free(bushes);
    
  } else if (strcmp(algoritmo, "iTAPAS") == 0) {
    // printf("Using iTAPAS's algorithm for shortest paths.\n");
    // iTAPAS(BPR_PARAMETERS, OD, Grafo, solucao);
  } else {
    printf("Error: Unknown algorithm '%s'\n", algoritmo);
  }
}



static inline void simulate_simple(const char *arquivoEDGES,
                                   const char *arquivoOD,
                                   const char *algoritmo) {

  struct PARAMETERS BPR_PARAMETERS;
  struct OD_MATRIX OD_MATRIX;

  igraph_vector_int_t edges;
  igraph_vector_int_init(&edges, 0);
  init_parameters(&BPR_PARAMETERS, &edges, arquivoEDGES);

  load_OD_from_file(arquivoOD, &OD_MATRIX);

  igraph_t Grafo;
  igraph_empty(&Grafo, BPR_PARAMETERS.N, IGRAPH_DIRECTED);
  igraph_add_edges(&Grafo, &edges, NULL);

  igraph_vector_t solucao;
  init_simulate(&BPR_PARAMETERS, &OD_MATRIX, &solucao, &Grafo, algoritmo);

  igraph_destroy(&Grafo);
  igraph_vector_destroy(&solucao);
  igraph_vector_int_destroy(&edges);
}

static inline void simulate_percolation(const char *arquivoEDGES,
                                        const char *arquivoOD,
                                        const char *algoritmo) {

  struct PARAMETERS BPR_PARAMETERS;
  struct OD_MATRIX OD_MATRIX;

  igraph_vector_int_t edges;
  igraph_vector_int_init(&edges, 0);
  init_parameters(&BPR_PARAMETERS, &edges, arquivoEDGES);

  load_OD_from_file(arquivoOD, &OD_MATRIX);
  // print_OD_matrix(&OD_MATRIX);
  igraph_t Grafo;
  igraph_empty(&Grafo, BPR_PARAMETERS.N, IGRAPH_DIRECTED);
  igraph_add_edges(&Grafo, &edges, NULL);
  PercolationState state;
  PercolationState prev_state;
  init_percolation_state(&state, BPR_PARAMETERS.N, BPR_PARAMETERS.L);
  init_percolation_state(&prev_state, BPR_PARAMETERS.N, BPR_PARAMETERS.L);

  struct BUSH *bushes = NULL;
  igraph_vector_t solucao;
  bool warm_start = false;
  //struct WARM_START WS;
  //initialize_warm_start(&WS, OD_MATRIX.size);

  // Inbounds salvos da iteração anterior (árvore de caminhos mínimos)
  igraph_vector_int_t *saved_inbounds =
      (igraph_vector_int_t *)malloc(OD_MATRIX.size * sizeof(igraph_vector_int_t));
  for (int j = 0; j < OD_MATRIX.size; j++)
    igraph_vector_int_init(&saved_inbounds[j], 0);

  // Custos SP salvos por par OD (para GAP aproximado)
  double *saved_sp_costs = (double *)calloc(OD_MATRIX.n_elements, sizeof(double));

  double demand_increment = 1.0;
  for (int i = 0; i < 563; i++) {

    if (warm_start) {
      // Incrementar volumes ANTES para que volumes == fluxo total real
      for (int j = 0; j < OD_MATRIX.size; j++) {
        int n_alvos = igraph_vector_int_size(&OD_MATRIX.Elementos[j].alvos);
        for (int k = 0; k < n_alvos; k++) {
          VECTOR(OD_MATRIX.Elementos[j].volumes)[k] += demand_increment;
          VECTOR(OD_MATRIX.Elementos[j].warm_volumes)[k] = demand_increment;
        }
      }
      // Carrega warm_volumes usando os inbounds salvos (sem Dijkstra)
      load_warm_volumes_from_inbounds(&OD_MATRIX, &Grafo, &solucao, saved_inbounds);
    }
    printf("Starting iteration %d with total demand: %f (increment: %.0f)\n",
           i + 1, VECTOR(OD_MATRIX.Elementos[0].volumes)[0], demand_increment);

    //iTAPAS(&BPR_PARAMETERS, &OD_MATRIX, &Grafo, &solucao, &WS);
    //Dial(&Grafo, &OD_MATRIX, &BPR_PARAMETERS, &solucao, &bushes, warm_start, true);
    leblanc(&BPR_PARAMETERS, &OD_MATRIX, &Grafo, &solucao, warm_start, saved_inbounds, saved_sp_costs);
    double p = 0.;
    for (int j = 0; j < BPR_PARAMETERS.L; j++) {
      double time =
          single_BPR(VECTOR(solucao)[j], VECTOR(BPR_PARAMETERS.cost_time)[j],
                     VECTOR(BPR_PARAMETERS.capacidade)[j]);
      if (time > 3 * VECTOR(BPR_PARAMETERS.cost_time)[j]) {
        p++;
      }
    }
    p /= BPR_PARAMETERS.L;
    state.probabilty_congested_edges = p;
    compute_strongly_connected_components(&state,&prev_state, &Grafo,NULL,0,0,1);

    warm_start = true;
    demand_increment += 3.0;
  }

  // Cleanup
  for (int i = 0; i < OD_MATRIX.size; i++)
    igraph_vector_int_destroy(&saved_inbounds[i]);
  free(saved_inbounds);
  free(saved_sp_costs);
  if (bushes != NULL) {
    for (int i = 0; i < OD_MATRIX.size; i++) {
      free_bush(&bushes[i]);
    }
    free(bushes);
  }
  igraph_vector_destroy(&solucao);
}