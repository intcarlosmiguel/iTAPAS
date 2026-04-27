#pragma once

#include "define.h"
#include <ctype.h>
#include <igraph/igraph.h>

#include <dirent.h>
#include <errno.h>
#include <limits.h>
#include <linux/limits.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

/* Conta quantos arquivos regulares existem dentro da pasta 'path'.
 * Retorna o número de arquivos (>=0) em caso de sucesso, ou -1 em erro
 * (errno é preservado conforme opendir/stat).
 */
static inline int count_files_in_dir(const char *path) {
  DIR *d = opendir(path);
  if (!d) {
    printf("Erro ao abrir o diretório '%s': %s\n", path, strerror(errno));
    exit(1);
    return -1;
  }

  struct dirent *entry;
  struct stat st;
  char fullpath[PATH_MAX];
  int count = 0;

  while ((entry = readdir(d)) != NULL) {
    if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
      continue;

#if defined(DT_REG)
    if (entry->d_type == DT_REG) {
      count++;
      continue;
    }
#endif

    /* Se o tipo não for informado (DT_UNKNOWN) ou DT_REG não estiver
       disponível, usamos stat para determinar se é um arquivo regular. */
    if (snprintf(fullpath, PATH_MAX, "%s/%s", path, entry->d_name) >=
        PATH_MAX) {
      /* caminho muito longo: considera erro */
      closedir(d);
      errno = ENAMETOOLONG;
      return -1;
    }

    if (stat(fullpath, &st) == 0) {
      if (S_ISREG(st.st_mode))
        count++;
    } else {
      /* ignora entradas que não puderam ser 'stat'-adas */
      continue;
    }
  }

  closedir(d);
  return count;
}

static inline void print_vector_igraph(igraph_vector_int_t *vetor) {
  int N = igraph_vector_int_size(vetor);
  for (int i = 0; i < N; i++) {
    if (i != N - 1)
      printf("%ld,", VECTOR(*vetor)[i]);
    else
      printf("%ld\n", VECTOR(*vetor)[i]);
  }
}

static inline void BPR(igraph_vector_t *tempo, struct PARAMETERS *BPR_PARAMETERS,
         igraph_vector_t *fluxo) {
  double s;
  for (int i = 0; i < BPR_PARAMETERS->L; i++) {
    s = pow(VECTOR(*fluxo)[i] / VECTOR(BPR_PARAMETERS->capacidade)[i], BETA);
    VECTOR(*tempo)
    [i] = (VECTOR(BPR_PARAMETERS->cost_time)[i]) * (1 + ALPHA * s);
  }
}

static inline double single_BPR(double fluxo, double free_flow_time, double capacidade) {
  double s = pow(fluxo / capacidade, BETA);
  if (isnan(s)) {
    printf("Warning: NaN detected in single_BPR\n");
    exit(0);
    s = 0.0;
  }
  return free_flow_time * (1 + ALPHA * s);
}

static inline double single_BPR_derivate(double fluxo, double free_flow_time,
                           double capacidade) {
  double s = pow(fluxo / capacidade, BETA - 1);
  if (isnan(s)) {
    printf("Warning: NaN detected in single_BPR_derivate\n");
    exit(0);
    s = 0.0;
  }
  return free_flow_time * ALPHA * BETA * s / capacidade;
}

/* OPT-2: Fused BPR cost and derivative calculation
 * Computes both values with a single pow() call, reducing math operations.
 * Uses pow(x, BETA-1) then multiplies by x to get pow(x, BETA).
 */
static inline void single_BPR_fused(double fluxo, double free_flow_time,
                                    double capacidade, double *tempo_out,
                                    double *derivada_out) {
  double ratio = fluxo / capacidade;
  double s_beta_minus_1 = pow(ratio, BETA - 1); /* Single pow() call */
  double s_beta = s_beta_minus_1 * ratio;

  if (isnan(s_beta)) {
    *tempo_out = free_flow_time;
    *derivada_out = 0.0;
    return;
  }

  *tempo_out = free_flow_time * (1.0 + ALPHA * s_beta);
  *derivada_out = free_flow_time * ALPHA * BETA * s_beta_minus_1 / capacidade;
}

static inline void BPR_derivate(igraph_vector_t *derivate, struct PARAMETERS *BPR_PARAMETERS,
                  igraph_vector_t *fluxo) {
  double s;
  for (int i = 0; i < BPR_PARAMETERS->L; i++) {
    s = pow(VECTOR(*fluxo)[i] / VECTOR(BPR_PARAMETERS->capacidade)[i],
            BETA - 1);
    if (isnan(s)) {
      printf("Warning: NaN detected in BPR_derivate\n");
      exit(0);
      s = 0.0;
    }
    VECTOR(*derivate)
    [i] = (VECTOR(BPR_PARAMETERS->cost_time)[i]) * ALPHA * BETA * s /
          VECTOR(BPR_PARAMETERS->capacidade)[i];
  }
}

static inline void BPR_gradient(igraph_vector_t *gradiente, struct PARAMETERS *BPR_PARAMETERS,
                  igraph_vector_t *fluxo, double *total_time) {
  double s;
  if (total_time != NULL)
    *total_time = 0;
  for (int i = 0; i < BPR_PARAMETERS->L; i++) {
    s = pow(VECTOR(*fluxo)[i] / VECTOR(BPR_PARAMETERS->capacidade)[i], BETA);
    if (isnan(s)) {
      printf("Warning: NaN detected in BPR_derivate\n");
      exit(0);
      s = 0.0;
    }
    if (total_time != NULL)
      *total_time +=
          (VECTOR(BPR_PARAMETERS->cost_time)[i] * VECTOR(*fluxo)[i]) *
          (1.0 + 0.03 * s);
    VECTOR(*gradiente)
    [i] = (VECTOR(BPR_PARAMETERS->cost_time)[i]) * (1 + ALPHA * (1 + BETA) * s);
  }
}

static inline void load_OD_from_file(const char *filename, struct OD_MATRIX *OD_MATRIX) {
  FILE *file = fopen(filename, "r");
  if (!file) {
    printf("Error opening file: %s\n", filename);
    return;
  }
  int current_source = -1;

  int source, target;
  double flow;
  OD_MATRIX->size = 0;
  OD_MATRIX->Elementos =
      (struct ElementOD *)malloc(sizeof(struct ElementOD) * 0);
  OD_MATRIX->n_elements = 0;
  while (fscanf(file, "%d %d %lf", &source, &target, &flow) == 3) {
    OD_MATRIX->n_elements++;
    if (current_source != source - 1) {
      OD_MATRIX->Elementos = (struct ElementOD *)realloc(
          OD_MATRIX->Elementos,
          sizeof(struct ElementOD) * (OD_MATRIX->size + 1));
      OD_MATRIX->Elementos[OD_MATRIX->size].fonte = source - 1;

      igraph_vector_int_init(&OD_MATRIX->Elementos[OD_MATRIX->size].alvos, 0);
      igraph_vector_init(&OD_MATRIX->Elementos[OD_MATRIX->size].volumes, 0);
      igraph_vector_init(&OD_MATRIX->Elementos[OD_MATRIX->size].warm_volumes,
                         0);

      igraph_vector_int_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size].alvos,
                                  target - 1);
      igraph_vector_push_back(&OD_MATRIX->Elementos[OD_MATRIX->size].volumes,
                              flow);
      igraph_vector_push_back(
          &OD_MATRIX->Elementos[OD_MATRIX->size].warm_volumes, 0.0);
      OD_MATRIX->size++;
      current_source = source - 1;
    } else {
      igraph_vector_int_push_back(
          &OD_MATRIX->Elementos[OD_MATRIX->size - 1].alvos, target - 1);
      igraph_vector_push_back(
          &OD_MATRIX->Elementos[OD_MATRIX->size - 1].volumes, flow);
      igraph_vector_push_back(
          &OD_MATRIX->Elementos[OD_MATRIX->size - 1].warm_volumes, 0.0);
    }
  }

  fclose(file);
}

static inline void print_OD_matrix(struct OD_MATRIX *OD_MATRIX) {
  printf("Origin-Destination Matrix:\n");
  for (int i = 0; i < OD_MATRIX->size; i++) {
    printf("Origin %d -> Destinations: ", OD_MATRIX->Elementos[i].fonte + 1);
    for (int j = 0; j < igraph_vector_int_size(&OD_MATRIX->Elementos[i].alvos);
         j++) {
      printf("(D:%ld, V:%f) ", VECTOR(OD_MATRIX->Elementos[i].alvos)[j] + 1,
             VECTOR(OD_MATRIX->Elementos[i].volumes)[j]);
    }
    printf("\n");
  }
}

static inline void free_bush(struct BUSH *bush) {
  igraph_vector_int_destroy(&bush->topological_order);
  igraph_vector_int_destroy(&bush->topological_value);
  igraph_vector_destroy(&bush->flow_per_origin);
  free(bush->is_ingraph);
  igraph_destroy(&bush->Grafo);
  igraph_vector_destroy(&bush->paths.dist_longest_local);
  igraph_vector_destroy(&bush->paths.dist_shortest_local);
  igraph_vector_int_destroy(&bush->paths.max_edges);
  igraph_vector_int_destroy(&bush->paths.min_edges);

  // Liberar merge nodes
  for (int i = 0; i < bush->numMerges; i++) {
    free(bush->merges[i].approaches);
    free(bush->merges[i].approachFlow);
  }
  free(bush->merges);
  free(bush->nodeToMerge);
}

static inline void print_edges(int fonte, int alvo, igraph_vector_int_t *edges,
                 igraph_t *Grafo) {
  int k = alvo, edge_id;
  printf("%d -> ", k);
  while (k != fonte) {
    edge_id = VECTOR(*edges)[k];
    printf("%ld -> ", IGRAPH_FROM(Grafo, edge_id));
    k = IGRAPH_FROM(Grafo, edge_id);
  }
  printf("\n");
}

static inline void print_flow(igraph_vector_t *flow, igraph_t *Grafo, char *filename) {
  if (igraph_ecount(Grafo) != igraph_vector_size(flow)) {
    printf(
        "Error: Number of edges (%ld) does not match flow vector size (%ld)\n",
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
  } else {
    FILE *file = fopen(filename, "w");
    for (long i = 0; i < igraph_ecount(Grafo); i++) {
      igraph_integer_t from, to;
      igraph_edge(Grafo, i, &from, &to);
      fprintf(file, "%ld %ld %f\n", from, to, VECTOR(*flow)[i]);
    }
    fclose(file);
  }
}

static inline void init_path(struct min_max_bush *path, int N) {
  igraph_vector_int_init(&path->min_edges, N);
  igraph_vector_int_init(&path->max_edges, N);

  igraph_vector_init(&path->dist_shortest_local, N);
  igraph_vector_fill(&path->dist_shortest_local, DBL_MAX);
  igraph_vector_init(&path->dist_longest_local, N);
  igraph_vector_fill(&path->dist_longest_local,
                     -DBL_MAX); // "Infinito negativo"
}

static inline void erase_path(struct min_max_bush *path) {
  igraph_vector_int_destroy(&path->min_edges);
  igraph_vector_int_destroy(&path->max_edges);
}

static inline void print_vetor(void *array, int N, int check) {
  if (check == sizeof(int)) {
    int *intArray = (int *)array;
    for (int i = 0; i < N; i++) {
      if (i != N - 1)
        printf("%d ", intArray[i]);
      else
        printf("%d\n", intArray[i]);
    }
  }
  if (check == sizeof(double)) {
    double *doubleArray = (double *)array;
    for (int i = 0; i < N; i++) {
      if (i != N - 1)
        printf("%.2f ", doubleArray[i]);
      else
        printf("%.2f\n", doubleArray[i]);
    }
  }
}

static inline int contarLinhasNoArquivo(const char *nomeArquivo) {
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

static inline double **lerArquivo(const char *nomeArquivo, int nColunas, int *size) {
  int N = contarLinhasNoArquivo(nomeArquivo);
  double **data = (double **)malloc(N * sizeof(double *));

  FILE *arquivo;
  char buffer[1024];

  arquivo = fopen(nomeArquivo, "r");
  if (!arquivo) {
    perror("Erro ao abrir o arquivo");
    exit(0);
  }

  // Lê o arquivo linha por linha
  int linha = 0, i;
  while (fgets(buffer, 1024, arquivo) != NULL) {
    // Processa cada linha do arquivo aqui. Neste exemplo, vamos assumir que os
    // dados são números inteiros.
    char *token = strtok(
        buffer, " "); // Supondo que os dados sejam separados por espaços.
                      // Ajuste o delimitador conforme necessário.
    data[linha] = (double *)malloc(nColunas * sizeof(double));
    for (i = 0; i < nColunas && token != NULL; i++) {
      // Converte o token (string) para o tipo de dado desejado, neste caso,
      // int.
      sscanf(token, "%lf", &data[linha][i]);

      // Processa o dado da coluna aqui. Neste exemplo, estamos apenas
      // imprimindo.
      // printf("Dado da coluna %d: %f\n", i + 1, data[linha][i]);

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
static inline void init_parameters(struct PARAMETERS *BPR_PARAMETERS,
                     igraph_vector_int_t *edges, const char *nomeDoArquivo) {

  int N = contarLinhasNoArquivo(nomeDoArquivo);

  igraph_vector_init(&BPR_PARAMETERS->capacidade, 0);
  igraph_vector_init(&BPR_PARAMETERS->cost_time, 0);

  BPR_PARAMETERS->L = N;
  BPR_PARAMETERS->N = 0;

  FILE *arquivo = fopen(nomeDoArquivo, "r");
  if (!arquivo) {
    perror("Erro ao abrir o arquivo");
    exit(1);
  }

  char buffer[1024];
  int linha = 0;
  while (fgets(buffer, sizeof(buffer), arquivo) != NULL && linha < N) {
    /* Skip blank lines */
    int blank = 1;
    for (int j = 0; buffer[j] != '\0'; j++) {
      if (!isspace((unsigned char)buffer[j])) { blank = 0; break; }
    }
    if (blank) continue;

    int e1, e2;
    double capacidade, tempo, comprimento = 0.0;

    /* Try 5 columns first, then 4 (mirrors Python's flexible column handling) */
    int ncols = sscanf(buffer, "%d %d %lf %lf %lf", &e1, &e2, &capacidade,
                       &tempo, &comprimento);
    if (ncols < 4) {
      fprintf(stderr, "Erro: linha %d do arquivo '%s' tem %d colunas (mín. 4)\n",
              linha + 1, nomeDoArquivo, ncols);
      continue;
    }

    igraph_vector_int_push_back(edges, e1 - 1);
    igraph_vector_int_push_back(edges, e2 - 1);
    igraph_vector_push_back(&BPR_PARAMETERS->capacidade, capacidade);
    igraph_vector_push_back(&BPR_PARAMETERS->cost_time, tempo);
    if (BPR_PARAMETERS->N < e1)
      BPR_PARAMETERS->N = e1;
    if (BPR_PARAMETERS->N < e2)
      BPR_PARAMETERS->N = e2;

    linha++;
  }

  BPR_PARAMETERS->L = linha; /* Actual number of edges parsed */
  fclose(arquivo);
}
static inline double relative_gap(igraph_vector_t *flow, igraph_t *Grafo,
                    struct PARAMETERS *BPR_PARAMETERS, struct OD_MATRIX *OD,
                    igraph_vector_int_t *out_inbounds,
                    double *out_sp_costs) {

  double gap = 0.0;
  double total_cost = 0.0;
  double total_flow = 0.0;
  int i, j, alvo, edge_id, fonte, size;
  int flat_idx = 0;

  igraph_vector_t time;
  igraph_vector_init(&time, BPR_PARAMETERS->L);
  BPR(&time, BPR_PARAMETERS, flow); // Calcula o tempo de cada aresta
  for (i = 0; i < BPR_PARAMETERS->L; i++)
    total_flow += VECTOR(*flow)[i] * VECTOR(time)[i];

  for (i = 0; i < OD->size; i++) {
    fonte = OD->Elementos[i].fonte;
    size = igraph_vector_int_size(&OD->Elementos[i].alvos);

    igraph_vector_int_t inbound;
    igraph_vector_int_init(&inbound, 0);
    igraph_get_shortest_paths_dijkstra(Grafo, NULL, NULL, fonte,
                                       igraph_vss_all(), &time, IGRAPH_OUT,
                                       NULL, &inbound);

    for (j = 0; j < size; j++) {
      alvo = VECTOR(OD->Elementos[i].alvos)[j];
      double path_cost = 0.0;
      int current = alvo;
      while (current != fonte) {
        edge_id = VECTOR(inbound)[current];
        path_cost += VECTOR(time)[edge_id];
        total_cost += VECTOR(time)[edge_id] * VECTOR(OD->Elementos[i].volumes)[j];
        current = IGRAPH_FROM(Grafo, edge_id);
      }
      if (out_sp_costs != NULL) {
        out_sp_costs[flat_idx] = path_cost;
      }
      flat_idx++;
    }

    // Salva ou destrói o inbound
    if (out_inbounds != NULL) {
      igraph_vector_int_update(&out_inbounds[i], &inbound);
    }
    igraph_vector_int_destroy(&inbound);
  }
  if (total_flow > 0)
    gap = 1 - total_cost / total_flow;
  else
    gap =
        1.0; // Evita divisão por zero, assume que o gap é 1 se não houver fluxo
  igraph_vector_destroy(&time);

  return gap;
}
static inline double Beckman_function(igraph_vector_t *flow, igraph_t *Grafo,
                        struct PARAMETERS *BPR_PARAMETERS) {
  double total_cost = 0.0;
  int i;

  igraph_vector_t time;
  igraph_vector_init(&time, BPR_PARAMETERS->L);
  BPR(&time, BPR_PARAMETERS, flow); // Calcula o tempo de cada aresta
  for (i = 0; i < BPR_PARAMETERS->L; i++)
    total_cost += VECTOR(*flow)[i] * VECTOR(time)[i];

  igraph_vector_destroy(&time);
  return total_cost;
}
static inline void get_percolation_mode_name(int mode, char *buffer) {
  switch (mode) {
  case PERCOLATION_MODE_RECENT:
    strcpy(buffer, "recent");
    break;
  case PERCOLATION_MODE_RATIO:
    strcpy(buffer, "ratio");
    break;
  case PERCOLATION_MODE_FLOW:
    strcpy(buffer, "flow");
    break;
  default:
    strcpy(buffer, "none");
    break;
  }
}


/**
 * GAP aproximado (upper bound) usando custos de SP salvos da iteração anterior.
 * Como BPR é monotoníco, tempos reais >= tempos antigos, logo SP_real >= SP_salvo.
 * SPTT_lower = Σ d_rs_new * SP_salvo(r,s) <= SPTT_real
 * GAP_approx = 1 - SPTT_lower / TSTT_new >= GAP_real
 * Se GAP_approx < ε, seguro pular. Se >= ε, rodar leblanc.
 */
static inline double relative_gap_approximate(igraph_vector_t *flow,
                                              double *saved_sp_costs,
                                              struct PARAMETERS *BPR_PARAMETERS,
                                              struct OD_MATRIX *OD) {
  double sptt_lower = 0.0;
  double tstt = 0.0;
  int i, j;
  int flat_idx = 0;

  igraph_vector_t time;
  igraph_vector_init(&time, BPR_PARAMETERS->L);
  BPR(&time, BPR_PARAMETERS, flow);

  for (i = 0; i < BPR_PARAMETERS->L; i++)
    tstt += VECTOR(*flow)[i] * VECTOR(time)[i];

  for (i = 0; i < OD->size; i++) {
    int n_alvos = igraph_vector_int_size(&OD->Elementos[i].alvos);
    for (j = 0; j < n_alvos; j++) {
      sptt_lower += VECTOR(OD->Elementos[i].volumes)[j] * saved_sp_costs[flat_idx];
      flat_idx++;
    }
  }

  igraph_vector_destroy(&time);

  if (tstt > 0)
    return 1.0 - (sptt_lower / tstt);
  else
    return 1.0;
}
