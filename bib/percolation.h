#pragma once

#include "calc.h"
#include "define.h"
#include <igraph/igraph.h>

static inline int minimo(int a, int b) { return (a < b) ? a : b; }

static inline void tarjan_dfs(const igraph_t *G, igraph_integer_t v,
                              tarjan_ctx_t *ctx) {

  VECTOR(ctx->index)[v] = ctx->index_counter;
  VECTOR(ctx->low)[v] = ctx->index_counter;
  ctx->index_counter++;

  igraph_vector_int_push_back(&ctx->stack, v);
  VECTOR(ctx->onstack)[v] = 1;

  igraph_vector_int_t edges;
  igraph_vector_int_init(&edges, 0);

  igraph_incident(G, &edges, v, IGRAPH_OUT, IGRAPH_NO_LOOPS);

  for (igraph_integer_t i = 0; i < igraph_vector_int_size(&edges); i++) {

    igraph_integer_t eid = VECTOR(edges)[i];

    /* aresta removida */
    if (!VECTOR(*(ctx->use_edge))[eid]) {
      continue;
    }

    igraph_integer_t from, to;
    igraph_edge(G, eid, &from, &to);

    igraph_integer_t w = to;

    if (VECTOR(ctx->index)[w] == -1) {
      tarjan_dfs(G, w, ctx);
      VECTOR(ctx->low)[v] = minimo(VECTOR(ctx->low)[v], VECTOR(ctx->low)[w]);
    } else if (VECTOR(ctx->onstack)[w]) {
      VECTOR(ctx->low)[v] = minimo(VECTOR(ctx->low)[v], VECTOR(ctx->index)[w]);
    }
  }

  igraph_vector_int_destroy(&edges);

  /* v é raiz de uma CFC */
  if (VECTOR(ctx->low)[v] == VECTOR(ctx->index)[v]) {

    igraph_integer_t w;
    do {
      w = igraph_vector_int_pop_back(&ctx->stack);
      VECTOR(ctx->onstack)[w] = 0;
      VECTOR(ctx->comp)[w] = ctx->comp_counter;
    } while (w != v);

    ctx->comp_counter++;
  }
}

static inline igraph_error_t tarjan_scc(const igraph_t *G,
                                        const igraph_vector_bool_t *use_edge,
                                        igraph_vector_int_t *comp_out,
                                        int *num_components) {

  igraph_integer_t N = igraph_vcount(G);
  igraph_integer_t M = igraph_ecount(G);

  /* verificação de consistência */
  if (igraph_vector_bool_size(use_edge) != M) {
    return IGRAPH_EINVAL;
  }

  tarjan_ctx_t ctx;

  igraph_vector_int_init(&ctx.index, N);
  igraph_vector_int_init(&ctx.low, N);
  igraph_vector_bool_init(&ctx.onstack, N);
  igraph_vector_int_init(&ctx.comp, N);
  igraph_vector_int_init(&ctx.stack, 0);

  for (igraph_integer_t v = 0; v < N; v++) {
    VECTOR(ctx.index)[v] = -1;
    VECTOR(ctx.low)[v] = -1;
    VECTOR(ctx.comp)[v] = -1;
    VECTOR(ctx.onstack)[v] = 0;
  }

  ctx.index_counter = 0;
  ctx.comp_counter = 0;
  ctx.use_edge = use_edge;

  for (igraph_integer_t v = 0; v < N; v++) {
    if (VECTOR(ctx.index)[v] == -1) {
      tarjan_dfs(G, v, &ctx);
    }
  }

  igraph_vector_int_init_copy(comp_out, &ctx.comp);

  if (num_components != NULL) {
    *num_components = ctx.comp_counter;
  }

  igraph_vector_int_destroy(&ctx.index);
  igraph_vector_int_destroy(&ctx.low);
  igraph_vector_bool_destroy(&ctx.onstack);
  igraph_vector_int_destroy(&ctx.stack);
  igraph_vector_int_destroy(&ctx.comp);

  return IGRAPH_SUCCESS;
}

static inline void init_percolation_state(PercolationState *state, int N,
                                          int L) {
  state->congested_edges = (bool *)calloc(L, sizeof(bool));
  state->nodes_in_giant_component = (bool *)calloc(N, sizeof(bool));
  state->giant_component_id = -1;
  state->giant_component_fraction = 0.0;
  state->second_component_fraction = 0.0;
  state->finite_cluster_mean = 0.0;
  state->susceptibility = 0.0;
}

static inline int compara_Sorting(const void *a, const void *b) {
  Sorting *sa = (Sorting *)a;
  Sorting *sb = (Sorting *)b;
  if (sa->value < sb->value)
    return 1;
  else if (sa->value > sb->value)
    return -1;
  else
    return 0;
}

static inline void compute_strongly_connected_components(
    PercolationState *state, PercolationState *prev_state, igraph_t *graph, const char *folder_path, int mode, int random_seed,int phase) {
  igraph_vector_int_t neighbors;
  igraph_vector_int_init(&neighbors, 0);
  igraph_neighbors(graph, &neighbors, 180, IGRAPH_OUT, IGRAPH_NO_LOOPS,
                   IGRAPH_NO_MULTIPLE);

  igraph_vector_int_destroy(&neighbors);
  
  if (state->probabilty_congested_edges > 0) {
    int i;
    int total_edge_count = igraph_ecount(graph);
    int total_node_count = igraph_vcount(graph);
    igraph_vector_int_t components;
    int number_of_components;
    igraph_vector_int_init(&components, 0);
    igraph_vector_bool_t ligacoes_ativas;
    igraph_vector_bool_init(&ligacoes_ativas, total_edge_count);

    for (i = 0; i < total_edge_count; i++) {
      VECTOR(ligacoes_ativas)[i] = !state->congested_edges[i];
    }
    for (i = 0; i < total_node_count; i++) {
      state->nodes_in_giant_component[i] = false;
    }

    tarjan_scc(graph, &ligacoes_ativas, &components, &number_of_components);

    /* Step 4: Calculate cluster metrics */
    if (number_of_components > 0) {

      Sorting *clusters =
          (Sorting *)malloc(number_of_components * sizeof(Sorting));

      state->susceptibility = 0.0;
      state->finite_cluster_mean = 0.0;

      for (i = 0; i < number_of_components; i++) {
        clusters[i].id = i;
        clusters[i].value = 0;
      }

      for (i = 0; i < total_node_count; i++) {
        int comp_id = VECTOR(components)[i];
        clusters[comp_id].value += 1;
      }

      qsort(clusters, number_of_components, sizeof(Sorting), compara_Sorting);
      state->giant_component_fraction =
          (double)clusters[0].value / total_node_count;

      for (i = 0; i < number_of_components; i++)
        state->finite_cluster_mean += clusters[i].value;

      state->finite_cluster_mean /= number_of_components;

      if (number_of_components > 1) {
        state->second_component_fraction =
            (double)clusters[1].value / total_node_count;

        for (i = 1; i < number_of_components; i++)
          state->susceptibility += clusters[i].value * clusters[i].value;

        state->susceptibility /= number_of_components;
        state->susceptibility /= state->finite_cluster_mean;
      } else
        state->second_component_fraction = 0.0;

      bool is_critical_point = false;
      double prev_scc2 = prev_state->susceptibility;
      /* Check for peak conditions */
      if (prev_scc2 > 0 && state->susceptibility / prev_scc2 > 10)
        is_critical_point = state->giant_component_fraction < 0.9;

      if (random_seed == SEED_DEBUG_ID && phase == PHASE_ANALYSIS && is_critical_point) {
        printf("Saving files...\n");
        char debug_map_path[1024];
        char mode_name[32];
        get_percolation_mode_name(mode, mode_name);
        snprintf(debug_map_path, 1024, "%s/node_cluster_%s.txt", folder_path,
                 mode_name);

        FILE *df = fopen(debug_map_path, "w");
        if (df != NULL) {
          for (int v = 0; v < total_node_count; v++)
            fprintf(df, "%d %ld\n", v, VECTOR(components)[v]);
          fclose(df);
        }
        snprintf(debug_map_path, 1024, "%s/congested_edges_%s.txt", folder_path,
                 mode_name);
        FILE *ef = fopen(debug_map_path, "w");
        if (ef != NULL) {
          for (int e = 0; e < total_edge_count; e++)
            fprintf(ef, "%d %d\n", e, state->congested_edges[e]);
          fclose(ef);
        }
      }
      igraph_vector_bool_destroy(&ligacoes_ativas);
      igraph_vector_int_destroy(&components);
      free(clusters);
    }
  } else {
    state->giant_component_fraction = 1.0;
    state->second_component_fraction = 0.0;
    state->finite_cluster_mean = 0.0;
    state->susceptibility = 0.0;
  }
}