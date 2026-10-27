/*
  Derived from the fast-loaded-dice-roller-experiments repository of
  F. A. Saad et al. (https://github.com/probsys/fast-loaded-dice-roller-experiments),
  Copyright (C) 2019-2020 MIT Probabilistic Computing Project, F. A. Saad and C. E. Freer,
  released under the Apache License 2.0 (see LICENSES/Apache-2.0.txt).

  The construction is a C port, to arbitrary-precision integers, of the alias
  table of the Rust crate weighted_rand (https://github.com/ichi-h/weighted_rand),
  Copyright (c) 2020 Ippee, released under the MIT License
  (see LICENSES/MIT-weighted_rand.txt).
  Modified in 2026 by Oriane Crouzet, Antoine Genitrini and Mehdi Naima: port of weighted_rand to GMP integers;
  extracted from construct.c and readio.c into a stand-alone file.
  This file is distributed under the GNU General Public License v3.0 or later
  as part of Alias_experiments (see LICENSE).
  SPDX-License-Identifier: GPL-3.0-or-later
*/

/* alias_rust_gmp.c -- construction of the Rust-Alias table on GMP integers. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alias_rust_gmp.h"

// *********************************************************************************
//              ALIAS FROM RUST - GMP (arbitrary-size integers)
// *********************************************************************************

WeightedError weighted_alias_new_gmp(
    alias_rust_gmp_s *out,
    const mpz_t *weights,
    unsigned int n
)
{
    if (n == 0) return WEIGHTED_NO_ITEM;

    vector_init(&out->aliases);
    vector_resize_zero(&out->aliases, n);
    vector_init(&out->small);
    vector_init(&out->large);
    
    vector_mpz_init(&out->prob);
    out->prob.data = malloc(n * sizeof(mpz_t));
    out->prob.size = n;
    out->prob.capacity = n;

    for (unsigned int i = 0; i < n; i++) {
        mpz_init(out->prob.data[i]);
    }

    mpz_init(out->weight_sum);

    mpz_t sum;
    mpz_init_set_ui(sum, 0);
    for (unsigned int i = 0; i < n; i++) {
        mpz_add(sum, sum, weights[i]);
    }

    if (mpz_cmp_ui(sum, 0) == 0) {
        for (unsigned int i = 0; i < n; i++) {
            mpz_clear(out->prob.data[i]);
        }
        free(out->prob.data);
        vector_free(&out->aliases);
        vector_free(&out->small);
        vector_free(&out->large);
        mpz_clear(out->weight_sum);
        mpz_clear(sum);
        return WEIGHTED_ALL_ZERO;
    }

    mpz_set(out->weight_sum, sum);

    mpz_t p, tmp;
    mpz_inits(p, tmp, NULL);

    for (unsigned int i = 0; i < n; i++) {
        mpz_mul_ui(p, weights[i], n);
        mpz_set(out->prob.data[i], p);

        if (mpz_cmp(p, sum) < 0)
            vector_push(&out->small, i);
        else
            vector_push(&out->large, i);
    }

    while (out->small.size > 0 && out->large.size > 0) {
        unsigned int s = vector_get(&out->small, out->small.size - 1);
        out->small.size--;

        unsigned int l = vector_get(&out->large, out->large.size - 1);
        out->large.size--;

        vector_set(&out->aliases, s, l);

        mpz_set(tmp, out->prob.data[l]);
        mpz_add(tmp, tmp, out->prob.data[s]);
        mpz_sub(tmp, tmp, sum);
        mpz_set(out->prob.data[l], tmp);

        if (mpz_cmp(tmp, sum) < 0)
            vector_push(&out->small, l);
        else
            vector_push(&out->large, l);
    }

    while (out->large.size > 0) {
        unsigned int l = vector_get(&out->large, out->large.size - 1);
        out->large.size--;
        mpz_set(out->prob.data[l], sum);
    }

    while (out->small.size > 0) {
        unsigned int s = vector_get(&out->small, out->small.size - 1);
        out->small.size--;
        mpz_set(out->prob.data[s], sum);
    }

    out->n = n;

    mpz_clears(p, tmp, sum, NULL);
    return WEIGHTED_OK;
}


struct sample_alias_rust_gmp_s read_sample_alias_rust_gmp(char *fname) {
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
        perror(fname);
        exit(EXIT_FAILURE);
    }

    char *line = NULL;
    size_t len = 0;
    if (getline(&line, &len, fp) == -1) {
        fprintf(stderr, "read_sample_alias_rust_gmp: cannot read Z in %s\n", fname);
        fclose(fp);
        exit(EXIT_FAILURE);
    }
    free(line);

    int n = 0;
    if (fscanf(fp, "%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "read_sample_alias_rust_gmp: cannot read n in %s\n", fname);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    mpz_t *weights = malloc(n * sizeof(mpz_t));
    if (weights == NULL) {
        perror("malloc");
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    char buffer[2048];
    for (int i = 0; i < n; ++i) {
        mpz_init(weights[i]);
        if (fscanf(fp, "%2047s", buffer) != 1 || mpz_set_str(weights[i], buffer, 10) != 0) {
            fprintf(stderr, "read_sample_alias_rust_gmp: invalid weight '%s' at index %d\n", buffer, i);
            for (int j = 0; j <= i; ++j) mpz_clear(weights[j]);
            free(weights);
            fclose(fp);
            exit(EXIT_FAILURE);
        }
    }
    fclose(fp);

    struct sample_alias_rust_gmp_s sampler;
    WeightedError err = weighted_alias_new_gmp(
        (alias_rust_gmp_s *)&sampler, weights, (unsigned int)n);

    for (int i = 0; i < n; ++i) {
        mpz_clear(weights[i]);
    }
    free(weights);

    if (err != WEIGHTED_OK) {
        fprintf(stderr, "weighted_alias_new_gmp failed for %s (error %d)\n", fname, err);
        exit(EXIT_FAILURE);
    }

    return sampler;
}

void free_sample_alias_rust_gmp_s(struct sample_alias_rust_gmp_s x) {
    vector_free(&x.aliases);
    vector_free(&x.small);
    vector_free(&x.large);
    vector_mpz_free(&x.prob);
    mpz_clear(x.weight_sum);
}
