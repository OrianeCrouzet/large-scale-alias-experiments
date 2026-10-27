/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
*/

/*
  measure_preprocess -- preprocessing time of one sampler on one distribution.

    measure_preprocess SAMPLER FILE.dist

  The weights are read first, outside the timing; then the construction of the
  sampling structure alone is repeated (at least 5 times and 50 ms) and the
  smallest time is printed, in seconds.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#include "alias_integers_opt.h"
#include "alias_integers_v2.h"
#include "fldr_opt.h"
#include "aldr_opt.h"
#include "alias_rust_opt.h"
#include "alias_fractions_opt.h"

unsigned long NUM_RNG_CALLS = 0;
unsigned long REJET = 0;

#define MIN_TOTAL_SECONDS 0.05
#define MIN_REPEATS 5
#define MAX_REPEATS 1000

static mpz_t *read_weights(const char *fname, uint32_t *n_out, mpz_t Z) {
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
        perror(fname);
        exit(EXIT_FAILURE);
    }
    int n = 0;
    if (mpz_inp_str(Z, fp, 10) == 0 || fscanf(fp, "%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "[Error] invalid Z or n in %s\n", fname);
        exit(EXIT_FAILURE);
    }
    mpz_t *w = malloc((n + 1) * sizeof(mpz_t));
    for (int i = 0; i < n; i++) {
        mpz_init(w[i]);
        if (mpz_inp_str(w[i], fp, 10) == 0) {
            fprintf(stderr, "[Error] cannot read weight %d in %s\n", i, fname);
            exit(EXIT_FAILURE);
        }
    }
    fclose(fp);
    *n_out = (uint32_t)n;
    return w;
}

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

// Copy of the weights with room for the virtual object (Integer-Alias modifies D).
static mpz_t *copy_weights(mpz_t *w, uint32_t n) {
    mpz_t *D = malloc((n + 1) * sizeof(mpz_t));
    for (uint32_t i = 0; i < n; i++) {
        mpz_init_set(D[i], w[i]);
    }
    return D;
}

// Integer-Alias initializes D[n] only when there is a virtual object.
static void free_weights(mpz_t *D, uint32_t n, int with_virtual) {
    for (uint32_t i = 0; i < n + (with_virtual ? 1 : 0); i++) {
        mpz_clear(D[i]);
    }
    free(D);
}

// Time of one construction: copy (Integer-Alias only) outside the timing.
static double time_once(const char *sampler, mpz_t *w, uint32_t n, const mpz_t Z) {
    double t0, t1;
    if (strcmp(sampler, "alias.integers_opt") == 0) {
        mpz_t *D = copy_weights(w, n);
        t0 = now();
        struct sample_alias_integers_opt_s s = preprocess_alias_integers_opt(D, n);
        t1 = now();
        free_weights(D, n, s.virtual_obj != ALIAS_OPT_NONE);
        free_sample_alias_integers_opt(s);
    } else if (strcmp(sampler, "alias.integers_v2") == 0) {
        mpz_t *D = copy_weights(w, n);
        t0 = now();
        struct sample_alias_integers_v2_s s = preprocess_alias_integers_v2(D, n);
        t1 = now();
        free_weights(D, n, s.virtual_obj != ALIAS_OPT_NONE);
        free_sample_alias_integers_v2(s);
    } else if (strcmp(sampler, "rej.enc_opt") == 0) {
        t0 = now();
        struct sample_fldr_opt_s s = preprocess_fldr_opt(w, n, Z);
        t1 = now();
        free_sample_fldr_opt(s);
    } else if (strcmp(sampler, "aldr_opt") == 0) {
        t0 = now();
        struct sample_aldr_opt_s s = preprocess_aldr_opt(w, (int)n);
        t1 = now();
        free_sample_aldr_opt(s);
    } else if (strcmp(sampler, "alias.rust_opt") == 0) {
        t0 = now();
        struct sample_alias_rust_opt_s s = preprocess_alias_rust_opt(w, n);
        t1 = now();
        free_sample_alias_rust_opt(s);
    } else if (strcmp(sampler, "alias.fractions_opt") == 0) {
        t0 = now();
        struct sample_alias_fractions_opt_s s = preprocess_alias_fractions_opt(w, n);
        t1 = now();
        free_sample_alias_fractions_opt(s);
    } else {
        fprintf(stderr, "[Error] unknown sampler: %s\n", sampler);
        exit(EXIT_FAILURE);
    }
    return t1 - t0;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: %s SAMPLER FILE.dist\n", argv[0]);
        return EXIT_FAILURE;
    }
    uint32_t n;
    mpz_t Z;
    mpz_init(Z);
    mpz_t *w = read_weights(argv[2], &n, Z);

    double best = 1e300, total = 0;
    int reps = 0;
    while (reps < MIN_REPEATS || (total < MIN_TOTAL_SECONDS && reps < MAX_REPEATS)) {
        double t = time_once(argv[1], w, n, Z);
        if (t < best) {
            best = t;
        }
        total += t;
        reps++;
    }
    printf("%.9f\n", best);

    for (uint32_t i = 0; i < n; i++) {
        mpz_clear(w[i]);
    }
    free(w);
    mpz_clear(Z);
    return EXIT_SUCCESS;
}
