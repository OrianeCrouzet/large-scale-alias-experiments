/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
*/

/*
  check_samplers -- goodness-of-fit test of the six samplers.

    check_samplers FILE.dist [SAMPLES]

  For each sampler, draws SAMPLES objects (default 10^6) and compares the counts
  with the exact probabilities w_i / sum(w) by a chi-square test (objects with an
  expected count below 5 are pooled). The p-value is approximated with the
  Wilson-Hilferty transform; the test fails for z > 4.5 (p < 4e-6) or if an
  index >= n is returned. Exit status 0 if every sampler passes.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "aldr_opt.h"
#include "alias_fractions_opt.h"
#include "alias_integers_opt.h"
#include "alias_integers_v2.h"
#include "alias_rust_opt.h"
#include "fldr_opt.h"
#include "flip_opt.h"

unsigned long NUM_RNG_CALLS = 0;
unsigned long REJET = 0;

static uint32_t n_obj;
static double *prob;
static uint64_t *count;

/* Exact probabilities, in double precision. */
static void read_probabilities(const char *path) {
    FILE *fp = fopen(path, "r");
    if (fp == NULL) {
        perror(path);
        exit(EXIT_FAILURE);
    }
    mpz_t z, total;
    mpz_inits(z, total, NULL);
    int n = 0;
    if (mpz_inp_str(z, fp, 10) == 0 || fscanf(fp, "%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "[Error] invalid Z or n in %s\n", path);
        exit(EXIT_FAILURE);
    }
    n_obj = (uint32_t)n;
    mpz_t *w = malloc((size_t)n * sizeof(mpz_t));
    for (int i = 0; i < n; i++) {
        mpz_init(w[i]);
        if (mpz_inp_str(w[i], fp, 10) == 0) {
            fprintf(stderr, "[Error] cannot read weight %d in %s\n", i, path);
            exit(EXIT_FAILURE);
        }
        mpz_add(total, total, w[i]);
    }
    fclose(fp);
    prob = malloc((size_t)n * sizeof(double));
    long et;
    double mt = mpz_get_d_2exp(&et, total);
    for (int i = 0; i < n; i++) {
        long ei;
        double mi = mpz_get_d_2exp(&ei, w[i]);
        prob[i] = mpz_sgn(w[i]) ? ldexp(mi / mt, (int)(ei - et)) : 0.0;
        mpz_clear(w[i]);
    }
    free(w);
    mpz_clears(z, total, NULL);
}

/* Chi-square statistic and Wilson-Hilferty z-score of the counts. */
static int report(const char *name, uint64_t samples, int bad_index) {
    double chi2 = 0, pooled_obs = 0, pooled_exp = 0;
    unsigned classes = 0;
    for (uint32_t i = 0; i < n_obj; i++) {
        double e = prob[i] * (double)samples;
        if (e < 5) {
            pooled_obs += (double)count[i];
            pooled_exp += e;
            continue;
        }
        chi2 += ((double)count[i] - e) * ((double)count[i] - e) / e;
        classes++;
    }
    if (pooled_exp > 0) {
        chi2 += (pooled_obs - pooled_exp) * (pooled_obs - pooled_exp) / pooled_exp;
        classes++;
    }
    double df = classes > 1 ? classes - 1 : 1;
    double z = (cbrt(chi2 / df) - (1 - 2 / (9 * df))) / sqrt(2 / (9 * df));
    int ok = !bad_index && z < 4.5;
    printf("  %-20s chi2 %10.1f  df %7.0f  z %+6.2f  %s\n", name, chi2, df, z, ok ? "ok" : "FAIL");
    return ok;
}

/* first: index of the first object returned by the sampler (FLDR numbers the
 * objects from 1, as the leaf labels of Saad's encoding). */
#define CHECK(key, first, struct_name, func_read, func_sample, func_free)  \
    do {                                                                   \
        struct struct_name s = func_read(path);                            \
        memset(count, 0, n_obj * sizeof(uint64_t));                        \
        int bad = 0;                                                       \
        for (uint64_t k = 0; k < samples; k++) {                           \
            uint32_t x = (uint32_t)func_sample(&s) - (first);              \
            if (x >= n_obj) { bad = 1; break; }                            \
            count[x]++;                                                    \
        }                                                                  \
        func_free(s);                                                      \
        passed &= report(key, samples, bad);                               \
    } while (0)

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s FILE.dist [SAMPLES]\n", argv[0]);
        return EXIT_FAILURE;
    }
    char *path = argv[1];
    uint64_t samples = argc > 2 ? strtoull(argv[2], NULL, 10) : 1000000;
    read_probabilities(path);
    count = malloc(n_obj * sizeof(uint64_t));
    printf("%s (n = %u, %llu samples per sampler)\n", path, n_obj, (unsigned long long)samples);

    int passed = 1;
    fo_start();
    CHECK("alias.integers_opt", 0, sample_alias_integers_opt_s, read_sample_alias_integers_opt,
          sample_alias_integers_opt, free_sample_alias_integers_opt);
    CHECK("alias.integers_v2", 0, sample_alias_integers_v2_s, read_sample_alias_integers_v2,
          sample_alias_integers_v2, free_sample_alias_integers_v2);
    CHECK("rej.enc_opt", 1, sample_fldr_opt_s, read_sample_fldr_opt, sample_fldr_opt,
          free_sample_fldr_opt);
    CHECK("aldr_opt", 0, sample_aldr_opt_s, read_sample_aldr_opt, sample_aldr_opt,
          free_sample_aldr_opt);
    CHECK("alias.rust_opt", 0, sample_alias_rust_opt_s, read_sample_alias_rust_opt,
          sample_alias_rust_opt, free_sample_alias_rust_opt);
    CHECK("alias.fractions_opt", 0, sample_alias_fractions_opt_s, read_sample_alias_fractions_opt,
          sample_alias_fractions_opt, free_sample_alias_fractions_opt);
    fo_stop();

    free(count);
    free(prob);
    return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
