/*
  Derived from main.c and macros.c of the fast-loaded-dice-roller-experiments
  repository of F. A. Saad et al.
  (https://github.com/probsys/fast-loaded-dice-roller-experiments),
  Copyright (C) 2020 Feras A. Saad and C. E. Freer,
  released under the Apache License 2.0 (see LICENSES/Apache-2.0.txt).

  Modified in 2026 by Oriane Crouzet, Antoine Genitrini and Mehdi Naima:
  restricted to the six samplers of the paper, which read .dist files and use
  the random-bit stream of flip_opt.h.
  This file is distributed under the GNU General Public License v3.0 or later
  as part of Alias_experiments (see LICENSE).
  SPDX-License-Identifier: GPL-3.0-or-later
*/

/*
  measure_runtime -- sampling time of one sampler on one distribution.

    measure_runtime SEED STEPS SAMPLER FILE.dist

  Builds the sampler from FILE.dist (not timed), then draws STEPS objects and
  measures the time of the sampling loop with clock() (CPU time of the
  process, including the thread that produces the random bits).

  Output: a first line repeating the arguments, then
    SAMPLER SECONDS WORDS REJECTIONS
  where WORDS is the number of 32-bit random words consumed and REJECTIONS the
  number of virtual-object rejections (Integer-Alias only).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "aldr_opt.h"
#include "alias_fractions_opt.h"
#include "alias_integers_opt.h"
#include "alias_integers_v2.h"
#include "alias_rust_opt.h"
#include "fldr_opt.h"
#include "flip_opt.h"

unsigned long NUM_RNG_CALLS = 0;
unsigned long REJET = 0;

#define CACHE_SIZE_BYTES (32 * 1024 * 1024 * 4)    /* 128 MB */

/* Writes a buffer larger than the caches, so that every run starts cold. */
static void flush_cache(void) {
    volatile char *dummy = malloc(CACHE_SIZE_BYTES);
    if (!dummy) {
        return;
    }
    for (size_t i = 0; i < CACHE_SIZE_BYTES; i++) {
        dummy[i] = i % 256;
    }
    free((void *)dummy);
}

/* Reads and builds the sampler (not timed), then times the sampling loop. */
#define READ_SAMPLE_TIME(key, struct_name, func_read, func_sample, func_free) \
    if (strcmp(sampler, key) == 0) {                                          \
        struct struct_name s = func_read(path);                               \
        t = clock();                                                          \
        for (int i = 0; i < steps; i++) {                                     \
            x += func_sample(&s);                                             \
        }                                                                     \
        t = clock() - t;                                                      \
        func_free(s);                                                         \
    }

int main(int argc, char **argv) {
    if (argc != 5) {
        printf("usage: %s seed steps sampler path\n", argv[0]);
        return EXIT_FAILURE;
    }
    int seed = atoi(argv[1]);
    int steps = atoi(argv[2]);
    char *sampler = argv[3];
    char *path = argv[4];

    printf("%d %d %s %s\n", seed, steps, sampler, path);
    srand(seed);
    flush_cache();
    fo_start();

    long x = 0;
    clock_t t = 0;
    READ_SAMPLE_TIME("alias.integers_opt", sample_alias_integers_opt_s,
                     read_sample_alias_integers_opt, sample_alias_integers_opt,
                     free_sample_alias_integers_opt)
    else READ_SAMPLE_TIME("alias.integers_v2", sample_alias_integers_v2_s,
                          read_sample_alias_integers_v2, sample_alias_integers_v2,
                          free_sample_alias_integers_v2)
    else READ_SAMPLE_TIME("rej.enc_opt", sample_fldr_opt_s,
                          read_sample_fldr_opt, sample_fldr_opt, free_sample_fldr_opt)
    else READ_SAMPLE_TIME("aldr_opt", sample_aldr_opt_s,
                          read_sample_aldr_opt, sample_aldr_opt, free_sample_aldr_opt)
    else READ_SAMPLE_TIME("alias.rust_opt", sample_alias_rust_opt_s,
                          read_sample_alias_rust_opt, sample_alias_rust_opt,
                          free_sample_alias_rust_opt)
    else READ_SAMPLE_TIME("alias.fractions_opt", sample_alias_fractions_opt_s,
                          read_sample_alias_fractions_opt, sample_alias_fractions_opt,
                          free_sample_alias_fractions_opt)
    else {
        printf("Unknown sampler: %s\n", sampler);
        fo_stop();
        return EXIT_FAILURE;
    }

    double e = ((double)t) / CLOCKS_PER_SEC;
    printf("%s %1.5f %ld %ld\n", sampler, e, NUM_RNG_CALLS, REJET);
    fo_stop();
    return (x == -1) ? EXIT_FAILURE : EXIT_SUCCESS;    /* x keeps the loop alive */
}
