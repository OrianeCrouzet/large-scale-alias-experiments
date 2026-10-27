/*
  Derived from aldr.c of the amplified-loaded-dice-roller repository
  (https://github.com/probsys/amplified-loaded-dice-roller),
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab,
  released under the Apache License 2.0 (see LICENSES/Apache-2.0.txt).

  Modified in 2026 by Oriane Crouzet, Antoine Genitrini and Mehdi Naima:
  ALDR on arbitrary-precision weights (GMP), 32-bit storage of the tree,
  walk inside the current word of the random-bit stream.
  This file is distributed under the GNU General Public License v3.0 or later
  as part of Alias_experiments (see LICENSE).
  SPDX-License-Identifier: GPL-3.0-or-later
*/
/*
  Name:     aldr_opt.c
  Purpose:  Optimized ALDR sampler (see aldr_opt.h).
*/

#include <stdio.h>
#include <stdlib.h>

#include "aldr_opt.h"
#include "flip_opt.h"

// Same tree as preprocess_aldr_flat_k_gmp(a, n, 2), with c * a[i] computed once.
struct sample_aldr_opt_s preprocess_aldr_opt(mpz_t *a, int n) {
    mpz_t m, c, r, pow2K;
    mpz_inits(m, c, r, pow2K, NULL);
    for (int i = 0; i < n; i++) {
        mpz_add(m, m, a[i]);
    }
    // k = ceil(log2(m)), K = 2k
    size_t k = mpz_sizeinbase(m, 2);
    mpz_sub_ui(r, m, 1);
    if (mpz_sizeinbase(r, 2) < k) {             // same test as preprocess_aldr_flat_k_gmp
        k--;
    }
    size_t K = 2 * k;
    mpz_setbit(pow2K, K);
    mpz_fdiv_qr(c, r, pow2K, m);                // c = 2^K / m, r = 2^K % m

    mpz_t *Q = malloc(n * sizeof(mpz_t));
    size_t num_leaves = mpz_popcount(r);
    for (int i = 0; i < n; i++) {
        mpz_init(Q[i]);
        mpz_mul(Q[i], c, a[i]);
        num_leaves += mpz_popcount(Q[i]);
    }
    if (num_leaves >= UINT32_MAX || K + 1 >= UINT32_MAX) {
        fprintf(stderr, "[Fatal error] ALDR tree too large.\n");
        exit(EXIT_FAILURE);
    }

    struct sample_aldr_opt_s s;
    s.length_breadths = (uint32_t)(K + 1);
    s.length_leaves = (uint32_t)num_leaves;
    s.breadths = calloc(K + 1, sizeof(uint32_t));
    s.leaves_flat = malloc((num_leaves + 1) * sizeof(uint32_t));
    uint32_t location = 0;
    for (size_t j = 0; j <= K; j++) {
        size_t bit = K - j;
        if (mpz_tstbit(r, bit)) {
            s.leaves_flat[location++] = 0;
            s.breadths[j]++;
        }
        for (int i = 0; i < n; i++) {
            if (mpz_tstbit(Q[i], bit)) {
                s.leaves_flat[location++] = (uint32_t)i + 1;
                s.breadths[j]++;
            }
        }
    }

    for (int i = 0; i < n; i++) {
        mpz_clear(Q[i]);
    }
    free(Q);
    mpz_clears(m, c, r, pow2K, NULL);
    return s;
}

// Same file format as aldr_gmp: Z, then n w_1 ... w_n, then the entropy.
struct sample_aldr_opt_s read_sample_aldr_opt(char *fname) {
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
        perror(fname);
        exit(EXIT_FAILURE);
    }
    mpz_t Z;
    mpz_init(Z);
    int n = 0;
    if (mpz_inp_str(Z, fp, 10) == 0 || fscanf(fp, "%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "[Error] invalid Z or n in %s\n", fname);
        exit(EXIT_FAILURE);
    }
    mpz_t *a = malloc(n * sizeof(mpz_t));
    for (int i = 0; i < n; i++) {
        mpz_init(a[i]);
        if (mpz_inp_str(a[i], fp, 10) == 0) {
            fprintf(stderr, "[Error] cannot read weight %d in %s\n", i, fname);
            exit(EXIT_FAILURE);
        }
    }
    fclose(fp);
    struct sample_aldr_opt_s s = preprocess_aldr_opt(a, n);
    for (int i = 0; i < n; i++) {
        mpz_clear(a[i]);
    }
    free(a);
    mpz_clear(Z);
    return s;
}

void free_sample_aldr_opt(struct sample_aldr_opt_s x) {
    free(x.breadths);
    free(x.leaves_flat);
}

size_t get_size_sample_aldr_opt_s(struct sample_aldr_opt_s *s) {
    return sizeof(*s) + ((size_t)s->length_breadths + s->length_leaves) * sizeof(uint32_t);
}

// Same walk as sample_aldr_flat_gmp:
// at depth d, if val < breadths[d] the leaf leaves_flat[location + val] is reached
// (0 = reject, restart from the root); otherwise location += breadths[d] and
// val = ((val - breadths[d]) << 1) | bit.
int sample_aldr_opt(struct sample_aldr_opt_s *x) {
    const uint32_t *breadths = x->breadths;
    const uint32_t *leaves = x->leaves_flat;
    uint32_t depth = 0;
    uint64_t location = 0, val = 0;
    uint32_t avail = 0, used = 0;
    uint64_t w = 0;
    for (;;) {
        uint64_t b = breadths[depth];
        if (val < b) {
            uint32_t leaf = leaves[location + val];
            if (leaf) {
                fo_skip(used);
                return (int)leaf - 1;
            }
            depth = 0;
            location = 0;
            val = 0;
            continue;
        }
        if (used == avail) {
            fo_skip(used);
            w = fo_peek(&avail);
            used = 0;
        }
        location += b;
        val = ((val - b) << 1) | (w & 1);
        w >>= 1;
        used++;
        depth++;
    }
}
