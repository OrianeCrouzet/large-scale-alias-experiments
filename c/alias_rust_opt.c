/*
  The table is the one of the Rust crate weighted_rand
  (https://github.com/ichi-h/weighted_rand), Copyright (c) 2020 Ippee,
  released under the MIT License (see LICENSES/MIT-weighted_rand.txt),
  built by alias_rust_gmp.c.

  Modified in 2026 by Oriane Crouzet, Antoine Genitrini and Mehdi Naima:
  compact storage of the table and sampling step on GMP integers.
  This file is distributed under the GNU General Public License v3.0 or later
  as part of Alias_experiments (see LICENSE).
  SPDX-License-Identifier: GPL-3.0-or-later
*/
/*
  Name:     alias_rust_opt.c
  Purpose:  Optimized Rust-Alias sampler (see alias_rust_opt.h).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alias_rust_opt.h"
#include "alias_rust_gmp.h"
#include "flip_opt.h"

// The table is built by the original constructor, then stored compactly.
static struct sample_alias_rust_opt_s from_gmp(struct sample_alias_rust_gmp_s g) {
    struct sample_alias_rust_opt_s s;
    s.n = g.n;
    s.n_nbits = fo_nbits(g.n);
    mpz_init_set(s.W, g.weight_sum);
    if (mpz_cmp_ui(s.W, 1) > 0) {
        mpz_t w1;
        mpz_init(w1);
        mpz_sub_ui(w1, s.W, 1);
        s.W_nbits = mpz_sizeinbase(w1, 2);
        mpz_clear(w1);
    } else {
        s.W_nbits = 0;
    }
    s.limbs = mpz_size(s.W);
    s.aliases = malloc(s.n * sizeof(uint32_t));
    s.prob = calloc((size_t)s.n * s.limbs, sizeof(mp_limb_t));
    for (uint32_t i = 0; i < s.n; i++) {
        s.aliases[i] = (uint32_t)g.aliases.data[i];
        // prob[i] <= W, hence at most `limbs` limbs
        mpz_export(s.prob + (size_t)i * s.limbs, NULL, -1, sizeof(mp_limb_t), 0, 0, g.prob.data[i]);
    }
    mpz_init2(s.x, s.W_nbits + 2 * GMP_NUMB_BITS);
    mpz_init2(s.bound, s.W_nbits + 2 * GMP_NUMB_BITS);
    free_sample_alias_rust_gmp_s(g);
    return s;
}

struct sample_alias_rust_opt_s read_sample_alias_rust_opt(char *fname) {
    return from_gmp(read_sample_alias_rust_gmp(fname));
}

struct sample_alias_rust_opt_s preprocess_alias_rust_opt(mpz_t *w, uint32_t n) {
    struct sample_alias_rust_gmp_s g;
    if (weighted_alias_new_gmp((alias_rust_gmp_s *)&g, w, n) != WEIGHTED_OK) {
        fprintf(stderr, "[Error] weighted_alias_new_gmp failed.\n");
        exit(EXIT_FAILURE);
    }
    return from_gmp(g);
}

void free_sample_alias_rust_opt(struct sample_alias_rust_opt_s x) {
    free(x.aliases);
    free(x.prob);
    mpz_clears(x.W, x.x, x.bound, NULL);
}

size_t get_size_sample_alias_rust_opt_s(struct sample_alias_rust_opt_s *s) {
    size_t total = sizeof(*s);
    total += (size_t)s->n * sizeof(uint32_t);
    total += (size_t)s->n * s->limbs * sizeof(mp_limb_t);
    total += 3 * (size_t)s->limbs * sizeof(mp_limb_t);     // W and the two scratch integers
    return total;
}

// Fast Dice Roller on [0, W) in arbitrary precision, same steps as uniform_with_gmp.
static inline void uniform_W(struct sample_alias_rust_opt_s *s) {
    uint32_t nb = s->W_nbits;
    if (nb == 0) {
        mpz_set_ui(s->x, 0);
        return;
    }
    // x = nb random bits, written limb by limb
    uint32_t L = (nb + 63) / 64;
    mp_limb_t *p = mpz_limbs_write(s->x, L);
    for (uint32_t k = 0; k + 1 < L; k++) {
        p[k] = fo_bits(64);
    }
    p[L - 1] = fo_bits(nb - 64 * (L - 1));
    mpz_limbs_finish(s->x, L);
    mpz_set_ui(s->bound, 0);
    mpz_setbit(s->bound, nb);
    for (;;) {
        if (mpz_cmp(s->bound, s->W) >= 0) {
            if (mpz_cmp(s->x, s->W) < 0) {
                return;
            }
            mpz_sub(s->bound, s->bound, s->W);
            mpz_sub(s->x, s->x, s->W);
        }
        mpz_mul_2exp(s->bound, s->bound, 1);
        mpz_mul_2exp(s->x, s->x, 1);
        if (fo_bit()) {
            mpz_add_ui(s->x, s->x, 1);
        }
    }
}

uint32_t sample_alias_rust_opt(struct sample_alias_rust_opt_s *x) {
    uint32_t i = fo_uniform(x->n, x->n_nbits);
    uniform_W(x);
    mpz_t p;
    mpz_roinit_n(p, x->prob + (size_t)i * x->limbs, x->limbs);
    return (mpz_cmp(x->x, p) < 0) ? i : x->aliases[i];
}
