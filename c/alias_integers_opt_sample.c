/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
  This file is part of Alias_experiments, distributed under the GNU General
  Public License v3.0 or later (see LICENSE).
*/
/*
  Name:     alias_integers_opt_sample.c
  Purpose:  Sampling step of the optimized Integer-Alias sampler.

  Same algorithm as sample_gmp_alias_integers (Lumbroso uniform + Bernoulli,
  virtual-object rejection), but:
    - the uniform over the q cells is done on machine integers;
    - the Bernoulli(Threshold/cs) uses the 31 precomputed binary digits of
      Threshold/cs: it reads coins until the first 1 (ctz on the current word)
      and returns the corresponding digit; GMP is only used when the first
      31 coins are all 0 (probability 2^-31);
    - bits come from the optimized stream of flip_opt.h.
*/

#include "alias_integers_opt.h"
#include "flip_opt.h"
#include "counters.h"


static __attribute__((noinline)) uint32_t bernoulli_tail(const struct sample_alias_integers_opt_s *x,
                                                         uint32_t i) {
    uint32_t start = x->cells[i].thr;
    uint32_t end = (i + 1 < x->m) ? x->cells[i + 1].thr : x->thr_size;
    mpz_t t;
    mpz_init(t);
    mpz_import(t, end - start, -1, 1, 0, 0, x->Threshold + start);
    uint32_t y = fo_bernoulli_tail(t, x->cs);
    mpz_clear(t);
    return y;
}

uint32_t sample_alias_integers_opt(struct sample_alias_integers_opt_s *x) {
    uint32_t result;

    REJET -= 1;

    do {
        uint32_t i = fo_uniform(x->q, x->q_nbits);
        const struct alias_opt_cell *c = &x->cells[i];
        uint32_t b = 1;
        if (c->alias != ALIAS_OPT_NONE) {
            b = fo_bernoulli_digits(c->bits);
            if (__builtin_expect(b == 2, 0)) {
                b = bernoulli_tail(x, i);
            }
        }
        result = b ? c->primary : c->alias;

        REJET += 1;
    } while (result == x->virtual_obj);

    return result;
}
