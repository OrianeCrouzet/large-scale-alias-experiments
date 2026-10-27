/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
  This file is part of Alias_experiments, distributed under the GNU General
  Public License v3.0 or later (see LICENSE).
*/
/*
  Name:     alias_integers_v2_sample.c
  Purpose:  Sampling step of the Integer-Alias sampler, version 2 (see
            alias_integers_v2.h): one 64-bit word per attempt, Lemire's uniform
            on the high half, comparison of the low half with 32 digits.
*/

#include "alias_integers_v2.h"
#include "flip_opt.h"
#include "counters.h"


// Next 64-bit word of the stream. Fast path when the stream is aligned on a word,
// which is the case except after a slow path that read single bits.
static inline uint64_t next_word(void) {
    uint32_t avail;
    uint64_t w = fo_peek(&avail);
    if (__builtin_expect(avail == 64, 1)) {
        fo_skip(64);
        return w;
    }
    return fo_bits(64);
}

// lo = D_i: V = (lo + U) / 2^32 with U uniform on [0, 1), and Pr(V < T_i / cs | lo = D_i)
// = frac(T_i 2^32 / cs) = (T_i 2^32 mod cs) / cs, drawn with Lumbroso's Bernoulli.
static __attribute__((noinline)) uint32_t bernoulli_tie(const struct sample_alias_integers_v2_s *x,
                                                        uint32_t i) {
    const struct alias_v2_cell *c = &x->cells[i];
    if (c->alias == c->primary) {
        return 1;
    }
    uint32_t end = (i + 1 < x->base.m) ? x->cells[i + 1].thr : x->base.thr_size;
    mpz_t num;
    mpz_init(num);
    mpz_import(num, end - c->thr, -1, 1, 0, 0, x->base.Threshold + c->thr);
    mpz_mul_2exp(num, num, 32);
    mpz_fdiv_r(num, num, x->base.cs);
    uint32_t y = 0;
    if (mpz_sgn(num) != 0) {
        for (;;) {
            mpz_mul_2exp(num, num, 1);
            int cmp = mpz_cmp(num, x->base.cs);
            if (cmp == 0) {
                y = fo_bit();
                break;
            }
            y = (cmp > 0);
            if (y) {
                mpz_sub(num, num, x->base.cs);
            }
            if (fo_bit()) {
                break;
            }
        }
    }
    mpz_clear(num);
    return y;
}

// Lemire's rejection (probability < q / 2^32): fresh 32-bit words until accepted.
static __attribute__((noinline)) uint32_t uniform_retry(const struct sample_alias_integers_v2_s *x,
                                                        uint64_t m) {
    while ((uint32_t)m < x->q_reject) {
        m = fo_bits(32) * (uint64_t)x->q;
    }
    return (uint32_t)(m >> 32);
}

uint32_t sample_alias_integers_v2(struct sample_alias_integers_v2_s *x) {
    uint32_t result;

    REJET -= 1;

    do {
        uint64_t w = next_word();
        uint64_t m = (w >> 32) * (uint64_t)x->q;
        uint32_t i = (uint32_t)(m >> 32);
        if (__builtin_expect((uint32_t)m < x->q, 0)) {
            i = uniform_retry(x, m);
        }
        const struct alias_v2_cell *c = &x->cells[i];
        uint32_t lo = (uint32_t)w;
        uint32_t b = lo < c->d32;
        if (__builtin_expect(lo == c->d32, 0)) {
            b = bernoulli_tie(x, i);
        }
        result = b ? c->primary : c->alias;

        REJET += 1;
    } while (result == x->virtual_obj);

    return result;
}
