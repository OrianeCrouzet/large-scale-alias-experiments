/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
  This file is part of Alias_experiments, distributed under the GNU General
  Public License v3.0 or later (see LICENSE).
*/
/*
  Name:     alias_integers_v2.c
  Purpose:  Integer-Alias sampler, version 2: construction, I/O and size.
            The table is built by preprocess_alias_integers_opt (same cells as
            Algorithm 1), then each cell gets the 32 first digits of T_i / cs.
*/

#include <stdio.h>
#include <stdlib.h>

#include "alias_integers_v2.h"


static struct sample_alias_integers_v2_s from_opt(struct sample_alias_integers_opt_s base) {
    struct sample_alias_integers_v2_s s;
    s.base = base;
    s.q = s.base.q;
    s.q_reject = (uint32_t)(-s.q) % s.q;
    s.virtual_obj = s.base.virtual_obj;
    s.cells = malloc(s.q * sizeof(struct alias_v2_cell));

    mpz_t t;
    mpz_init(t);
    for (uint32_t i = 0; i < s.q; i++) {
        const struct alias_opt_cell *c = &s.base.cells[i];
        struct alias_v2_cell *v = &s.cells[i];
        v->primary = c->primary;
        v->thr = c->thr;
        if (c->alias == ALIAS_OPT_NONE) {
            v->alias = c->primary;
            v->d32 = UINT32_MAX;
            continue;
        }
        v->alias = c->alias;
        uint32_t end = (i + 1 < s.base.m) ? s.base.cells[i + 1].thr : s.base.thr_size;
        mpz_import(t, end - c->thr, -1, 1, 0, 0, s.base.Threshold + c->thr);
        mpz_mul_2exp(t, t, 32);
        mpz_fdiv_q(t, t, s.base.cs);        // < 2^32 since T_i < cs
        v->d32 = (uint32_t)mpz_get_ui(t);
    }
    mpz_clear(t);

    // The cells of alias.integers_opt are no longer needed.
    free(s.base.cells);
    s.base.cells = NULL;
    return s;
}

struct sample_alias_integers_v2_s read_sample_alias_integers_v2(char *fname) {
    return from_opt(read_sample_alias_integers_opt(fname));
}

struct sample_alias_integers_v2_s preprocess_alias_integers_v2(mpz_t *D, uint32_t n) {
    return from_opt(preprocess_alias_integers_opt(D, n));
}

void free_sample_alias_integers_v2(struct sample_alias_integers_v2_s x) {
    free(x.cells);
    free(x.base.Threshold);
    mpz_clear(x.base.cs);
}

// Same accounting as get_size_sample_alias_integers_opt_s: structures + useful bytes.
size_t get_size_sample_alias_integers_v2_s(struct sample_alias_integers_v2_s *s) {
    size_t total = sizeof(*s);
    total += (size_t)s->q * sizeof(struct alias_v2_cell);
    total += s->base.thr_size;
    total += (mpz_sizeinbase(s->base.cs, 2) + 7) / 8;
    return total;
}
