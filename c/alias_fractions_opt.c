/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
  This file is part of Alias_experiments, distributed under the GNU General
  Public License v3.0 or later (see LICENSE).
*/
/*
  Name:     alias_fractions_opt.c
  Purpose:  Optimized Fraction-Alias sampler (see alias_fractions_opt.h).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alias_fractions_opt.h"
#include "alias_fractions_gmp.h"
#include "flip_opt.h"

static void append(struct sample_alias_fractions_opt_s *s, size_t *cap, const void *src, size_t len) {
    if ((size_t)s->frac_size + len > UINT32_MAX) {
        fprintf(stderr, "[Fatal error] fractions too large.\n");
        exit(EXIT_FAILURE);
    }
    if (s->frac_size + len > *cap) {
        *cap = 2 * (s->frac_size + len);
        s->Fractions = realloc(s->Fractions, *cap);
    }
    memcpy(s->Fractions + s->frac_size, src, len);
    s->frac_size += len;
}

static void append_mpz(struct sample_alias_fractions_opt_s *s, size_t *cap, const mpz_t z) {
    size_t len = (mpz_sizeinbase(z, 2) + 7) / 8;
    unsigned char *buf = malloc(len);
    mpz_export(buf, &len, -1, 1, 0, 0, z);
    append(s, cap, buf, len);
    free(buf);
}

// The table is built by the original constructor, then stored compactly.
static struct sample_alias_fractions_opt_s from_gmp(struct sample_alias_fractions_gmp_s g) {
    struct sample_alias_fractions_opt_s s;
    s.taille = (uint32_t)g.taille;
    s.nbits = fo_nbits(s.taille);
    s.m = 0;
    s.frac_size = 0;
    s.cells = malloc(s.taille * sizeof(struct frac_opt_cell));
    size_t cap = 1024;
    s.Fractions = malloc(cap);
    uint32_t back = s.taille;

    mpz_t tmp, rem;
    mpz_inits(tmp, rem, NULL);
    for (uint32_t t = 0; t < s.taille; t++) {
        struct AliasEntryGMP *e = &g.table[t];
        mpz_srcptr num = mpq_numref(e->prob), den = mpq_denref(e->prob);
        struct frac_opt_cell *c;
        if (e->j < 0 || mpz_cmp(num, den) >= 0) {           // Bernoulli = 1 without bits
            c = &s.cells[--back];
            c->i = (uint32_t)e->i;
            c->j = FRAC_OPT_NONE;
        } else if (mpz_sgn(num) == 0) {                     // Bernoulli = 0 without bits
            c = &s.cells[--back];
            c->i = (uint32_t)e->j;
            c->j = FRAC_OPT_NONE;
        } else {
            c = &s.cells[s.m++];
            c->i = (uint32_t)e->i;
            c->j = (uint32_t)e->j;
            c->bits = fo_digits(num, den, tmp, rem);
            c->frac = s.frac_size;
            uint32_t num_len = (mpz_sizeinbase(num, 2) + 7) / 8;
            append(&s, &cap, &num_len, sizeof(num_len));
            append_mpz(&s, &cap, num);
            append_mpz(&s, &cap, den);
            continue;
        }
        c->bits = 0;
        c->frac = 0;
    }
    mpz_clears(tmp, rem, NULL);
    s.Fractions = realloc(s.Fractions, s.frac_size + 1);
    free_sample_alias_fractions_gmp_s(g);
    return s;
}

struct sample_alias_fractions_opt_s read_sample_alias_fractions_opt(char *fname) {
    return from_gmp(read_sample_alias_fractions_gmp(fname));
}

struct sample_alias_fractions_opt_s preprocess_alias_fractions_opt(mpz_t *w, uint32_t n) {
    return from_gmp(preprocess_alias_fractions_gmp(w, (int)n));
}

void free_sample_alias_fractions_opt(struct sample_alias_fractions_opt_s x) {
    free(x.cells);
    free(x.Fractions);
}

size_t get_size_sample_alias_fractions_opt_s(struct sample_alias_fractions_opt_s *s) {
    return sizeof(*s) + (size_t)s->taille * sizeof(struct frac_opt_cell) + s->frac_size;
}

static __attribute__((noinline)) uint32_t bernoulli_tail(const struct sample_alias_fractions_opt_s *x,
                                                         uint32_t t) {
    uint32_t start = x->cells[t].frac;
    uint32_t end = (t + 1 < x->m) ? x->cells[t + 1].frac : x->frac_size;
    uint32_t num_len;
    memcpy(&num_len, x->Fractions + start, sizeof(num_len));
    start += sizeof(num_len);
    mpz_t num, den;
    mpz_inits(num, den, NULL);
    mpz_import(num, num_len, -1, 1, 0, 0, x->Fractions + start);
    mpz_import(den, end - start - num_len, -1, 1, 0, 0, x->Fractions + start + num_len);
    uint32_t y = fo_bernoulli_tail(num, den);
    mpz_clears(num, den, NULL);
    return y;
}

uint32_t sample_alias_fractions_opt(struct sample_alias_fractions_opt_s *x) {
    uint32_t t = fo_uniform(x->taille, x->nbits);
    const struct frac_opt_cell *c = &x->cells[t];
    if (c->j == FRAC_OPT_NONE) {
        return c->i;
    }
    uint32_t b = fo_bernoulli_digits(c->bits);
    if (__builtin_expect(b == 2, 0)) {
        b = bernoulli_tail(x, t);
    }
    return b ? c->i : c->j;
}
