/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
  This file is part of Alias_experiments, distributed under the GNU General
  Public License v3.0 or later (see LICENSE).
*/
/*
  Name:     alias_integers_opt.c
  Purpose:  Optimized Integer-Alias sampler: construction, I/O and size.
            The table is the one of preprocess_gmp_alias_integers (Algorithm 1),
            stored as an array of cells with precomputed Bernoulli digits.
*/

#include <stdio.h>
#include <stdlib.h>

#include "alias_integers_opt.h"
#include "flip_opt.h"


// Cells with an alias are filled from the front (index s->m), the others from the back.
static void set_cell_alone(struct sample_alias_integers_opt_s *s, uint32_t *back, uint32_t x) {
    struct alias_opt_cell *c = &s->cells[--*back];
    c->primary = x;
    c->alias = ALIAS_OPT_NONE;
    c->bits = 0;
    c->thr = 0;
}

// Cell (x2, x) where x2 has weight t < cs.
static void set_cell(struct sample_alias_integers_opt_s *s, uint32_t *back, size_t *capacity,
                     const mpz_t t, uint32_t x2, uint32_t x, mpz_t tmp, mpz_t rem) {
    if (mpz_sgn(t) == 0) {
        // Bernoulli(0/cs) = 0 without consuming bits: the cell always yields x.
        set_cell_alone(s, back, x);
        return;
    }
    struct alias_opt_cell *c = &s->cells[s->m];
    c->primary = x2;
    c->alias = x;
    size_t len = (mpz_sizeinbase(t, 2) + 7) / 8;
    if ((size_t)s->thr_size + len > UINT32_MAX) {
        fprintf(stderr, "[Fatal error] thresholds too large.\n");
        exit(EXIT_FAILURE);
    }
    if (s->thr_size + len > *capacity) {
        *capacity = 2 * (s->thr_size + len);
        s->Threshold = realloc(s->Threshold, *capacity);
    }
    mpz_export(s->Threshold + s->thr_size, &len, -1, 1, 0, 0, t);
    c->thr = s->thr_size;
    s->thr_size += len;
    c->bits = fo_digits(t, s->cs, tmp, rem);
    s->m++;
}

// D: n weights, with room for one more (virtual object). D is modified.
struct sample_alias_integers_opt_s preprocess_alias_integers_opt(mpz_t *D, uint32_t n) {
    struct sample_alias_integers_opt_s s;
    mpz_init(s.cs);

    mpz_t w, r, tmp, rem;
    mpz_inits(w, r, tmp, rem, NULL);

    for (uint32_t i = 0; i < n; i++) {
        mpz_add(w, w, D[i]);
    }
    mpz_fdiv_q_ui(s.cs, w, n);                  // cs = w / n
    if (mpz_sgn(s.cs) == 0) {
        fprintf(stderr, "[Fatal error] cs is 0 (total weight w < n).\n");
        exit(EXIT_FAILURE);
    }
    mpz_fdiv_qr(tmp, r, w, s.cs);               // q = w / cs, r = w % cs

    s.virtual_obj = ALIAS_OPT_NONE;
    if (mpz_sgn(r) > 0) {
        s.virtual_obj = n;
        mpz_init(D[n]);
        mpz_sub(D[n], s.cs, r);                 // virtual object of weight cs - r
        mpz_add_ui(tmp, tmp, 1);
        n++;
    }
    if (!mpz_fits_uint_p(tmp) || mpz_get_ui(tmp) >= ALIAS_OPT_NONE) {
        fprintf(stderr, "[Fatal error] alias table too large.\n");
        exit(EXIT_FAILURE);
    }
    s.q = (uint32_t)mpz_get_ui(tmp);
    s.q_nbits = (s.q > 1) ? 32 - __builtin_clz(s.q - 1) : 0;

    s.cells = malloc(s.q * sizeof(struct alias_opt_cell));
    s.m = 0;
    s.thr_size = 0;
    size_t capacity = (size_t)n * ((mpz_sizeinbase(s.cs, 2) + 7) / 8);
    s.Threshold = malloc(capacity + 1);
    uint32_t back = s.q;

    // Vose-like construction with stacks L (> cs) and S (<= cs), same order as cons_alias_gmp.
    uint32_t *small = malloc(n * sizeof(uint32_t));
    uint32_t *large = malloc(n * sizeof(uint32_t));
    uint32_t n_small = 0, n_large = 0;
    for (uint32_t i = 0; i < n; i++) {
        if (mpz_cmp(D[i], s.cs) > 0) {
            large[n_large++] = i;
        } else {
            small[n_small++] = i;
        }
    }

    while (n_large > 0) {
        uint32_t x = large[--n_large];
        if (n_small > 0) {
            uint32_t x2 = small[--n_small];
            if (mpz_cmp(D[x2], s.cs) != 0) {
                set_cell(&s, &back, &capacity, D[x2], x2, x, tmp, rem);
                mpz_sub(D[x], D[x], s.cs);      // D[x] -= cs - D[x2]
                mpz_add(D[x], D[x], D[x2]);
            } else {
                set_cell_alone(&s, &back, x2);
            }
        } else {
            set_cell_alone(&s, &back, x);
            mpz_sub(D[x], D[x], s.cs);
        }
        if (mpz_cmp(D[x], s.cs) > 0) {
            large[n_large++] = x;
        } else {
            small[n_small++] = x;
        }
    }
    while (n_small > 0) {
        set_cell_alone(&s, &back, small[--n_small]);
    }

    if (back != s.m) {
        fprintf(stderr, "[Fatal error] %u cells built, %u expected.\n",
                s.m + s.q - back, s.q);
        exit(EXIT_FAILURE);
    }
    s.Threshold = realloc(s.Threshold, s.thr_size + 1);

    free(small);
    free(large);
    mpz_clears(w, r, tmp, rem, NULL);
    return s;
}

struct sample_alias_integers_opt_s read_sample_alias_integers_opt(char *fname) {
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
        perror(fname);
        exit(EXIT_FAILURE);
    }

    // Line 1: Z (unused); line 2: n w_1 ... w_n; line 3: entropy (unused).
    mpz_t Z;
    mpz_init(Z);
    int n = 0;
    if (mpz_inp_str(Z, fp, 10) == 0 || fscanf(fp, "%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "[Error] invalid Z or n in %s\n", fname);
        exit(EXIT_FAILURE);
    }

    mpz_t *D = malloc((n + 1) * sizeof(mpz_t));
    for (int i = 0; i < n; i++) {
        mpz_init(D[i]);
        if (mpz_inp_str(D[i], fp, 10) == 0) {
            fprintf(stderr, "[Error] cannot read weight %d in %s\n", i, fname);
            exit(EXIT_FAILURE);
        }
    }
    fclose(fp);

    struct sample_alias_integers_opt_s s = preprocess_alias_integers_opt(D, n);

    int n_init = n + (s.virtual_obj != ALIAS_OPT_NONE);
    for (int i = 0; i < n_init; i++) {
        mpz_clear(D[i]);
    }
    free(D);
    mpz_clear(Z);
    return s;
}

void free_sample_alias_integers_opt(struct sample_alias_integers_opt_s x) {
    free(x.Threshold);
    free(x.cells);
    mpz_clear(x.cs);
}

// Same accounting as get_size_sample_gmp_alias_integers_s: structures + useful bytes.
size_t get_size_sample_alias_integers_opt_s(struct sample_alias_integers_opt_s *s) {
    size_t total = sizeof(*s);
    total += (size_t)s->q * sizeof(struct alias_opt_cell);
    total += s->thr_size;
    total += (mpz_sizeinbase(s->cs, 2) + 7) / 8;
    return total;
}
