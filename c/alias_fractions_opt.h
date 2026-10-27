/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
  This file is part of Alias_experiments, distributed under the GNU General
  Public License v3.0 or later (see LICENSE).
*/
/*
  Name:     alias_fractions_opt.h
  Purpose:  Optimized Fraction-Alias sampler (alias.fractions_opt).
            Same table as alias.fractions_gmp (Vose construction on exact
            rational probabilities, prob = num/den as reduced GMP fractions)
            and same sampling (uniform cell, Lumbroso Bernoulli(num/den)).
            The fractions are kept exactly; the first 31 binary digits of
            each fraction are precomputed so that GMP is only used with
            probability 2^-31 (see fo_bernoulli_digits in flip_opt.h).
*/

#ifndef ALIAS_FRACTIONS_OPT_H
#define ALIAS_FRACTIONS_OPT_H

#include <stddef.h>
#include <stdint.h>
#include <gmp.h>

#define FRAC_OPT_NONE UINT32_MAX

struct frac_opt_cell {
    uint32_t bits;      // digits of num/den (fo_digits)
    uint32_t i;
    uint32_t j;         // FRAC_OPT_NONE: the cell always yields i (no bits)
    uint32_t frac;      // offset in Fractions: [num length (4 bytes)][num bytes][den bytes]
};

struct sample_alias_fractions_opt_s {
    uint32_t taille;
    uint32_t nbits;
    uint32_t m;                 // cells 0..m-1 have an alias
    uint32_t frac_size;
    struct frac_opt_cell *cells;
    unsigned char *Fractions;   // exact fractions of cells 0..m-1 (slow path only)
};

struct sample_alias_fractions_opt_s read_sample_alias_fractions_opt(char *fname);
// Construction from the weights w_1..w_n (w is not modified).
struct sample_alias_fractions_opt_s preprocess_alias_fractions_opt(mpz_t *w, uint32_t n);
void free_sample_alias_fractions_opt(struct sample_alias_fractions_opt_s x);
size_t get_size_sample_alias_fractions_opt_s(struct sample_alias_fractions_opt_s *s);
uint32_t sample_alias_fractions_opt(struct sample_alias_fractions_opt_s *x);

#endif
