/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
  This file is part of Alias_experiments, distributed under the GNU General
  Public License v3.0 or later (see LICENSE).
*/
/*
  Name:     alias_integers_opt.h
  Purpose:  Optimized Integer-Alias sampler (arbitrary-precision weights).
            Duplicate of alias.integers_gmp, with the same cells (Algorithm 1)
            and the same random-bit consumption, but a GMP-free sampling step
            in all but a 2^-31 fraction of the Bernoulli trials.
            Random bits come from the optimized stream of flip_opt.h.
            Cells with an alias are stored first (the order of the cells is
            irrelevant since the cell is drawn uniformly), so that the exact
            thresholds are only stored for them.
*/

#ifndef ALIAS_INTEGERS_OPT_H
#define ALIAS_INTEGERS_OPT_H

#include <stddef.h>
#include <stdint.h>
#include <gmp.h>

#define ALIAS_OPT_NONE UINT32_MAX

// One cell of the alias table (16 bytes).
// bits: digits of Threshold/cs packed by fo_digits (flip_opt.h).
// alias == ALIAS_OPT_NONE: the cell holds `primary` alone (Bernoulli = 1, no bits).
// thr: offset in Threshold of the exact threshold (cells with an alias only).
struct alias_opt_cell {
    uint32_t bits;
    uint32_t primary;
    uint32_t alias;
    uint32_t thr;
};

struct sample_alias_integers_opt_s {
    struct alias_opt_cell *cells;
    uint32_t q;            // number of cells
    uint32_t q_nbits;      // bit length of q-1 (first block of the FDR uniform)
    uint32_t virtual_obj;  // index of the virtual object, ALIAS_OPT_NONE if absent
    uint32_t m;            // number of cells with an alias (cells 0..m-1)
    uint32_t thr_size;     // bytes used in Threshold
    unsigned char *Threshold; // exact thresholds of cells 0..m-1, little-endian bytes (slow path only)
    mpz_t cs;
};

struct sample_alias_integers_opt_s preprocess_alias_integers_opt(mpz_t *D, uint32_t n);
struct sample_alias_integers_opt_s read_sample_alias_integers_opt(char *fname);
void free_sample_alias_integers_opt(struct sample_alias_integers_opt_s x);
size_t get_size_sample_alias_integers_opt_s(struct sample_alias_integers_opt_s *s);

uint32_t sample_alias_integers_opt(struct sample_alias_integers_opt_s *x);

#endif
