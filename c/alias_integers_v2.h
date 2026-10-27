/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
  This file is part of Alias_experiments, distributed under the GNU General
  Public License v3.0 or later (see LICENSE).
*/
/*
  Name:     alias_integers_v2.h
  Purpose:  Integer-Alias sampler, version 2: same table as alias.integers_opt
            (Algorithm 1), but the uniform and the Bernoulli of Lumbroso are
            replaced by faster exact algorithms that consume more random bits.

            Each attempt reads one 64-bit word W of the stream (flip_opt.h):
              - the cell is i = floor(hi * q / 2^32), hi = W >> 32 (Lemire's
                multiply-shift, exact thanks to a rejection of probability
                < q / 2^32, which draws fresh 32-bit words);
              - the Bernoulli(T_i / cs) compares lo = W mod 2^32 with the first
                32 binary digits D_i = floor(T_i 2^32 / cs): lo < D_i gives the
                cell's object, lo > D_i its alias, and lo = D_i (probability
                2^-32) is decided exactly with GMP and fresh bits.
            About 64 bits per attempt, against 17.6 for alias.integers_opt.
*/

#ifndef ALIAS_INTEGERS_V2_H
#define ALIAS_INTEGERS_V2_H

#include <stddef.h>
#include <stdint.h>
#include <gmp.h>

#include "alias_integers_opt.h"

// One cell (16 bytes). For a cell without alias, alias == primary and d32 = UINT32_MAX.
struct alias_v2_cell {
    uint32_t d32;          // first 32 binary digits of T_i / cs
    uint32_t primary;
    uint32_t alias;
    uint32_t thr;          // offset of the exact threshold in base.Threshold (cells with an alias)
};

struct sample_alias_integers_v2_s {
    struct alias_v2_cell *cells;
    uint32_t q;            // number of cells
    uint32_t q_reject;     // (2^32 - q) mod q: Lemire rejects when the low part is below it
    uint32_t virtual_obj;  // index of the virtual object, ALIAS_OPT_NONE if absent
    // Table of alias.integers_opt; only Threshold, m, thr_size and cs are used (slow path).
    struct sample_alias_integers_opt_s base;
};

struct sample_alias_integers_v2_s read_sample_alias_integers_v2(char *fname);
// D: n weights, with room for one more (virtual object). D is modified, as in
// preprocess_alias_integers_opt.
struct sample_alias_integers_v2_s preprocess_alias_integers_v2(mpz_t *D, uint32_t n);
void free_sample_alias_integers_v2(struct sample_alias_integers_v2_s x);
size_t get_size_sample_alias_integers_v2_s(struct sample_alias_integers_v2_s *s);

uint32_t sample_alias_integers_v2(struct sample_alias_integers_v2_s *x);

#endif
