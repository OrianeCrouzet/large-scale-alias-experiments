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
  Name:     aldr_opt.h
  Purpose:  Optimized ALDR sampler (aldr_opt).
            Same flattened ALDR tree as aldr_gmp (built from arbitrary-precision
            weights, depth K = 2k), but the breadths and leaves are stored as
            32-bit integers (the number of nodes of a level is at most n + 1),
            and the tree is walked inside the current word of the random-bit
            stream (flip_opt.h).
*/

#ifndef ALDR_OPT_H
#define ALDR_OPT_H

#include <stddef.h>
#include <stdint.h>
#include <gmp.h>

struct sample_aldr_opt_s {
    uint32_t length_breadths;  // K + 1
    uint32_t length_leaves;
    uint32_t *breadths;        // number of leaves at each depth
    uint32_t *leaves_flat;     // leaf labels, 0 = reject, i + 1 = object i
};

struct sample_aldr_opt_s preprocess_aldr_opt(mpz_t *a, int n);
struct sample_aldr_opt_s read_sample_aldr_opt(char *fname);
void free_sample_aldr_opt(struct sample_aldr_opt_s x);
size_t get_size_sample_aldr_opt_s(struct sample_aldr_opt_s *s);
int sample_aldr_opt(struct sample_aldr_opt_s *x);

#endif
