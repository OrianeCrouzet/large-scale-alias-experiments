/*
  Derived from the fast-loaded-dice-roller-experiments repository of
  F. A. Saad et al. (https://github.com/probsys/fast-loaded-dice-roller-experiments),
  Copyright (C) 2019-2020 MIT Probabilistic Computing Project, F. A. Saad and C. E. Freer,
  released under the Apache License 2.0 (see LICENSES/Apache-2.0.txt).

  The construction is a C port, to arbitrary-precision integers, of the alias
  table of the Rust crate weighted_rand (https://github.com/ichi-h/weighted_rand),
  Copyright (c) 2020 Ippee, released under the MIT License
  (see LICENSES/MIT-weighted_rand.txt).
  Modified in 2026 by Oriane Crouzet, Antoine Genitrini and Mehdi Naima: port of weighted_rand to GMP integers;
  extracted from construct.c and readio.c into a stand-alone file.
  This file is distributed under the GNU General Public License v3.0 or later
  as part of Alias_experiments (see LICENSE).
  SPDX-License-Identifier: GPL-3.0-or-later
*/

/* alias_rust_gmp.h -- Rust-Alias table (weighted_rand) on GMP integers. */

#ifndef ALIAS_RUST_GMP_H
#define ALIAS_RUST_GMP_H

#include <stdint.h>
#include <gmp.h>

#include "vector_int.h"
#include "vector_mpz.h"

typedef enum {
    WEIGHTED_OK = 0,
    WEIGHTED_NO_ITEM,
    WEIGHTED_TOO_MANY,
    WEIGHTED_INVALID_WEIGHT,
    WEIGHTED_ALL_ZERO
} WeightedError;

typedef struct {
    VectorInt aliases;
    VectorMpz prob;
    VectorInt small;
    VectorInt large;
    unsigned int n;
    mpz_t weight_sum;
} alias_rust_gmp_s;

struct sample_alias_rust_gmp_s {
    VectorInt aliases;
    VectorMpz prob;
    VectorInt small;
    VectorInt large;
    unsigned int n;
    mpz_t weight_sum;
};

WeightedError weighted_alias_new_gmp(alias_rust_gmp_s *out, const mpz_t *weights, uint32_t n);
struct sample_alias_rust_gmp_s read_sample_alias_rust_gmp(char *fname);
void free_sample_alias_rust_gmp_s(struct sample_alias_rust_gmp_s x);

#endif
