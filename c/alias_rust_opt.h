/*
  The table is the one of the Rust crate weighted_rand
  (https://github.com/ichi-h/weighted_rand), Copyright (c) 2020 Ippee,
  released under the MIT License (see LICENSES/MIT-weighted_rand.txt),
  built by alias_rust_gmp.c.

  Modified in 2026 by Oriane Crouzet, Antoine Genitrini and Mehdi Naima:
  compact storage of the table and sampling step on GMP integers.
  This file is distributed under the GNU General Public License v3.0 or later
  as part of Alias_experiments (see LICENSE).
  SPDX-License-Identifier: GPL-3.0-or-later
*/
/*
  Name:     alias_rust_opt.h
  Purpose:  Optimized Rust-Alias sampler (alias.rust_opt).
            Same table as alias.rust_gmp (weighted_alias_new_gmp: weights scaled
            by n, cells of capacity W = sum of the weights) and same sampling:
            a uniform cell i, then a uniform integer r in [0, W) drawn with the
            Fast Dice Roller in arbitrary precision, and i if r < prob[i], its
            alias otherwise. Optimizations: no allocation per sample
            (preallocated GMP scratch), random limbs written directly, the
            thresholds stored as fixed-width limbs, the uniform over the cells
            on machine integers, and the bit stream of flip_opt.h.
*/

#ifndef ALIAS_RUST_OPT_H
#define ALIAS_RUST_OPT_H

#include <stddef.h>
#include <stdint.h>
#include <gmp.h>

struct sample_alias_rust_opt_s {
    uint32_t n;
    uint32_t n_nbits;      // bit length of n - 1
    uint32_t W_nbits;      // bit length of W - 1
    uint32_t limbs;        // limbs per threshold (= limbs of W)
    uint32_t *aliases;
    mp_limb_t *prob;       // n thresholds of `limbs` limbs each
    mpz_t W;               // weight_sum
    mpz_t x, bound;        // scratch of the uniform on [0, W)
};

struct sample_alias_rust_opt_s read_sample_alias_rust_opt(char *fname);
// Construction from the weights w_1..w_n (w is not modified).
struct sample_alias_rust_opt_s preprocess_alias_rust_opt(mpz_t *w, uint32_t n);
void free_sample_alias_rust_opt(struct sample_alias_rust_opt_s x);
size_t get_size_sample_alias_rust_opt_s(struct sample_alias_rust_opt_s *s);
uint32_t sample_alias_rust_opt(struct sample_alias_rust_opt_s *x);

#endif
