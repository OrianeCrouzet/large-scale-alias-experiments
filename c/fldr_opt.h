/*
  Derived from the fast-loaded-dice-roller-experiments repository of
  F. A. Saad et al. (https://github.com/probsys/fast-loaded-dice-roller-experiments),
  Copyright (C) 2019-2020 MIT Probabilistic Computing Project, F. A. Saad and C. E. Freer,
  released under the Apache License 2.0 (see LICENSES/Apache-2.0.txt).

  Modified in 2026 by Oriane Crouzet, Antoine Genitrini and Mehdi Naima:
  FLDR on arbitrary-precision weights; C port of the Python construction of
  the DDG tree (rejection.py, matrix.py, tree.py, packing.py); walk of the
  encoding inside the current word of the random-bit stream.
  This file is distributed under the GNU General Public License v3.0 or later
  as part of Alias_experiments (see LICENSE).
  SPDX-License-Identifier: GPL-3.0-or-later
*/
/*
  Name:     fldr_opt.h
  Purpose:  Optimized FLDR sampler (rej.enc_opt).
            DDG-tree encoding of Saad et al. (pack_tree), built in C from the
            arbitrary-precision weights by preprocess_fldr_opt, and walked as in
            Saad's sample_rejection_encoding.
            The encoding holds memory locations and labels, stored as 32-bit
            integers; the tree is walked inside the current word of the
            random-bit stream (flip_opt.h).
*/

#ifndef FLDR_OPT_H
#define FLDR_OPT_H

#include <stddef.h>
#include <stdint.h>
#include <gmp.h>

struct sample_fldr_opt_s {
    uint32_t n;          // labels >= n are rejected (same header as rej.enc_gmp)
    uint32_t length;     // length of the encoding
    int32_t *encoding;   // pack_tree: leaf = [-label], internal node = [loc(left), loc(right)]
};

// Reads a .dist file and builds the encoding (preprocess_fldr_opt).
struct sample_fldr_opt_s read_sample_fldr_opt(char *fname);
// Construction from the probabilities w_i / Z (w is not modified): same encoding as
// construct_sample_rejection_encoding_gmp + pack_tree of the original Python code,
// written in C. Z is the first line of the .dist file; it equals sum(w) except for
// the uniform distribution of the experiments, where sum(w) exceeds it slightly.
struct sample_fldr_opt_s preprocess_fldr_opt(mpz_t *w, uint32_t n, const mpz_t Z);
void free_sample_fldr_opt(struct sample_fldr_opt_s x);
size_t get_size_sample_fldr_opt_s(struct sample_fldr_opt_s *s);
int sample_fldr_opt(struct sample_fldr_opt_s *x);

#endif
