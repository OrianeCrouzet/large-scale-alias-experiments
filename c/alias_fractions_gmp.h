/*
  Derived from the fast-loaded-dice-roller-experiments repository of
  F. A. Saad et al. (https://github.com/probsys/fast-loaded-dice-roller-experiments),
  Copyright (C) 2019-2020 MIT Probabilistic Computing Project, F. A. Saad and C. E. Freer,
  released under the Apache License 2.0 (see LICENSES/Apache-2.0.txt).

  Modified in 2026 by Oriane Crouzet, Antoine Genitrini and Mehdi Naima: Fraction-Alias sampler (Walker-Vose on
  exact GMP rationals) by Oriane Crouzet; extracted from construct.c and readio.c.
  This file is distributed under the GNU General Public License v3.0 or later
  as part of Alias_experiments (see LICENSE).
  SPDX-License-Identifier: GPL-3.0-or-later
*/

/* alias_fractions_gmp.h -- Fraction-Alias table: Walker-Vose on exact rationals. */

#ifndef ALIAS_FRACTIONS_GMP_H
#define ALIAS_FRACTIONS_GMP_H

#include <gmp.h>

typedef struct {
    mpq_t pdsCase;
    int* S0;        // index where distrib[n] > pdsCase  (heavy)
    int* S1;        // index where distrib[n] <= pdsCase (light)
    int lenS0;
    int lenS1;
} PileResult;

struct AliasEntryGMP {
    int i;
    int j;          // -1 if NULL
    mpq_t prob;     // = k / pdsCase
};

struct sample_alias_fractions_gmp_s {
    int taille;                      // number of entries
    struct AliasEntryGMP* table;     // the alias table
};

PileResult piles(mpq_t* distrib, int N, int affiche);
struct AliasEntryGMP* algo_alias_fractions_gmp(mpq_t* distrib, int N, int affiche);
struct sample_alias_fractions_gmp_s preprocess_alias_fractions_gmp(mpz_t* a, int n);
struct sample_alias_fractions_gmp_s read_sample_alias_fractions_gmp(char *fname);
void free_sample_alias_fractions_gmp_s(struct sample_alias_fractions_gmp_s x);

#endif
