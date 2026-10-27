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

/* alias_fractions_gmp.c -- construction of the Fraction-Alias table. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alias_fractions_gmp.h"

// *********************************************************************************
//              ALIAS FRACTIONS - GMP (arbitrary-size integers)
// *********************************************************************************

PileResult piles(mpq_t* distrib, int N, int affiche) {
    PileResult res;

    // pdsCase = 1 / N
    mpq_init(res.pdsCase);
    mpq_set_ui(res.pdsCase, 1, N);

    res.S0 = malloc(N * sizeof(int));
    res.S1 = malloc(N * sizeof(int));
    for (int i = 0; i < N; i++) {
        res.S0[i] = -1;
        res.S1[i] = -1;
    }

    res.lenS0 = 0;
    res.lenS1 = 0;

    for (int n = 0; n < N; n++) {
        if (mpq_cmp(distrib[n], res.pdsCase) > 0) {
            res.S0[res.lenS0++] = n;
        } else {
            res.S1[res.lenS1++] = n;
        }
    }

    if (affiche) {
        printf("S0 (heavy):\n");
        for (int i = 0; i < res.lenS0; i++) {
            int idx = res.S0[i];
            gmp_printf("(%d, %Qd)\n", idx, distrib[idx]);
        }
    }

    return res;
}

struct AliasEntryGMP* algo_alias_fractions_gmp(mpq_t* distrib, int N, int affiche) {
    // piles(distrib) -> pdsCase, S0, S1
    PileResult res = piles(distrib, N, 0);
    mpq_t pdsCase;
    mpq_init(pdsCase);
    mpq_set(pdsCase, res.pdsCase);  // pdsCase = 1/N

    int* S0 = res.S0;
    int lenS0 = res.lenS0;
    int* S1 = res.S1;
    int lenS1 = res.lenS1;

    if (affiche) {
        gmp_printf("pdsCase = %Qd\n", pdsCase);
        printf("heavy (S0): ");
        for (int i = 0; i < lenS0; i++) printf("%d ", S0[i]);
        printf("\nlight (S1): ");
        for (int i = 0; i < lenS1; i++) printf("%d ", S1[i]);
        printf("\n");
    }

    struct AliasEntryGMP* T = malloc(N * sizeof(struct AliasEntryGMP));

    int idx0 = 0;  // for S0
    int idx1 = 0;  // for S1

    for (int t = 0; t < N; t++) {
         if ((lenS0 - idx0 == 0) && (lenS1 - idx1 == 0)) {
            fprintf(stderr, "ERROR2: no heavy and no light objects left at step %d!\n", t);
            exit(EXIT_FAILURE);
        }

        if (lenS1 - idx1 > 0) {
            int moins = S1[idx1];

            if (mpq_cmp(distrib[moins], pdsCase) == 0) {
                // T[t] = (light, None, pdsCase)
                T[t].i = moins;
                T[t].j = -1;
                mpq_init(T[t].prob);
                mpq_set(T[t].prob, pdsCase);
                idx1++;
            } else {
                if (idx0 >= lenS0) {
                    fprintf(stderr, "[ERROR] idx0=%d exceeds lenS0=%d at step %d!\n", idx0, lenS0, t);
                    exit(EXIT_FAILURE);
                }
                int plus = S0[idx0];

                T[t].i = moins;
                T[t].j = plus;
                mpq_init(T[t].prob);
                mpq_set(T[t].prob, distrib[moins]);

                idx1++;  // S1 = S1[1:]

                // distrib[heavy] -= (pdsCase - distrib[light])
                mpq_t diff;
                mpq_init(diff);
                mpq_sub(diff, pdsCase, distrib[moins]);
                mpq_sub(distrib[plus], distrib[plus], diff);
                mpq_clear(diff);

                // distrib[light] = 0
                mpq_set_ui(distrib[moins], 0, 1);

                if (mpq_cmp(distrib[plus], pdsCase) <= 0) {
                    if (idx1 <= 0) {
                        fprintf(stderr, "[ERROR] idx1=%d cannot be decremented at step %d!\n", idx1, t);
                        exit(EXIT_FAILURE);
                    }
                    S1[--idx1] = plus;  // S1 = [plus] + S1
                    idx0++; // S0 = S0[1:]
                }
            }
        }else {
            // Case : no more light ones, but still heavy ones
            if (idx0 >= lenS0) {
                fprintf(stderr, "[ERROR] idx0=%d exceeds lenS0=%d: no heavy object left!\n", idx0, lenS0);
                exit(1);
            }

            int only_heavy = S0[idx0];
            T[t].i = only_heavy;
            T[t].j = -1;
            mpq_init(T[t].prob);
            mpq_set(T[t].prob, pdsCase);

            mpq_sub(distrib[only_heavy], distrib[only_heavy], pdsCase);

            if (mpq_cmp(distrib[only_heavy], pdsCase) <= 0) {
                if (idx1 <= 0) {
                    fprintf(stderr, "[ERROR] idx1=%d cannot be decremented in the only_heavy case at step %d!\n", idx1, t);
                    exit(EXIT_FAILURE);
                }
                S1[--idx1] = only_heavy;  // S1 = [only_heavy] + S1
                idx0++;
            }

        }
    }

    // prob = k / pdsCase
    for (int t = 0; t < N; t++) {
        mpq_div(T[t].prob, T[t].prob, pdsCase);
    }

    // Cleaning
    mpq_clear(pdsCase);
    mpq_clear(res.pdsCase);
    free(res.S0);
    free(res.S1);

    if (affiche) {
        printf("T:\n");
        for (int t = 0; t < N; t++) {
            gmp_printf("(%d, %d, %Qd)\n", T[t].i, T[t].j, T[t].prob);
        }
    }

    return T;
}

struct sample_alias_fractions_gmp_s preprocess_alias_fractions_gmp(mpz_t* a, int n) {
    // Rational distribution a[i] / sum(a)
    mpq_t* distrib = malloc(n * sizeof(mpq_t));
    mpq_t total;
    mpq_init(total);
    mpq_set_ui(total, 0, 1);

    for (int i = 0; i < n; i++) {
        mpq_init(distrib[i]);
        mpq_set_z(distrib[i], a[i]);  // distrib[i] = a[i] / 1
        mpq_add(total, total, distrib[i]);
    }

    // Normaliser : distrib[i] = distrib[i] / total
    for (int i = 0; i < n; i++) {
        mpq_div(distrib[i], distrib[i], total);
    }

    mpq_clear(total);

    // Appel de l'algo principal
    struct AliasEntryGMP* T = algo_alias_fractions_gmp(distrib, n, 0);    // affichage = False

    // Free the distribution
    for (int i = 0; i < n; i++) {
        mpq_clear(distrib[i]);
    }
    free(distrib);

    // Build and return the structure
    struct sample_alias_fractions_gmp_s result;
    result.taille = n;
    result.table = T;

    return result;
}

struct sample_alias_fractions_gmp_s read_sample_alias_fractions_gmp(char *fname){
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
        perror(fname);
        exit(EXIT_FAILURE);
    }
    mpz_t Z;
    mpz_init(Z);

    char *line = NULL;
    size_t len = 0;

    if (getline(&line, &len, fp) == -1) {
        fprintf(stderr, "[Error] empty file or cannot read Z in %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Converts the first line (base 10) into Z
    if (mpz_set_str(Z, line, 10) != 0) {
        fprintf(stderr, "[Error] cannot convert Z in %s\n", fname);
        exit(EXIT_FAILURE);
    }

    int n = 0;
    if (fscanf(fp, "%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "[Error] invalid n in %s\n", fname);
        exit(EXIT_FAILURE);
    }

    mpz_t* array = malloc(n * sizeof(mpz_t));
    char str_buf[1024]; // buffer for each big integer

    for (int i = 0; i < n; ++i) {
        mpz_init(array[i]);
        // %1023s reads the next word, skipping spaces and newlines
        if (fscanf(fp, "%1023s", str_buf) != 1) {
            fprintf(stderr, "[Error] cannot read weight %d\n", i);
            exit(EXIT_FAILURE);
        }
        if (mpz_set_str(array[i], str_buf, 10) != 0) {
            fprintf(stderr, "[Error] cannot convert weight %d to mpz_t\n", i);
            exit(EXIT_FAILURE);
        }
    }

    free(line);
    fclose(fp);

    struct sample_alias_fractions_gmp_s sampler = preprocess_alias_fractions_gmp(array, n);

    mpz_clear(Z);
    for (int i = 0; i < n; ++i) {
        mpz_clear(array[i]);
    }
    free(array);

    return sampler;
}

void free_sample_alias_fractions_gmp_s(struct sample_alias_fractions_gmp_s x) {
   for (int i = 0; i < x.taille; i++) {
        mpq_clear(x.table[i].prob);  // frees the mpq_t
    }
    free(x.table);  // frees the table
}
