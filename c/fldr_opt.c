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
  Name:     fldr_opt.c
  Purpose:  Optimized FLDR sampler (see fldr_opt.h).
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "fldr_opt.h"
#include "flip_opt.h"

// .dist file: Z, then n w_1 ... w_n, then the entropy. The DDG tree is built in C
// by preprocess_fldr_opt from the probabilities w_i / Z (in the original experiments,
// it was built in Python and read from a separate file).
struct sample_fldr_opt_s read_sample_fldr_opt(char *fname) {
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
        perror(fname);
        exit(EXIT_FAILURE);
    }
    mpz_t Z;
    mpz_init(Z);
    int n = 0;
    if (mpz_inp_str(Z, fp, 10) == 0 || fscanf(fp, "%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "[Error] invalid Z or n in %s\n", fname);
        exit(EXIT_FAILURE);
    }
    mpz_t *w = malloc((size_t)n * sizeof(mpz_t));
    for (int i = 0; i < n; i++) {
        mpz_init(w[i]);
        if (mpz_inp_str(w[i], fp, 10) == 0) {
            fprintf(stderr, "[Error] cannot read weight %d in %s\n", i, fname);
            exit(EXIT_FAILURE);
        }
    }
    fclose(fp);
    struct sample_fldr_opt_s x = preprocess_fldr_opt(w, (uint32_t)n, Z);
    for (int i = 0; i < n; i++) {
        mpz_clear(w[i]);
    }
    free(w);
    mpz_clear(Z);
    return x;
}

// *********************************************************************************
//   Construction (port of src/rejection.py, matrix.py, tree.py, packing.py)
// *********************************************************************************
//
// The probabilities are p_i = w_i / Z. With g = gcd(Z, w_i), their common denominator
// is Zc = Z / g, the rows are M_i = w_i / g (i < n) and the rejection row is
// M_n = 2^k - sum(M_i), with k = ceil(log2 Zc); sum(M) = 2^k. reduce_fractions
// then halves all the M while they are even. The DDG tree of a dyadic matrix has no
// back edge: at depth d (1..k), the children of the internal nodes of depth d-1 are,
// from left to right, the internal nodes of depth d, then the leaves of the rows whose
// bit k-d is 1, by decreasing row (make_leaf_table fills the level from the right).
// pack_tree writes it in preorder, left child first: a leaf is [-(row+1)], an internal
// node is [loc(left), loc(right)].

struct fldr_build {
    uint32_t k;
    uint32_t *n_internal;     // internal nodes per depth (depth 0 = root)
    uint32_t **rows;          // rows[d]: rows with a leaf at depth d, increasing
    uint32_t *n_rows;
    int32_t *enc;
    uint32_t pos;
};

// Writes the node at position j of depth d (preorder) and returns nothing; b->pos advances.
static void pack_node(struct fldr_build *b, uint32_t d, uint32_t j) {
    if (j >= b->n_internal[d]) {
        // Leaf: positions n_internal[d] .. width-1 hold the rows in decreasing order.
        uint32_t width = b->n_internal[d] + b->n_rows[d];
        uint32_t row = b->rows[d][width - 1 - j];
        b->enc[b->pos++] = -(int32_t)(row + 1);
        return;
    }
    uint32_t here = b->pos;
    b->pos += 2;
    b->enc[here] = (int32_t)b->pos;
    pack_node(b, d + 1, 2 * j);
    b->enc[here + 1] = (int32_t)b->pos;
    pack_node(b, d + 1, 2 * j + 1);
}

struct sample_fldr_opt_s preprocess_fldr_opt(mpz_t *w, uint32_t n, const mpz_t Zfile) {
    uint32_t N = n + 1;                       // rows, rejection included
    mpz_t *M = malloc(N * sizeof(mpz_t));
    mpz_t g, Z, t;
    mpz_inits(g, Z, t, NULL);
    mpz_set(g, Zfile);
    for (uint32_t i = 0; i < n; i++) {
        mpz_gcd(g, g, w[i]);
    }
    mpz_divexact(Z, Zfile, g);                // common denominator of the w_i / Z
    // get_rejection_precision: log2(Z) if Z is a power of 2, else bit_length(Z).
    uint32_t k = (uint32_t)mpz_sizeinbase(Z, 2);
    if (mpz_popcount(Z) == 1) {
        k--;
    }
    mpz_set_ui(t, 0);
    for (uint32_t i = 0; i < n; i++) {
        mpz_init(M[i]);
        mpz_divexact(M[i], w[i], g);
        mpz_add(t, t, M[i]);
    }
    mpz_init(M[n]);
    mpz_setbit(M[n], k);
    mpz_sub(M[n], M[n], t);                   // rejection: 2^k - sum(M_i)
    if (mpz_sgn(M[n]) < 0) {
        fprintf(stderr, "[Fatal error] FLDR: sum of the probabilities > 2^k.\n");
        exit(EXIT_FAILURE);
    }

    struct sample_fldr_opt_s x;
    x.n = N;
    // reduce_fractions, case l = k (dyadic).
    for (;;) {
        int one = 0, all_even = 1, all_equal = 1;
        for (uint32_t i = 0; i < N; i++) {
            if (mpz_sizeinbase(M[i], 2) == k + 1 && mpz_popcount(M[i]) == 1) {
                one = 1;                      // M_i = 2^k
            }
            if (mpz_odd_p(M[i])) {
                all_even = 0;
            }
            if (mpz_cmp(M[i], M[0]) != 0) {
                all_equal = 0;
            }
        }
        if (one) {
            // (k', l') = (1, 0) and P = [[1]]: a single leaf of label 1.
            x.n = 1;
            x.length = 1;
            x.encoding = malloc(sizeof(int32_t));
            x.encoding[0] = -1;
            for (uint32_t i = 0; i < N; i++) {
                mpz_clear(M[i]);
            }
            free(M);
            mpz_clears(g, Z, t, NULL);
            return x;
        }
        if (k == 0) {
            break;
        }
        if (all_even) {
            for (uint32_t i = 0; i < N; i++) {
                mpz_fdiv_q_2exp(M[i], M[i], 1);
            }
            k--;
            continue;
        }
        if (all_equal && mpz_popcount(M[0]) == 1) {
            // 2^k / M_0 = 2^k' : every row becomes 1 over 2^k'.
            uint32_t kp = k - (uint32_t)(mpz_sizeinbase(M[0], 2) - 1);
            for (uint32_t i = 0; i < N; i++) {
                mpz_set_ui(M[i], 1);
            }
            k = kp;
        }
        break;
    }

    struct fldr_build b;
    b.k = k;
    b.n_internal = calloc(k + 1, sizeof(uint32_t));
    b.n_rows = calloc(k + 1, sizeof(uint32_t));
    b.rows = calloc(k + 1, sizeof(uint32_t *));
    size_t total = 0;
    for (uint32_t d = 1; d <= k; d++) {
        uint32_t cnt = 0;
        for (uint32_t i = 0; i < N; i++) {
            cnt += mpz_tstbit(M[i], k - d);
        }
        b.rows[d] = malloc((cnt + 1) * sizeof(uint32_t));
        b.n_rows[d] = 0;
        for (uint32_t i = 0; i < N; i++) {
            if (mpz_tstbit(M[i], k - d)) {
                b.rows[d][b.n_rows[d]++] = i;
            }
        }
    }
    b.n_internal[0] = 1;
    for (uint32_t d = 1; d <= k; d++) {
        uint64_t width = 2 * (uint64_t)b.n_internal[d - 1];
        if (b.n_rows[d] > width) {
            fprintf(stderr, "[Fatal error] inconsistent FLDR tree (depth %u).\n", d);
            exit(EXIT_FAILURE);
        }
        b.n_internal[d] = (uint32_t)(width - b.n_rows[d]);
        total += 2 * (size_t)b.n_internal[d - 1];
    }
    if (b.n_internal[k] != 0 || total + 2 >= INT32_MAX) {
        fprintf(stderr, "[Fatal error] inconsistent or too large FLDR tree.\n");
        exit(EXIT_FAILURE);
    }
    // Encoding: 2 cells per internal node, 1 per leaf.
    size_t length = 0;
    for (uint32_t d = 0; d < k; d++) {
        length += 2 * (size_t)b.n_internal[d];
    }
    for (uint32_t d = 1; d <= k; d++) {
        length += b.n_rows[d];
    }
    b.enc = malloc(length * sizeof(int32_t));
    b.pos = 0;
    pack_node(&b, 0, 0);

    x.length = b.pos;
    x.encoding = b.enc;
    for (uint32_t d = 1; d <= k; d++) {
        free(b.rows[d]);
    }
    free(b.rows);
    free(b.n_rows);
    free(b.n_internal);
    for (uint32_t i = 0; i < N; i++) {
        mpz_clear(M[i]);
    }
    free(M);
    mpz_clears(g, Z, t, NULL);
    return x;
}

void free_sample_fldr_opt(struct sample_fldr_opt_s x) {
    free(x.encoding);
}

size_t get_size_sample_fldr_opt_s(struct sample_fldr_opt_s *s) {
    return sizeof(*s) + (size_t)s->length * sizeof(int32_t);
}

// Walk of the encoding of pack_tree (Saad et al.), as in sample_rejection_encoding:
// c = 0; repeat { c = enc[c + bit]; if enc[c] < 0 the label is -enc[c] },
// restarting from the root when the label is >= n.
// Note: sample_ky_encoding_gmp tests enc[c + bit] < 0 instead, which reads one extra
// bit at each leaf and then the cell that follows the leaf: it is biased.
int sample_fldr_opt(struct sample_fldr_opt_s *x) {
    if (x->length == 1) {
        return 1;                     // as sample_ky_encoding
    }
    const int32_t *enc = x->encoding;
    int32_t c = 0;
    for (;;) {
        uint32_t avail;
        uint64_t w = fo_peek(&avail);
        uint32_t used = 0;
        while (used < avail) {
            c = enc[c + (int32_t)(w & 1)];
            w >>= 1;
            used++;
            int32_t e = enc[c];
            if (e < 0) {
                if ((uint32_t)(-e) < x->n) {
                    fo_skip(used);
                    return -e;
                }
                c = 0;
            }
        }
        fo_skip(used);
    }
}
