/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
  This file is part of Alias_experiments, distributed under the GNU General
  Public License v3.0 or later (see LICENSE).
*/
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <stdio.h>
#include "vector_mpz.h"

void vector_mpz_init(VectorMpz* v) {
    v->capacity = 16;
    v->size = 0;
    v->data = malloc(v->capacity * sizeof(mpz_t));
    if (!v->data) {
        fprintf(stderr, "malloc error in vector_mpz_init\n");
        exit(EXIT_FAILURE);
    }
}

void vector_mpz_push(VectorMpz* vec, const mpz_t val) {
    if (vec->size >= vec->capacity) {
        size_t new_cap = vec->capacity ? vec->capacity * 2 : 4;
        mpz_t *new_data = realloc(vec->data, new_cap * sizeof(mpz_t));
        if (!new_data) {
            fprintf(stderr, "realloc error in vector_mpz_push\n");
            exit(EXIT_FAILURE);
        }
        vec->data = new_data;
        vec->capacity = new_cap;
    }

    // Initializes and copies only the NEW element, at index vec->size
    mpz_init(vec->data[vec->size]);
    mpz_set(vec->data[vec->size], val);
    vec->size++;
}

void vector_mpz_get(mpz_t dest, const VectorMpz* v, unsigned int index) {
    mpz_set(dest, v->data[index]);
}

void vector_mpz_set(VectorMpz* v, unsigned int index, const mpz_t val) {
    mpz_set(v->data[index], val);
}

void vector_mpz_swap(VectorMpz* v, unsigned int i, unsigned int j) {
    mpz_swap(v->data[i], v->data[j]);
}

void vector_mpz_free(VectorMpz* v) {
    if (!v->data) return;
    
    // Clears (mpz_clear) ONLY the inserted elements (0 to v->size - 1)
    for (size_t i = 0; i < v->size; i++) {
        mpz_clear(v->data[i]);
    }
    free(v->data);
    v->data = NULL;
    v->size = 0;
    v->capacity = 0;
}

void vector_mpz_reserve(VectorMpz *vec, const mpz_t min_capacity)
{
    if (!mpz_fits_ulong_p(min_capacity)) {
        return;
    }

    size_t min = (size_t) mpz_get_ui(min_capacity);
    if (vec->capacity >= min) return;

    size_t new_capacity = (vec->capacity == 0) ? 1 : vec->capacity;
    while (new_capacity < min) {
        new_capacity *= 2;
    }

    mpz_t *temp = realloc(vec->data, new_capacity * sizeof(mpz_t));
    if (temp == NULL) return;

    vec->data = temp;
    vec->capacity = new_capacity;
}