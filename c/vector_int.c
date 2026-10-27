/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
  This file is part of Alias_experiments, distributed under the GNU General
  Public License v3.0 or later (see LICENSE).
*/
/*
  Name:     vector_int.c
  Purpose:  Vector implementation for integers.
  Author:   Oriane Crouzet
  Copyright (C) 2026 Oriane Crouzet, All Rights Reserved.
*/

#include <stdlib.h>
#include <stdio.h>

#include "vector_int.h"

void vector_init(VectorInt* v) {
    v->capacity = 16;
    v->size = 0;
    v->data = malloc(v->capacity * sizeof(int));
}

void vector_push(VectorInt* v, int value) {
    if (v->size >= v->capacity) {
        v->capacity *= 2;
        v->data = realloc(v->data, v->capacity * sizeof(int));
        
    }
    v->data[v->size++] = value;
}

int vector_get(VectorInt* v, unsigned int index) {
    return v->data[index];
}

void vector_set(VectorInt* v, unsigned int index, unsigned int val) {
    v->data[index] = val;
}

void vector_free(VectorInt* v) {
    free(v->data);
    v->data = NULL;
    v->size = 0;
    v->capacity = 0;
}

void vector_resize_zero(VectorInt* vec, unsigned int new_size) {
    if (new_size > vec->capacity) {
        vec->capacity = new_size * 2; // some room for the next pushes
        vec->data = realloc(vec->data, vec->capacity * sizeof(int));
    }

    for (unsigned int i = vec->size; i < new_size; i++) {
        vec->data[i] = 0;
    }

    vec->size = new_size;
}

void vector_swap(VectorInt* v, unsigned int i, unsigned int j) {
    int temp = v->data[i];
    v->data[i] = v->data[j];
    v->data[j] = temp;
}

// Preallocates a minimum capacity
void vector_reserve(VectorInt* vec, unsigned int min_capacity) {
    if (vec->capacity < min_capacity) {
        unsigned int new_capacity = vec->capacity;
        while (new_capacity < min_capacity) new_capacity *= 2;
        
        int* temp = realloc(vec->data, new_capacity * sizeof(int));
        
        vec->data = temp;
        vec->capacity = new_capacity;
    }
}

