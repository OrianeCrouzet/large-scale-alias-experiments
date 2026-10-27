/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
  This file is part of Alias_experiments, distributed under the GNU General
  Public License v3.0 or later (see LICENSE).
*/
/*
  Name:     vector_int.h
  Purpose:  Vector implementation for integers.
  Author:   Oriane Crouzet
  Copyright (C) 2026 Oriane Crouzet, All Rights Reserved.
*/

#ifndef VECTOR_INT_H
#define VECTOR_INT_H

typedef struct {
    int* data;
    unsigned int size;
    unsigned int capacity;
} VectorInt;

void vector_init(VectorInt* v);
void vector_push(VectorInt* v, int value);
int vector_get(VectorInt* v, unsigned int index);
void vector_free(VectorInt* v);
void vector_resize_zero(VectorInt* vec, unsigned int new_size);
void vector_swap(VectorInt* v, unsigned int i, unsigned int j);
void vector_reserve(VectorInt* vec, unsigned int min_capacity);
void vector_set(VectorInt* v, unsigned int index, unsigned int val);

#endif
