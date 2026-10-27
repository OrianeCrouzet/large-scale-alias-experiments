/*
  Derived from the fast-loaded-dice-roller-experiments repository of
  F. A. Saad et al. (https://github.com/probsys/fast-loaded-dice-roller-experiments),
  Copyright (C) 2019-2020 MIT Probabilistic Computing Project, F. A. Saad and C. E. Freer,
  released under the Apache License 2.0 (see LICENSES/Apache-2.0.txt).

  Modified in 2026 by Oriane Crouzet, Antoine Genitrini and Mehdi Naima:
  new producer (ring buffer of 64-bit words, filled by blocks; the producer
  sleeps when the buffer is full), following the principle of flip.c.
  This file is distributed under the GNU General Public License v3.0 or later
  as part of Alias_experiments (see LICENSE).
  SPDX-License-Identifier: GPL-3.0-or-later
*/
/*
  Name:     flip_opt.c
  Purpose:  Producer thread of the optimized random-bit stream (see flip_opt.h).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

#include "flip_opt.h"

fo_queue_t fo_q;
static pthread_t fo_thread;

static void *fo_producer(void *arg) {
    (void)arg;
    while (fo_q.running) {
        uint64_t t = fo_q.tail;
        uint64_t h = __atomic_load_n(&fo_q.head, __ATOMIC_ACQUIRE);
        if (FO_WORDS - (t - h) >= FO_CHUNK) {
            // FO_WORDS is a multiple of FO_CHUNK: a chunk never wraps around.
            arc4random_buf(&fo_q.data[t & FO_MASK], FO_CHUNK * sizeof(uint64_t));
            __atomic_store_n(&fo_q.tail, t + FO_CHUNK, __ATOMIC_RELEASE);
        } else {
            // Full buffer (~2 million bits): sleep instead of spinning.
            usleep(50);
        }
    }
    return NULL;
}

void fo_start(void) {
    memset(&fo_q, 0, sizeof(fo_q));
    fo_q.running = 1;
    NUM_RNG_CALLS = 0;
    if (pthread_create(&fo_thread, NULL, fo_producer, NULL) != 0) {
        perror("fo_start: pthread_create");
        exit(EXIT_FAILURE);
    }
}

void fo_stop(void) {
    fo_q.running = 0;
    pthread_join(fo_thread, NULL);
}

// Called when every known word has been consumed.
void fo_wait_slow(void) {
    __atomic_store_n(&fo_q.head, fo_q.c_head, __ATOMIC_RELEASE);
    fo_q.c_publish = 0;
    for (;;) {
        fo_q.c_tail = __atomic_load_n(&fo_q.tail, __ATOMIC_ACQUIRE);
        if (fo_q.c_head != fo_q.c_tail) {
            return;
        }
        sched_yield();
    }
}


// *********************************************************************************
//   Bernoulli with precomputed digits
// *********************************************************************************

uint32_t fo_digits(const mpz_t t, const mpz_t d, mpz_t tmp, mpz_t rem) {
    // digits = floor(t * 2^FO_DIGITS / d) < 2^FO_DIGITS, rem = t * 2^FO_DIGITS mod d
    mpz_mul_2exp(tmp, t, FO_DIGITS);
    mpz_fdiv_qr(tmp, rem, tmp, d);
    return ((uint32_t)mpz_get_ui(tmp) << (32 - FO_DIGITS)) | (mpz_sgn(rem) == 0);
}

uint32_t fo_bernoulli_tail(const mpz_t t, const mpz_t d) {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mul_2exp(tmp, t, FO_DIGITS);
    mpz_fdiv_r(tmp, tmp, d);
    uint32_t y;
    for (;;) {
        mpz_mul_2exp(tmp, tmp, 1);
        if (mpz_cmp(tmp, d) == 0) {
            y = fo_bit();
            break;
        }
        y = (mpz_cmp(tmp, d) > 0);
        if (y) {
            mpz_sub(tmp, tmp, d);
        }
        if (fo_bit()) {
            break;
        }
    }
    mpz_clear(tmp);
    return y;
}
