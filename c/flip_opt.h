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
  Name:     flip_opt.h
  Purpose:  Optimized random-bit stream for the *_opt samplers.

  Same principle as flip.c (a background thread fills a buffer of
  cryptographic random bits with arc4random), but:
    - 64-bit words in a 256 KB ring buffer, filled by blocks of 512 words
      (one arc4random_buf call per 4 KB instead of one call per 32 bits);
    - head and tail on separate cache lines, published with release/acquire;
    - the producer sleeps when the buffer is full instead of spinning;
    - the consumer reads bits directly in the current word (inline).
  Bits are consumed from the least significant bit of each word, as in flip.c.
  NUM_RNG_CALLS is incremented by 2 per 64-bit word, so that it keeps counting
  32-bit words as for the other samplers.
*/

#ifndef FLIP_OPT_H
#define FLIP_OPT_H

#include <stdint.h>
#include <stddef.h>
#include <sched.h>
#include <gmp.h>

#define FO_WORDS (1u << 15)     // 32768 words of 64 bits = 256 KB
#define FO_MASK (FO_WORDS - 1)
#define FO_CHUNK 512            // words written per arc4random_buf call

extern unsigned long NUM_RNG_CALLS;

typedef struct {
    _Alignas(128) uint64_t tail;     // written by the producer (release)
    _Alignas(128) uint64_t head;     // published by the consumer (release)
    _Alignas(128) uint64_t c_head;   // consumer: current word
    uint64_t c_tail;                 // consumer: last tail seen
    uint32_t c_off;                  // consumer: bits already read in the current word
    uint32_t c_publish;              // consumer: words consumed since last publication
    _Alignas(128) volatile int running;
    _Alignas(128) uint64_t data[FO_WORDS];
} fo_queue_t;

extern fo_queue_t fo_q;

void fo_start(void);   // initialize the queue and start the producer thread
void fo_stop(void);    // stop and join the producer thread

void fo_wait_slow(void);

// Ensures the current word is available.
static inline void fo_ready(void) {
    if (__builtin_expect(fo_q.c_head == fo_q.c_tail, 0)) {
        fo_wait_slow();
    }
}

// Current word shifted so that bit 0 is the next unread bit; *avail = unread bits in it.
static inline uint64_t fo_peek(uint32_t *avail) {
    fo_ready();
    *avail = 64 - fo_q.c_off;
    return fo_q.data[fo_q.c_head & FO_MASK] >> fo_q.c_off;
}

// Consumes nb <= avail bits of the current word.
static inline void fo_skip(uint32_t nb) {
    fo_q.c_off += nb;
    if (fo_q.c_off == 64) {
        fo_q.c_off = 0;
        fo_q.c_head++;
        NUM_RNG_CALLS += 2;
        if (++fo_q.c_publish == 64) {
            fo_q.c_publish = 0;
            __atomic_store_n(&fo_q.head, fo_q.c_head, __ATOMIC_RELEASE);
        }
    }
}

// nb random bits, 1 <= nb <= 64.
static inline uint64_t fo_bits(uint32_t nb) {
    uint32_t avail;
    uint64_t w = fo_peek(&avail);
    if (nb <= avail) {
        fo_skip(nb);
        return (nb == 64) ? w : (w & ((1ull << nb) - 1));
    }
    uint64_t lo = w;                 // avail < 64 bits
    fo_skip(avail);
    uint32_t rest = nb - avail;
    uint64_t hi = fo_bits(rest);
    return lo | (hi << avail);
}

static inline uint32_t fo_bit(void) {
    uint32_t avail;
    uint64_t w = fo_peek(&avail);
    fo_skip(1);
    return (uint32_t)(w & 1);
}

// Reads coins until the first 1, at most `cap` coins.
// Returns the index (1..cap) of the first 1, or 0 if the `cap` coins are all 0.
static inline uint32_t fo_coins_until_one(uint32_t cap) {
    uint32_t consumed = 0;
    for (;;) {
        uint32_t avail;
        uint64_t w = fo_peek(&avail);
        uint32_t look = avail;
        if (look > cap - consumed) {
            look = cap - consumed;
        }
        if (w != 0) {
            uint32_t tz = (uint32_t)__builtin_ctzll(w);
            if (tz < look) {
                fo_skip(tz + 1);
                return consumed + tz + 1;
            }
        }
        fo_skip(look);
        consumed += look;
        if (consumed == cap) {
            return 0;
        }
    }
}

// Fast Dice Roller (Lumbroso) on [0, n), nbits = bit length of n-1 (precomputed).
static inline uint32_t fo_uniform(uint32_t n, uint32_t nbits) {
    if (n <= 1) {
        return 0;
    }
    uint64_t bound = 1ull << nbits;
    uint64_t x = fo_bits(nbits);
    while (bound < n || x >= n) {
        if (bound >= n) {
            bound -= n;
            x -= n;
        }
        bound <<= 1;
        x = (x << 1) | fo_bit();
    }
    return (uint32_t)x;
}

static inline uint32_t fo_nbits(uint32_t n) {
    return (n > 1) ? 32 - __builtin_clz(n - 1) : 0;
}


// *********************************************************************************
//   Bernoulli(t/d) with precomputed binary digits (Lumbroso's algorithm, exact)
// *********************************************************************************

// Number of precomputed digits (1..31); smaller values only for testing.
#ifndef FO_DIGITS
#define FO_DIGITS 31
#endif

// Packs digits d_1..d_FO_DIGITS of t/d (0 < t < d) at bits 31..(32-FO_DIGITS),
// bit 0 = 1 if the expansion is finite within these digits. tmp, rem: scratch.
uint32_t fo_digits(const mpz_t t, const mpz_t d, mpz_t tmp, mpz_t rem);

// Continuation of Lumbroso's Bernoulli(t/d) after FO_DIGITS steps.
uint32_t fo_bernoulli_tail(const mpz_t t, const mpz_t d);

// Result of the fast part: 0 or 1, or 2 when the slow path (probability 2^-31) is needed.
// Step k: if the expansion ends at digit k, return a coin; else if coin = 1 return d_k.
static inline uint32_t fo_bernoulli_digits(uint32_t bits) {
    uint32_t finite = bits & 1;
    uint32_t cap = finite ? 32 - __builtin_ctz(bits & ~1u) : FO_DIGITS;
    uint32_t k = fo_coins_until_one(cap);
    if (k) {
        return (bits >> (32 - k)) & 1;
    }
    return finite ? 0 : 2;
}

#endif
