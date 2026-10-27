/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
*/
/* counters.h -- global counters of the measurement programs. */

#ifndef COUNTERS_H
#define COUNTERS_H

/* Words of 32 random bits consumed (incremented by flip_opt.h). */
extern unsigned long NUM_RNG_CALLS;
/* Rejections of the virtual object (Integer-Alias). */
extern unsigned long REJET;

#endif
