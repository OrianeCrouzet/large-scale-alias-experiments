# Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Measurement programs of the experiments.
#
#   make          builds measure_runtime, measure_sizes, measure_preprocess, check_samplers
#   make clean

CC     ?= cc
CFLAGS ?= -O3
CFLAGS += -std=c11 -Wall -Wno-unused-result
GMP_CFLAGS := $(shell pkg-config --cflags gmp 2>/dev/null)
GMP_LIBS   := $(shell pkg-config --libs gmp 2>/dev/null || echo -lgmp)
CFLAGS += $(GMP_CFLAGS)
LDLIBS += $(GMP_LIBS) -lpthread -lm

SAMPLERS := flip_opt.c vector_int.c vector_mpz.c \
            fldr_opt.c aldr_opt.c \
            alias_rust_gmp.c alias_rust_opt.c \
            alias_fractions_gmp.c alias_fractions_opt.c \
            alias_integers_opt.c alias_integers_opt_sample.c \
            alias_integers_v2.c alias_integers_v2_sample.c
HEADERS := $(wildcard *.h)
PROGRAMS := measure_runtime measure_sizes measure_preprocess check_samplers

all: $(PROGRAMS)

$(PROGRAMS): %: %.c $(SAMPLERS) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $< $(SAMPLERS) $(LDLIBS)

clean:
	rm -f $(PROGRAMS)

.PHONY: all clean
