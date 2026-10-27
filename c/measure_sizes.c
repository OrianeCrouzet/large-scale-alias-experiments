/*
  Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
  SPDX-License-Identifier: GPL-3.0-or-later
*/

/*
  measure_sizes -- memory footprint of the data structure of one sampler.

    measure_sizes SAMPLER FILE.dist OUT_FILE

  Builds the sampler from FILE.dist and writes, in OUT_FILE, the number of
  bytes of its data structure: the structure itself plus the useful bytes of
  its arrays and integers (see each get_size_* function).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aldr_opt.h"
#include "alias_fractions_opt.h"
#include "alias_integers_opt.h"
#include "alias_integers_v2.h"
#include "alias_rust_opt.h"
#include "fldr_opt.h"

unsigned long NUM_RNG_CALLS = 0;
unsigned long REJET = 0;

#define SIZE_OF(key, struct_name, func_read, func_size, func_free) \
    if (strcmp(sampler, key) == 0) {                               \
        struct struct_name s = func_read(path);                    \
        size = func_size(&s);                                      \
        func_free(s);                                              \
    }

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s sampler path out_file\n", argv[0]);
        return EXIT_FAILURE;
    }
    char *sampler = argv[1];
    char *path = argv[2];
    size_t size = 0;

    SIZE_OF("alias.integers_opt", sample_alias_integers_opt_s, read_sample_alias_integers_opt,
            get_size_sample_alias_integers_opt_s, free_sample_alias_integers_opt)
    else SIZE_OF("alias.integers_v2", sample_alias_integers_v2_s, read_sample_alias_integers_v2,
                 get_size_sample_alias_integers_v2_s, free_sample_alias_integers_v2)
    else SIZE_OF("rej.enc_opt", sample_fldr_opt_s, read_sample_fldr_opt,
                 get_size_sample_fldr_opt_s, free_sample_fldr_opt)
    else SIZE_OF("aldr_opt", sample_aldr_opt_s, read_sample_aldr_opt,
                 get_size_sample_aldr_opt_s, free_sample_aldr_opt)
    else SIZE_OF("alias.rust_opt", sample_alias_rust_opt_s, read_sample_alias_rust_opt,
                 get_size_sample_alias_rust_opt_s, free_sample_alias_rust_opt)
    else SIZE_OF("alias.fractions_opt", sample_alias_fractions_opt_s, read_sample_alias_fractions_opt,
                 get_size_sample_alias_fractions_opt_s, free_sample_alias_fractions_opt)
    else {
        fprintf(stderr, "Unknown sampler: %s\n", sampler);
        return EXIT_FAILURE;
    }

    FILE *out = fopen(argv[3], "w");
    if (out == NULL) {
        perror(argv[3]);
        return EXIT_FAILURE;
    }
    fprintf(out, "%zu\n", size);
    fclose(out);
    return EXIT_SUCCESS;
}
