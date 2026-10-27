# Notices

Alias_experiments is distributed under the GNU General Public License v3.0 or
later (see `LICENSE`). It contains code derived from the following projects,
whose licenses are compatible with the GPL v3 and whose notices are retained
in the derived files and in `LICENSES/`.

## fast-loaded-dice-roller-experiments

* Source: <https://github.com/probsys/fast-loaded-dice-roller-experiments>
  (experiments of F. A. Saad, C. E. Freer, M. C. Rinard, V. K. Mansinghka,
  *The Fast Loaded Dice Roller: A Near-Optimal Exact Sampler for Discrete
  Probability Distributions*, AISTATS 2020).
* Copyright (C) 2019–2020 MIT Probabilistic Computing Project, F. A. Saad and
  C. E. Freer.
* License: Apache License 2.0 (`LICENSES/Apache-2.0.txt`).
* Derived files and changes:
  * `c/measure_runtime.c`: from `c/main.c` and `c/macros.c` — restricted to
    the six samplers of the paper, reading `.dist` files.
  * `c/fldr_opt.{h,c}`: FLDR — C port, to arbitrary-precision weights, of the
    Python construction of the DDG tree (`src/rejection.py`, `matrix.py`,
    `tree.py`, `packing.py`), and a new walk of the encoding.
  * `c/flip_opt.{h,c}`: random-bit stream following the principle of
    `c/flip.c` (producer thread), rewritten.
  * `c/alias_rust_gmp.{h,c}`, `c/alias_fractions_gmp.{h,c}`: extracted from
    `c/construct.c` and `c/readio.c` (code added to these files by O. Crouzet).
  * `scripts/generate_distributions.py`: from `experiments/dists.py`,
    `experiments/create_distributions.py`, `src/entropy.py` and
    `src/utils.py` — arbitrary-precision total weight, uniform distribution
    added, only the `.dist` files written.

## amplified-loaded-dice-roller

* Source: <https://github.com/probsys/amplified-loaded-dice-roller>
  (T. L. Draper and F. A. Saad, *Efficient Rejection Sampling in the
  Entropy-Optimal Range*).
* Copyright (C) 2025 CMU Probabilistic Computing Systems Lab.
* License: Apache License 2.0 (`LICENSES/Apache-2.0.txt`).
* Derived files and changes: `c/aldr_opt.{h,c}` — ALDR on arbitrary-precision
  weights (GMP), 32-bit storage of the tree, new walk.

## weighted_rand

* Source: <https://github.com/ichi-h/weighted_rand> (Rust crate).
* Copyright (c) 2020 Ippee.
* License: MIT License or Apache License 2.0, at the user's option; the MIT
  License is used here (`LICENSES/MIT-weighted_rand.txt`).
* Derived files and changes: `c/alias_rust_gmp.{h,c}` (C port of the alias
  table construction, on GMP integers) and `c/alias_rust_opt.{h,c}`.

## Other files

All other files (`c/alias_integers_*`, `c/alias_fractions_opt.*`,
`c/vector_*`, `c/measure_sizes.c`, `c/measure_preprocess.c`,
`c/check_samplers.c`, `c/counters.h`, `scripts/pipeline.sh`,
`scripts/make_figures.py`) are Copyright (C) 2026 Oriane Crouzet, Antoine
Genitrini and Mehdi Naima, under the GPL v3 or later.

This notice describes the provenance of the code; it is not legal advice.
