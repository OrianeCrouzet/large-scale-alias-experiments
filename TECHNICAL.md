# Alias_experiments — technical description

This document describes how the experiments of the paper *Exact Integer Alias
Method: Robust and Error-Free Large-Scale Sampling* are produced: the
distributions, the six samplers, the measurement programs, the pipeline and
the figures. [README.md](README.md) explains how to run them.

## Contents

1. [Layout](#1-layout)
2. [Distributions](#2-distributions)
3. [Samplers](#3-samplers)
4. [Random-bit stream](#4-random-bit-stream)
5. [Measurement programs](#5-measurement-programs)
6. [Pipeline](#6-pipeline)
7. [Figures](#7-figures)
8. [Validation](#8-validation)
9. [Differences with the original experiments](#9-differences-with-the-original-experiments)

---

## 1. Layout

```
c/
  flip_opt.{h,c}                 random-bit stream (§4)
  fldr_opt.{h,c}                 FLDR
  aldr_opt.{h,c}                 ALDR
  alias_rust_gmp.{h,c}           Rust-Alias table construction (port of weighted_rand)
  alias_rust_opt.{h,c}           Rust-Alias storage and sampling
  alias_fractions_gmp.{h,c}      Fraction-Alias table construction (exact rationals)
  alias_fractions_opt.{h,c}      Fraction-Alias storage and sampling
  alias_integers_opt.{h,c}       Integer-Alias table construction (Algorithm 1)
  alias_integers_opt_sample.c    Integer-Alias sampling, alea mode (v1)
  alias_integers_v2.{h,c}        Integer-Alias, speed mode (v2): digits of the thresholds
  alias_integers_v2_sample.c     Integer-Alias sampling, speed mode
  vector_int.{h,c}, vector_mpz.{h,c}   dynamic arrays used by the constructions
  counters.h                     global counters (random words, rejections)
  measure_runtime.c              sampling time (§5.1)
  measure_sizes.c                memory of the data structure (§5.2)
  measure_preprocess.c           preprocessing time (§5.3)
  check_samplers.c               chi-square test of the samplers (§5.4)
scripts/
  generate_distributions.py      distributions (§2)
  pipeline.sh                    pipeline (§6)
  make_figures.py                figures (§7)
```

Every program reads the same input, a `.dist` file, and builds the data
structure of the sampler from the weights, in C.

---

## 2. Distributions

### 2.1 Format

```
Z                   common denominator of the probabilities
n M_1 ... M_n       object i (from 1 to n in the file) has probability M_i / Z
H                   Shannon entropy, 5 decimals
```

All integers are decimal and of arbitrary size.

### 2.2 Generation

`scripts/generate_distributions.py` follows the protocol of Saad et al.
(Figure 4 of the FLDR paper), adapted to arbitrary-precision total weights:

1. **Dirichlet parameters spread over the entropy range.** 1000 parameters
   `α` are taken on a logarithmic grid in `[10^-7, 5]`, one Dirichlet
   distribution of dimension `n` is drawn for each (NumPy `RandomState(1)`),
   and, for 1000 target entropies evenly spaced in `[0.001, log2 n]`, the `α`
   whose distribution has the nearest entropy is kept.
2. **Integer weights.** For one parameter out of five (200 distributions), a
   probability vector `p` is drawn from Dirichlet(`α`) and turned into integer
   weights `1 + Multinomial(Z − n, p)` (`RandomState(2)`). The multinomial is
   drawn as a sequence of binomials, with a Gaussian approximation when the
   number of trials exceeds 2^63, so that `Z` can have any size.
3. **Uniform distribution.** The uniform distribution is added, with weights
   computed by `normalize_to_Z`.
4. The 201 distributions are sorted by entropy and written as `d.00000.dist`
   (entropy 0) to `d.00200.dist` (uniform, entropy `log2 n ≈ 14.32`).

Parameters of the paper: `n = 20484`,
`Z = 467734548828719531605240429605504136894484614517198` (about `4.7·10^50`,
the order of magnitude of the combinatorial-chemistry case study), seeds 1 and
2, one parameter out of five. The script is deterministic: it reproduces the
files of the paper byte for byte (§8).

**Note on the uniform distribution.** `normalize_to_Z` rounds `(Z − n)/n` in
double precision; for this `Z`, the sum of the weights exceeds `Z` by about
`1.25·10^34` (`2.7·10^-17` in relative terms). The alias samplers and ALDR
normalize by the sum of the weights; FLDR uses the probabilities `M_i / Z` and
rejects the remaining mass. All of them sample the uniform distribution
exactly.

---

## 3. Samplers

All samplers are exact: the probability of every object is exactly
`M_i / Σ M_j` (FLDR: `M_i / Z`, see §2.2). They work on GMP integers during the
construction, and store their data structure in compact machine words.

| Sampler | Construction | Sampling step | Memory |
|---|---|---|---|
| FLDR (`rej.enc_opt`) | DDG tree of the dyadic rejection distribution, depth `k = ⌈log2 Z⌉` | walk of the tree, one random bit per level | `O(n k)` cells of 32 bits |
| ALDR (`aldr_opt`) | flattened DDG tree of depth `K = 2k` | walk of the tree, level by level | `O(n K)` cells of 32 bits |
| Rust-Alias (`alias.rust_opt`) | alias table of `weighted_rand` (weights scaled by `n`, cells of capacity `W = Σ w_i`) | uniform cell, uniform integer in `[0, W)`, comparison | `n` thresholds of the width of `W` |
| Fraction-Alias (`alias.fractions_opt`) | Walker–Vose on exact rationals | uniform cell, Bernoulli on a rational | numerator and denominator per cell |
| Integer-Alias alea (`alias.integers_opt`) | Algorithm 1 of the paper (cell size `cs = ⌊w/n⌋`, virtual object) | Fast Dice Roller for the cell, Lumbroso's Bernoulli on 31 precomputed digits | 16-byte cells + exact thresholds |
| Integer-Alias speed (`alias.integers_v2`) | Algorithm 1, then 32 digits per threshold | one 64-bit word: Lemire's multiply-shift for the cell, comparison of 32 bits with the digits | same as alea |

### 3.1 FLDR

The construction is a C port of the Python construction of Saad et al.
(`get_rejection_Ms_k`, `make_ddg_matrix`, `make_ddg_tree`, `pack_tree`). With
`g = gcd(Z, M_i)`, the rows are `M_i / g` plus a rejection row
`2^k − Σ M_i / g`, so that the rows sum to `2^k`; they are halved while they
are all even. Since the matrix is dyadic, the DDG tree has no back edge and is
described level by level: at depth `d`, the children of the internal nodes of
depth `d − 1` are, from left to right, the internal nodes, then the leaves of
the rows whose bit `k − d` is 1, by decreasing row. The tree is written in
preorder, as `pack_tree` does: a leaf is one cell `−label`, an internal node two
cells (positions of its children). The walk reads one bit per level inside the
current word of the random-bit stream, and restarts from the root on the
rejection leaf.

FLDR returns the objects numbered from 1 (leaf labels), the other samplers from
0; this only matters for `check_samplers`.

### 3.2 ALDR

Same flattened tree as the ALDR reference implementation (depth `K = 2k`,
factor `c = ⌊2^K / m⌋`, rows `c·a_i` and `2^K mod m`), computed with one GMP
multiplication per object. Level widths and leaves are stored as 32-bit
integers; the walk counters fit in machine integers because a DDG tree has at
most `n + 1` internal nodes per level.

### 3.3 Rust-Alias

`alias_rust_gmp.c` is a C port, to GMP integers, of the table construction of
the Rust crate `weighted_rand`: weights multiplied by `n`, cells of capacity
`W = Σ w_i`, LIFO stacks, objects of weight exactly `W` classified as large.
`alias_rust_opt.c` stores the thresholds as fixed-width limbs and samples as
the crate does: a uniform cell, then a uniform integer in `[0, W)` drawn with
the Fast Dice Roller on preallocated GMP integers (about 185 bits per sample
here), compared with the threshold.

### 3.4 Fraction-Alias

Walker–Vose construction on exact GMP rationals (`mpq_t`). The optimized
storage keeps, for each cell, the numerator and denominator of its threshold as
bytes, and samples the Bernoulli with Lumbroso's algorithm on 31 precomputed
binary digits, falling back to the exact fraction with probability `2^-31`.

### 3.5 Integer-Alias

`alias_integers_opt.c` implements Algorithm 1 of the paper: `cs = ⌊w/n⌋`, a
virtual object of weight `cs − (w mod cs)` when needed, and stacks `L`
(weight `> cs`) and `S` (weight `≤ cs`). Each cell is 16 bytes (digits,
object, alias, position of the exact threshold); the exact thresholds are
stored as bytes, for the cells with an alias only. An attempt that yields the
virtual object is rejected.

* **alea (v1)**: Fast Dice Roller (Lumbroso) for the cell, and Lumbroso's
  Bernoulli on the 31 precomputed digits of `T_i / cs` (first coin equal to 1
  found with `ctz`), with an exact GMP continuation with probability `2^-31`.
  It reads exactly the coins of Lumbroso's algorithms, and thus consumes
  exactly the number of bits analyzed in the paper.
* **speed (v2)**: one 64-bit word per attempt. The high 32 bits give the cell
  by Lemire's multiply-shift (exact, with a rejection of probability
  `< q / 2^32`); the low 32 bits are compared with `D_i = ⌊T_i 2^32 / cs⌋`,
  and equality (probability `2^-32`) is decided exactly with GMP. 64 bits per
  attempt.

This code requires `w ≥ n` (so that `cs ≥ 1`), which holds for the
distributions of the experiments; the stand-alone library of the sampler also
handles `w < n`.

---

## 4. Random-bit stream

`flip_opt.{h,c}`: a background thread fills a ring buffer of 32768 words of 64
bits (256 KiB) with `arc4random_buf` (macOS), by blocks of 512 words, and sleeps
when the buffer is full. The sampling thread reads bits from the current word
with inline functions (`fo_bits`, `fo_bit`, `fo_coins_until_one`, `fo_uniform`),
least significant bit first, and synchronizes with the producer only when it
has consumed the words it knows of (release / acquire on two counters placed on
separate cache lines).

The counter `NUM_RNG_CALLS` is incremented by 2 for each 64-bit word consumed,
so it counts 32-bit words; the number of random bits per sample is
`32 · NUM_RNG_CALLS / samples`.

On Linux, `arc4random_buf` is available in glibc 2.36 and later.

---

## 5. Measurement programs

### 5.1 `measure_runtime SEED STEPS SAMPLER FILE.dist`

1. Writes a 128 MB buffer, so that every run starts with cold caches.
2. Starts the random-bit stream.
3. Reads `FILE.dist` and builds the data structure (not timed).
4. Runs `STEPS` samples, timed with `clock()`.

`clock()` measures the CPU time of the process, including the producer thread;
the time per sample is `seconds / STEPS`. The output's last line is
`SAMPLER SECONDS WORDS REJECTIONS`.

### 5.2 `measure_sizes SAMPLER FILE.dist OUT`

Writes the memory of the data structure, in bytes: the size of the C structure
plus the useful part of its arrays and integers (`get_size_*` functions). For
FLDR, 4 bytes per encoding cell; for ALDR, 4 bytes per level width and per
leaf; for Rust-Alias, 4 bytes per alias and the fixed-width thresholds; for
Fraction-Alias, 16 bytes per cell and the bytes of the fractions; for
Integer-Alias, 16 bytes per cell and the bytes of the exact thresholds.

### 5.3 `measure_preprocess SAMPLER FILE.dist`

Reads the weights (not timed), then repeats the construction of the data
structure from the weights (at least 5 times and 50 ms) and prints the
smallest time, in seconds (`CLOCK_MONOTONIC_RAW`).

### 5.4 `check_samplers FILE.dist [SAMPLES]`

Draws `SAMPLES` objects (default `10^6`) with each sampler and compares the
counts with the exact probabilities by a chi-square test (expected counts below
5 pooled, Wilson–Hilferty approximation, failure for `z > 4.5`).

---

## 6. Pipeline

`scripts/pipeline.sh` runs the measurements **sequentially**: one process at a
time, each sampler once per distribution, the six samplers in a row for each
distribution. Parallel measurements, as in the original pipeline of Saad et
al., compete for cores, caches and memory bandwidth, and inflated the
measured times by 15 to 21 % (medians) in our tests.

| Command | Action | Output |
|---|---|---|
| `generate` | §2 | `results/distributions/d.XXXXX.dist` |
| `runtimes [STEPS]` | `measure_runtime` on every distribution and sampler (default `10^7` samples) | `results/raw/d.XXXXX.SAMPLER.runtime` |
| `sizes` | `measure_sizes` | `results/measures/SAMPLER.sizes` |
| `preprocess` | `measure_preprocess` | `results/measures/SAMPLER.preprocess` |
| `aggregate [STEPS]` | gathers the raw outputs | `results/measures/SAMPLER.{runtimes,calls,bits,rejet}` |
| `figures` | §7 | `results/figures/` |
| `all [STEPS]` | all of the above | |

Each file of `results/measures/` has one line per distribution, in the order of
the `.dist` files (increasing entropy). `runtimes` holds the seconds for
`STEPS` samples, `bits` the random bits per sample, `sizes` bytes,
`preprocess` seconds. The script sets `LC_ALL=C`, so that numbers are written
with a decimal point whatever the locale.

---

## 7. Figures

`scripts/make_figures.py` reads `results/distributions/` (entropies) and
`results/measures/`, and draws the six panels of Figures 1 and 2 of the paper
with the same scales, labels and colors (`Runtime-memory`, `Entropy-runtime`,
`Entropy-memory`, `Entropy-preprocessing`, `Entropy-bits`, `Entropy-bits2`), as
PDF and PNG. It uses LaTeX for the text when available (`--no-tex` to disable).

---

## 8. Validation

* **Distributions**: the 201 files written by `generate_distributions.py` are
  identical, byte for byte, to those of the paper.
* **FLDR**: the encoding built by the C port is identical to the one written by
  the original Python code, for the 201 distributions.
* **Integer-Alias**: the tables are identical, cell by cell, to those of the
  stand-alone library of the sampler (both modes).
* **Memory**: the sizes measured by this repository are identical to those of
  the paper.
* **Exactness**: `check_samplers` passes for the six samplers on the
  distributions of the experiments.

---

## 9. Differences with the original experiments

This repository keeps only what the paper needs. Compared with the repository
the experiments were developed in (itself derived from the experiments of
Saad et al.):

* only the six samplers of the paper (the floating-point samplers of Saad et
  al. and the non-optimized GMP versions are not included);
* one input file per distribution (`.dist`): every sampler builds its data
  structure in C. In particular, the FLDR tree, previously built in Python and
  read from a separate file, is built by the C port of §3.1, which gives the
  same encoding;
* sequential pipeline only, and a Python script instead of notebooks for the
  figures;
* in `Entropy-bits2`, the legend leaves out the two samplers outside the
  vertical range (Rust-Alias and Integer-Alias speed), as intended by the
  original notebook.
