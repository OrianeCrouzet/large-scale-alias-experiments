# Alias_experiments

Code and scripts to reproduce the experiments of the paper

> Oriane Crouzet, Antoine Genitrini, Mehdi Naima.
> *Exact Integer Alias Method: Robust and Error-Free Large-Scale Sampling*.
> 2026. Link to the archived version: to be added.

The experiments compare six **exact** samplers on distributions with
arbitrary-precision integer weights (`n = 20 484` objects, total weight
`w ≈ 4.7·10^50`): sampling time, memory of the data structure, preprocessing
time and number of random bits per sample, as functions of the entropy of the
distribution.

| Sampler | Name in the code | Origin |
|---|---|---|
| FLDR (Fast Loaded Dice Roller) | `rej.enc_opt` | Saad et al., AISTATS 2020 |
| ALDR (Amplified Loaded Dice Roller) | `aldr_opt` | Draper and Saad |
| Rust-Alias | `alias.rust_opt` | Rust crate `weighted_rand` |
| Fraction-Alias | `alias.fractions_opt` | Walker–Vose on exact rationals |
| Integer-Alias, *alea* (v1, fewest random bits) | `alias.integers_opt` | this paper |
| Integer-Alias, *speed* (v2, fastest) | `alias.integers_v2` | this paper |

All samplers are written in C with GMP and read random bits from the same
cryptographic stream. The stand-alone library of the Integer-Alias Sampler is
distributed separately; this repository keeps the exact code used for the
figures of the paper.

## Contents

```
c/                       samplers and measurement programs (C, GMP)
scripts/
  generate_distributions.py   the 201 distributions of the paper (.dist files)
  pipeline.sh                 sequential measurement pipeline
  make_figures.py             the six figures of the paper
README.md  TECHNICAL.md  NOTICE.md  LICENSE  LICENSES/
```

[TECHNICAL.md](TECHNICAL.md) describes the samplers, the measurements and the
pipeline in detail; [NOTICE.md](NOTICE.md) lists the third-party code.

## Requirements

* a C11 compiler (GCC or Clang) and [GMP](https://gmplib.org/);
* Python 3 with NumPy and Matplotlib (`pip install -r requirements.txt`);
* optionally LaTeX, for the fonts of the figures (as in the paper); without
  LaTeX, the figures use Matplotlib's default fonts.

The results of the paper were obtained on a MacBook Pro (Apple M4 Pro, 14-core
CPU, 48 GB, macOS 26.6), with Apple Clang 21, GMP 6.3.0, Python 3.13,
NumPy 2.3.4 and Matplotlib 3.10.

## Reproducing the experiments

Everything at once (about 15 minutes on the machine above):

```sh
scripts/pipeline.sh all
```

Step by step:

```sh
make -C c                          # builds the programs (done automatically by the pipeline)
scripts/pipeline.sh generate       # 201 distributions in results/distributions/ (~15 s)
scripts/pipeline.sh runtimes       # sampling time, 10^7 samples per sampler and distribution (~10 min)
scripts/pipeline.sh sizes          # memory of the data structures (~1 min)
scripts/pipeline.sh preprocess     # preprocessing time (~2 min)
scripts/pipeline.sh aggregate      # one file per sampler and quantity in results/measures/
scripts/pipeline.sh figures        # results/figures/*.{pdf,png}
```

The measurements run one at a time: running them in parallel inflates the
times (shared cores, caches and memory bandwidth). For stable times, keep the
machine otherwise idle and on mains power; on macOS, `caffeinate -i
scripts/pipeline.sh all` prevents sleep. Set `RESULTS=dir` to write the results
elsewhere.

**Distributions.** `generate` is deterministic: it writes the 201 `.dist` files
used for the paper, identical byte for byte.

**Figures.** `results/figures/` then contains:

| File | Figure of the paper |
|---|---|
| `Runtime-memory` | Fig. 1, left: sampling time against memory |
| `Entropy-runtime` | Fig. 1, right: sampling time against entropy |
| `Entropy-memory` | Fig. 2, top left: memory against entropy |
| `Entropy-preprocessing` | Fig. 2, top right: preprocessing time against entropy |
| `Entropy-bits` | Fig. 2, bottom left: random bits per sample |
| `Entropy-bits2` | Fig. 2, bottom right: the same, capped at 30 bits |

Memory and random bits are reproduced exactly. Times depend on the machine,
and vary by up to about 30 % between runs on the same machine; the relative
positions of the samplers are stable.

**Exactness check.** `c/check_samplers FILE.dist [SAMPLES]` draws from each
sampler and runs a chi-square test against the exact probabilities, for
example:

```sh
c/check_samplers results/distributions/d.00100.dist 3000000
```

## License

Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. See [LICENSE](LICENSE).

It includes code derived from the experiments of the Fast Loaded Dice Roller
(F. A. Saad et al., Apache License 2.0), from the Amplified Loaded Dice Roller
(CMU Probabilistic Computing Systems Lab, Apache License 2.0) and from the Rust
crate `weighted_rand` (Ippee, MIT License): see [NOTICE.md](NOTICE.md) and
[LICENSES/](LICENSES).
