#!/usr/bin/env python3
# Derived from experiments/dists.py, experiments/create_distributions.py,
# src/entropy.py and src/utils.py of the fast-loaded-dice-roller-experiments
# repository of F. A. Saad et al.
# (https://github.com/probsys/fast-loaded-dice-roller-experiments),
# Copyright 2019 MIT Probabilistic Computing Project, released under the
# Apache License 2.0 (see LICENSES/Apache-2.0.txt).
#
# Modified in 2026 by Oriane Crouzet, Antoine Genitrini and Mehdi Naima:
# arbitrary-precision total weight (multinomial with big integers), uniform
# distribution added, only the .dist files are written; gathered in one file.
# This file is distributed under the GNU General Public License v3.0 or later
# as part of Alias_experiments (see LICENSE).
# SPDX-License-Identifier: GPL-3.0-or-later
"""Generates the random distributions of the experiments (.dist files).

With the default parameters, writes the 201 distributions of the paper:
200 Dirichlet-multinomial distributions of n = 20484 objects with total weight
Z = 467734548828719531605240429605504136894484614517198 (about 4.7e50),
whose entropies cover [0, log2 n], plus the uniform distribution; sorted by
entropy, in DIR/d.00000.dist ... DIR/d.00200.dist.

The output is deterministic (NumPy RandomState, seeds 1 and 2): the files are
identical, byte for byte, to those used for the figures of the paper.

Format of a .dist file:
    Z                   common denominator of the probabilities
    n M_1 ... M_n       numerators: object i has probability M_i / Z
    H                   entropy (5 decimals)
Note: for the uniform distribution, the rounding of normalize_to_Z makes
sum(M_i) slightly larger than Z (by about 2.7e-17 in relative terms); the
samplers normalize by sum(M_i), FLDR by Z (with rejection), and all of them
sample the uniform distribution exactly.
"""

import argparse
import math
import multiprocessing
import os
from fractions import Fraction
from math import isinf, isnan, log2

import numpy
from numpy import lcm


def compute_entropy(ps):
    """Binary entropy of the probability vector ps."""
    return sum(-log2(p) * p for p in ps if p != 0)


def get_nearest_neighbors(k, values, x):
    """Indexes of the k items of values nearest to x."""
    values_idx = list(enumerate(values))
    values_sorted = sorted(values_idx, key=lambda v: abs(v[1] - x))
    return [r[0] for r in values_sorted[:k]]


def sample_dirichlet(alpha, N, rng):
    """Dirichlet sampler, retried on nan or inf."""
    for _ in range(1000):
        try:
            dist = rng.dirichlet([alpha] * N)
            if not (any(map(isnan, dist)) or any(map(isinf, dist))):
                return dist
        except ZeroDivisionError:
            continue
    raise AssertionError('Failed to generate dirichlet distribution.')


def get_alpha_entropies(n, maxalpha=5, numalpha=1000, parallel=None):
    """Dirichlet parameters giving entropies spread over [0, log2 n]."""
    rng = numpy.random.RandomState(1)
    low = -5 if n < 10000 else -7
    high = numpy.log10(maxalpha)
    alphas = numpy.logspace(low, high, numalpha)
    dists = [sample_dirichlet(alpha, n, rng) for alpha in alphas]
    if parallel:
        with multiprocessing.Pool() as pool:
            entropies_list = pool.map(compute_entropy, dists)
    else:
        entropies_list = list(map(compute_entropy, dists))
    entropies = numpy.asarray(entropies_list)
    entropies[numpy.isnan(entropies)] = 1e-10
    entropies_desired = numpy.linspace(0.001, log2(n), numalpha)
    idxs = [get_nearest_neighbors(1, entropies, e) for e in entropies_desired]
    return [alphas[i[0]] for i in idxs]


def sample_large_binomial(n, p, rng):
    """Binomial(n, p) for an arbitrarily large n."""
    if n <= 0 or p <= 0:
        return 0
    if p >= 1:
        return n
    if n <= 9223372036854775807:
        return int(rng.binomial(n, p))
    # Very large n: Gaussian approximation.
    mu = float(n) * p
    sigma = math.sqrt(float(n) * p * (1.0 - p))
    z = rng.standard_normal()
    val = round(mu + sigma * z)
    return max(0, min(int(val), n))


def sample_large_multinomial(total, probabilities, rng):
    """Multinomial(total, probabilities) for an arbitrarily large total."""
    remaining_total = total
    remaining_p_sum = float(sum(probabilities))
    counts = []
    for p in probabilities[:-1]:
        if remaining_total <= 0:
            counts.append(0)
            continue
        q = min(1.0, max(0.0, float(p) / remaining_p_sum))
        count = sample_large_binomial(remaining_total, q, rng)
        counts.append(count)
        remaining_total -= count
        remaining_p_sum -= float(p)
    counts.append(max(0, remaining_total))
    return counts


def sample_dirichlet_multinomial_positive(alpha, N, Z, rng):
    """Positive Dirichlet-multinomial distribution with integer weights summing to Z."""
    assert N <= Z
    probabilities = sample_dirichlet(alpha, N, rng)
    numerators = sample_large_multinomial(Z - N, probabilities, rng)
    numerators = [n + 1 for n in numerators]
    return [Fraction(n, Z) for n in numerators]


def normalize_to_Z(weights, Z):
    """Integer weights proportional to weights, each >= 1, summing to about Z."""
    weights = numpy.maximum(weights, 0)
    w_sum = numpy.sum(weights)
    if w_sum == 0:
        weights = numpy.ones(len(weights))
        w_sum = len(weights)
    n = len(weights)
    if Z < n:
        raise ValueError('Z must be at least n to assign at least 1 to each element')
    probs = weights / w_sum
    scaled = probs * float(Z - n)
    floored = numpy.floor(scaled)
    integers = [int(x) + 1 for x in floored]
    diff = int(Z - sum(integers))
    if diff > 0:
        frac = scaled - floored
        idx = numpy.argsort(-frac)
        for i in range(diff):
            integers[idx[i]] += 1
    return integers


def get_common_denominator(probabilities):
    """Least Z such that each probability is a multiple of 1/Z."""
    return lcm.reduce([p.denominator for p in probabilities])


def write_dist(dirname, idx, p_target, entropy):
    Z = get_common_denominator(p_target)
    Ms = [int(Z * p) for p in p_target]
    path = os.path.join(dirname, 'd.%05d.dist' % (idx,))
    with open(path, 'w') as f:
        f.write('%d\n' % (Z,))
        f.write('%d %s\n' % (len(p_target), ' '.join(map(str, Ms))))
        f.write('%1.5f\n' % (entropy,))
    return path


def _write(args):
    return write_dist(*args)


def generate(N, Z, seed, thin, offset, dirname):
    os.makedirs(dirname, exist_ok=True)
    rng = numpy.random.RandomState(seed)
    alphas = get_alpha_entropies(N, maxalpha=5, numalpha=1000, parallel=True)
    distributions = [sample_dirichlet_multinomial_positive(a, N, Z, rng)
                     for a in alphas[offset::thin]]
    # Last distribution: the uniform one.
    weights = normalize_to_Z(numpy.ones(N), Z)
    distributions.append([Fraction(int(w), Z) for w in weights])
    with multiprocessing.Pool() as pool:
        entropies = pool.map(compute_entropy, distributions)
    idxs = numpy.argsort(entropies)
    args = [(dirname, i, distributions[idx], entropies[idx]) for i, idx in enumerate(idxs)]
    with multiprocessing.Pool() as pool:
        for path in pool.map(_write, args):
            print(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=20484)
    parser.add_argument('--Z', type=int, default=467734548828719531605240429605504136894484614517198)
    parser.add_argument('--seed', type=int, default=2)
    parser.add_argument('--thin', type=int, default=5)
    parser.add_argument('--offset', type=int, default=0)
    parser.add_argument('--out', default='distributions')
    a = parser.parse_args()
    generate(a.n, a.Z, a.seed, a.thin, a.offset, a.out)
