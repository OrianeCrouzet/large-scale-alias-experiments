#!/usr/bin/env python3
# Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
# SPDX-License-Identifier: GPL-3.0-or-later
"""Draws the figures of the paper from the measurements of scripts/pipeline.sh.

    python3 scripts/make_figures.py [--results DIR] [--no-tex]

Reads DIR/distributions/*.dist (entropy on the last line) and
DIR/measures/SAMPLER.{runtimes,sizes,bits,preprocess}, and writes, in
DIR/figures/, as PDF and PNG:

    Runtime-memory          sampling time against memory (Fig. 1, left)
    Entropy-runtime         sampling time against entropy (Fig. 1, right)
    Entropy-memory          memory against entropy (Fig. 2, top left)
    Entropy-preprocessing   preprocessing time against entropy (Fig. 2, top right)
    Entropy-bits            random bits per sample against entropy (Fig. 2, bottom left)
    Entropy-bits2           the same, vertical axis capped at 30 bits (Fig. 2, bottom right)

LaTeX is used for the text when it is available (as in the paper), unless
--no-tex is given.
"""

import argparse
import glob
import os
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

STEPS = 10_000_000          # samples per measurement (scripts/pipeline.sh)

# Plotting order, labels and colors of the paper.
SELECTION = [
    'alias.integers_opt',
    'alias.integers_v2',
    'rej.enc_opt',
    'aldr_opt',
    'alias.rust_opt',
    'alias.fractions_opt',
]
LABELS = {
    'alias.integers_opt': 'Our Integer-Alias Sampler (alea opt.)',
    'alias.integers_v2': 'Our Integer-Alias Sampler (speed opt.)',
    'rej.enc_opt': 'FLDR',
    'aldr_opt': 'ALDR',
    'alias.rust_opt': 'Rust-Alias Sampler',
    'alias.fractions_opt': 'Fraction-Alias Sampler',
}
_TAB20 = plt.get_cmap('tab20').colors
COLORS = {
    'rej.enc_opt': _TAB20[0],
    'aldr_opt': _TAB20[2],
    'alias.rust_opt': _TAB20[4],
    'alias.fractions_opt': _TAB20[6],
    'alias.integers_opt': _TAB20[8],
    'alias.integers_v2': _TAB20[18],
}
# Out of the range of Entropy-bits2: Rust-Alias (about 185 bits) and v2 (64 bits).
OUT_OF_ZOOM = {'alias.rust_opt', 'alias.integers_v2'}


def load(results):
    dists = sorted(glob.glob(os.path.join(results, 'distributions', 'd.*.dist')))
    if not dists:
        raise SystemExit('no distribution in %s/distributions' % results)
    entropies = []
    for path in dists:
        with open(path) as f:
            entropies.append(float(f.read().split()[-1]))
    n = int(open(dists[0]).read().split()[1])
    data = {}
    for s in SELECTION:
        d = {}
        for q in ('runtimes', 'sizes', 'bits', 'preprocess'):
            path = os.path.join(results, 'measures', '%s.%s' % (s, q))
            if os.path.exists(path):
                d[q] = np.loadtxt(path, ndmin=1)
                if len(d[q]) != len(dists):
                    raise SystemExit('%s: %d values for %d distributions' % (path, len(d[q]), len(dists)))
        if d:
            data[s] = d
    return np.asarray(entropies), n, data


def annotate(ax, n, legend_keep=None, labelspacing=0.3, markerscale=None):
    ax.text(0.97, 0.03, r'$n=%d, w\approx10^{50}$' % n, transform=ax.transAxes,
            ha='right', va='bottom', fontsize=11,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    handles, labels = ax.get_legend_handles_labels()
    keep = [i for i in range(len(handles)) if legend_keep is None or legend_keep(labels[i])]
    kwargs = dict(loc='upper left', framealpha=0, ncol=2, fontsize=7,
                  labelspacing=labelspacing, handletextpad=.5, handlelength=.5)
    if markerscale is not None:
        kwargs['markerscale'] = markerscale
    ax.legend([handles[i] for i in keep], [labels[i] for i in keep], **kwargs)


def save(fig, out, name):
    fig.set_tight_layout(True)
    fig.set_size_inches(5.25, 3.5)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(out, '%s.%s' % (name, ext)))
    plt.close(fig)
    print(os.path.join(out, name + '.{pdf,png}'))


def against_entropy(entropies, data, quantity, scale, legend_size):
    fig, ax = plt.subplots()
    for s in SELECTION:
        if s not in data or quantity not in data[s]:
            continue
        y = data[s][quantity] / STEPS if quantity == 'runtimes' else data[s][quantity]
        sc = ax.scatter([], [], label=LABELS[s], color=COLORS[s], s=legend_size)
        ax.scatter(entropies, y, color=sc.get_facecolors()[0], marker='.', s=2)
    if scale == 'log2':
        ax.set_yscale('log', base=2)
    ax.set_xlabel('Entropy of Target Distribution', fontsize=11 if legend_size == 10 else 14)
    return fig, ax


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--results', default=os.path.join(os.path.dirname(__file__), '..', 'results'))
    parser.add_argument('--no-tex', action='store_true', help='do not use LaTeX for the text')
    args = parser.parse_args()

    if not args.no_tex and shutil.which('latex'):
        plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'sans-serif'

    entropies, n, data = load(args.results)
    out = os.path.join(args.results, 'figures')
    os.makedirs(out, exist_ok=True)

    # Runtime-memory
    fig, ax = plt.subplots()
    for s in SELECTION:
        if s in data and 'runtimes' in data[s] and 'sizes' in data[s]:
            ax.scatter(data[s]['runtimes'] / STEPS, data[s]['sizes'], label=LABELS[s],
                       color=COLORS[s], s=10, marker='.')
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)
    ax.set_xlabel('Sampler Runtime (Seconds per Sample)', fontsize=14)
    ax.set_ylabel('Sampler Memory (Bytes)', fontsize=14)
    ax.set_xlim([.9 * 2**-28, 2**-23])
    ax.set_ylim([2**17, 2**27])
    ax.xaxis.set_major_locator(ticker.LogLocator(base=2., subs=(1,), numticks=7))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.LogLocator(base=2., subs=(1,), numticks=11))
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    annotate(ax, n, labelspacing=0.2, markerscale=4)
    save(fig, out, 'Runtime-memory')

    # Entropy-runtime
    fig, ax = against_entropy(entropies, data, 'runtimes', 'log2', 10)
    ax.set_ylabel('Sampler Runtime (Seconds per Sample)', fontsize=11)
    ax.set_ylim([2**-29, 2**-22])
    annotate(ax, n, markerscale=2)
    save(fig, out, 'Entropy-runtime')

    # Entropy-memory
    fig, ax = against_entropy(entropies, data, 'sizes', 'log2', 50)
    ax.set_ylabel('Sampler Memory (Bytes)', fontsize=12)
    ax.set_ylim([2**17, 2**27])
    annotate(ax, n)
    save(fig, out, 'Entropy-memory')

    # Entropy-preprocessing
    fig, ax = against_entropy(entropies, data, 'preprocess', 'log2', 10)
    ax.set_ylabel('Preprocessing Time (Seconds)', fontsize=11)
    ax.set_ylim([2**-15, 2**-2])
    annotate(ax, n, markerscale=2)
    save(fig, out, 'Entropy-preprocessing')

    # Entropy-bits and Entropy-bits2
    fig, ax = against_entropy(entropies, data, 'bits', 'linear', 50)
    ax.set_ylabel('Average number of bits consumed', fontsize=12)
    ax.set_ylim([-20, 240])
    annotate(ax, n)
    save(fig, out, 'Entropy-bits')

    fig, ax = against_entropy(entropies, data, 'bits', 'linear', 50)
    ax.set_ylabel('Average number of bits consumed', fontsize=12)
    ax.set_ylim([0, 30])
    out_labels = {LABELS[s] for s in OUT_OF_ZOOM}
    annotate(ax, n, legend_keep=lambda label: label not in out_labels)
    save(fig, out, 'Entropy-bits2')


if __name__ == '__main__':
    main()
