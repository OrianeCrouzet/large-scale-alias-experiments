#!/bin/bash
# Copyright (C) 2026 Oriane Crouzet, Antoine Genitrini, Mehdi Naima
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Sequential measurement pipeline of the experiments.
#
# Every measurement runs alone (one process at a time), each sampler once per
# distribution, the samplers interleaved for each distribution, so that a drift
# of the machine affects all samplers alike.
#
#   scripts/pipeline.sh generate            writes the 201 distributions
#   scripts/pipeline.sh runtimes [STEPS]    sampling time (default 10^7 samples)
#   scripts/pipeline.sh sizes               memory of the data structures
#   scripts/pipeline.sh preprocess          preprocessing time
#   scripts/pipeline.sh aggregate [STEPS]   one file per sampler and quantity
#   scripts/pipeline.sh figures             the figures of the paper
#   scripts/pipeline.sh all [STEPS]         everything above
#
# Results go to $RESULTS (default: results/):
#   distributions/d.XXXXX.dist   input distributions
#   raw/d.XXXXX.SAMPLER.runtime  output of measure_runtime
#   measures/SAMPLER.{runtimes,calls,bits,rejet,sizes,preprocess}
#                                one line per distribution, in the order of the .dist files
#   figures/*.{pdf,png}

set -euo pipefail
export LC_ALL=C          # decimal point in every number written by awk and printf

ROOT=$(cd "$(dirname "$0")/.." && pwd)
RESULTS=${RESULTS:-$ROOT/results}
DISTS=$RESULTS/distributions
RAW=$RESULTS/raw
MEASURES=$RESULTS/measures
BIN=$ROOT/c

SAMPLERS='alias.integers_opt alias.integers_v2 rej.enc_opt aldr_opt alias.rust_opt alias.fractions_opt'
DEFAULT_STEPS=10000000

cmd=${1:-}
shift || true

list_dists() {
    ls "$DISTS"/d.*.dist 2>/dev/null || { echo "no distribution in $DISTS (run: $0 generate)" >&2; exit 1; }
}

build() {
    make -s -C "$BIN"
}

case "$cmd" in
generate)
    python3 "$ROOT/scripts/generate_distributions.py" --out "$DISTS" > /dev/null
    echo "$(ls "$DISTS"/d.*.dist | wc -l | tr -d ' ') distributions in $DISTS"
    ;;

runtimes)
    steps=${1:-$DEFAULT_STEPS}
    build
    mkdir -p "$RAW"
    dists=$(list_dists)
    total=$(echo "$dists" | wc -l | tr -d ' ')
    start=$(date +%s)
    k=0
    for d in $dists; do
        k=$((k + 1))
        base=$(basename "$d" .dist)
        for sampler in $SAMPLERS; do
            "$BIN/measure_runtime" 0 "$steps" "$sampler" "$d" > "$RAW/$base.$sampler.runtime"
        done
        elapsed=$(( $(date +%s) - start ))
        printf '[%d/%d] %s  (elapsed %dm%02ds, remaining ~%dm%02ds)\n' "$k" "$total" "$base" \
            $((elapsed / 60)) $((elapsed % 60)) \
            $((elapsed * (total - k) / k / 60)) $((elapsed * (total - k) / k % 60))
    done
    ;;

sizes)
    build
    mkdir -p "$MEASURES"
    tmp=$(mktemp)
    for sampler in $SAMPLERS; do
        rm -f "$MEASURES/$sampler.sizes"
        for d in $(list_dists); do
            "$BIN/measure_sizes" "$sampler" "$d" "$tmp"
            cat "$tmp" >> "$MEASURES/$sampler.sizes"
        done
        echo "wrote $MEASURES/$sampler.sizes"
    done
    rm -f "$tmp"
    ;;

preprocess)
    build
    mkdir -p "$MEASURES"
    for sampler in $SAMPLERS; do
        rm -f "$MEASURES/$sampler.preprocess"
    done
    for d in $(list_dists); do
        for sampler in $SAMPLERS; do
            "$BIN/measure_preprocess" "$sampler" "$d" >> "$MEASURES/$sampler.preprocess"
        done
    done
    echo "wrote $MEASURES/*.preprocess"
    ;;

aggregate)
    steps=${1:-$DEFAULT_STEPS}
    mkdir -p "$MEASURES"
    for sampler in $SAMPLERS; do
        rm -f "$MEASURES/$sampler".{runtimes,calls,bits,rejet}
        for d in $(list_dists); do
            line=$(tail -n1 "$RAW/$(basename "$d" .dist).$sampler.runtime")
            echo "$line" | awk '{print $2}' >> "$MEASURES/$sampler.runtimes"
            echo "$line" | awk '{print $3}' >> "$MEASURES/$sampler.calls"
            echo "$line" | awk '{print $4}' >> "$MEASURES/$sampler.rejet"
            # bits per sample = 32-bit words * 32 / steps
            echo "$line" | awk -v steps="$steps" '{print ($3 * 32) / steps}' >> "$MEASURES/$sampler.bits"
        done
        printf '%-20s mean bits per sample %s\n' "$sampler" \
            "$(awk '{s += $1} END {printf "%.3f", s / NR}' "$MEASURES/$sampler.bits")"
    done
    ;;

figures)
    python3 "$ROOT/scripts/make_figures.py" --results "$RESULTS"
    ;;

all)
    steps=${1:-$DEFAULT_STEPS}
    [ -n "$(ls "$DISTS"/d.*.dist 2>/dev/null)" ] || "$0" generate
    "$0" runtimes "$steps"
    "$0" sizes
    "$0" preprocess
    "$0" aggregate "$steps"
    "$0" figures
    ;;

*)
    sed -n '4,25p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
    ;;
esac
