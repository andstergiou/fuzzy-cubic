#!/bin/bash
# Run the full DMRG workflow locally for one (nm, w, h) parameter set.
#
# Usage:
#   ./run_dmrg.sh <nm> <w> [h] [n_states] [threads]
#
# Arguments:
#   nm        : system size (e.g. 12)
#   w         : cubic anisotropy coupling (e.g. 0.6)
#   h         : polarization parameter (e.g. 15.323); if omitted, the tabulated
#               optimum h_opt(nm, w) from ../hopt.jl is used
#   n_states  : number of states per sector (default: 8)
#   threads   : Julia threads to use (default: number of CPU cores)
#
# Examples:
#   ./run_dmrg.sh 12 0.6 15.323
#   ./run_dmrg.sh 12 0.6            # same thing, h looked up automatically
#
# The script runs the three steps in sequence:
#   1. Build MPO cache (skipped if cache already exists)
#   2. Run DMRG for each of the 4 symmetry sectors (Z=0/1, P2=0/1) sequentially
#   3. Combine sector results into a single spectrum file

set -e

if [[ $# -lt 2 ]]; then
    echo "Usage: ./run_dmrg.sh <nm> <w> [h] [n_states] [threads]"
    exit 1
fi

NM=$1
W=$2
H=$3
N_STATES=${4:-8}
THREADS=${5:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# If h was not given, resolve the tabulated optimum for this (nm, w). Doing it
# here rather than in each script keeps every step, including the combine step,
# working on exactly the same value.
if [[ -z "$H" ]]; then
    H=$(julia -e \
        'include(ARGS[1]); print(default_h(parse(Int, ARGS[2]), parse(Float64, ARGS[3])))' \
        "$SCRIPT_DIR/../hopt.jl" "$NM" "$W")
    if [[ -z "$H" ]]; then
        echo "ERROR: could not determine a default h for nm=$NM, w=$W; pass it explicitly."
        exit 1
    fi
    echo "h not supplied; using tabulated h_opt($NM, $W) = $H"
fi

echo "=============================================="
echo "Cubic DMRG — local run"
echo "  nm=$NM  w=$W  h=$H  n_states=$N_STATES  threads=$THREADS"
echo "=============================================="

# ── Step 1: Build MPO cache ────────────────────────────────────────────────────
# The build script derives the cache filename from (nm, w, h) itself and exits
# early if the cache already exists, so there is no filename formatting here.
echo ""
echo "Step 1/3: Building MPO cache (skipped if it already exists)..."
JULIA_NUM_THREADS=$THREADS julia --project=. -t $THREADS \
    cubic_spectrum_dmrg_build_mpo.jl $NM $W $H

# ── Step 2: Run DMRG for each sector ──────────────────────────────────────────
echo ""
echo "Step 2/3: Running DMRG sectors (Z=0/1, P2=0/1) sequentially..."
echo "  (This is the slow step — bond dimensions up to 1500+ depending on nm)"

for Z in 0 1; do
    for P2 in 0 1; do
        echo ""
        echo "  Sector Z=$Z P2=$P2..."
        JULIA_NUM_THREADS=$THREADS julia --project=. -t $THREADS \
            cubic_spectrum_dmrg.jl $Z $P2 $N_STATES $NM $W $H --resume
    done
done

# ── Step 3: Combine results ────────────────────────────────────────────────────
echo ""
echo "Step 3/3: Combining sector results..."
julia --project=. cubic_spectrum_dmrg_combine.jl $NM $H $W

echo ""
echo "=============================================="
echo "Done. See the combined output file printed above."
echo "=============================================="
