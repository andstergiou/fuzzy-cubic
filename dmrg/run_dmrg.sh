#!/bin/bash
# Run the full DMRG workflow locally for one (nm, w, h) parameter set.
#
# Usage:
#   ./run_dmrg.sh <nm> <w> <h> [n_states] [threads]
#
# Arguments:
#   nm        : system size (e.g. 12)
#   w         : cubic anisotropy coupling (e.g. 0.6)
#   h         : polarization parameter (e.g. 15.323)
#   n_states  : number of states per sector (default: 8)
#   threads   : Julia threads to use (default: number of CPU cores)
#
# Example:
#   ./run_dmrg.sh 12 0.6 15.323
#
# The script runs the three steps in sequence:
#   1. Build MPO cache (skipped if cache already exists)
#   2. Run DMRG for each of the 4 symmetry sectors (Z=0/1, P2=0/1) sequentially
#   3. Combine sector results into a single spectrum file

set -e

if [[ $# -lt 3 ]]; then
    echo "Usage: ./run_dmrg.sh <nm> <w> <h> [n_states] [threads]"
    exit 1
fi

NM=$1
W=$2
H=$3
N_STATES=${4:-8}
THREADS=${5:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}

# Format h and w for filename matching (replace . with p)
H_STR=$(echo "$H" | sed 's/\./p/g' | sed 's/p\([0-9]\{3\}\).*/p\1/')
W_STR=$(echo "$W" | sed 's/\./p/g')
MPO_CACHE="mpo_cache_nm${NM}_w${W_STR}_h${H_STR}.jls"

echo "=============================================="
echo "Cubic DMRG — local run"
echo "  nm=$NM  w=$W  h=$H  n_states=$N_STATES  threads=$THREADS"
echo "=============================================="

# ── Step 1: Build MPO cache ────────────────────────────────────────────────────
if [[ -f "$MPO_CACHE" ]]; then
    echo ""
    echo "Step 1/3: MPO cache already exists ($MPO_CACHE). Skipping build."
else
    echo ""
    echo "Step 1/3: Building MPO cache..."
    JULIA_NUM_THREADS=$THREADS julia --project=. -t $THREADS cubic_spectrum_dmrg_build_mpo.jl $NM
fi

# ── Step 2: Run DMRG for each sector ──────────────────────────────────────────
echo ""
echo "Step 2/3: Running DMRG sectors (Z=0/1, P2=0/1) sequentially..."
echo "  (This is the slow step — bond dimensions up to 1500+ depending on nm)"

for Z in 0 1; do
    for P2 in 0 1; do
        echo ""
        echo "  Sector Z=$Z P2=$P2..."
        JULIA_NUM_THREADS=$THREADS julia --project=. -t $THREADS \
            cubic_spectrum_dmrg.jl $Z $P2 $N_STATES $NM --resume
    done
done

# ── Step 3: Combine results ────────────────────────────────────────────────────
echo ""
echo "Step 3/3: Combining sector results..."
julia --project=. cubic_spectrum_dmrg_combine.jl $NM $H $W

echo ""
echo "=============================================="
echo "Done. Output: dmrg_combined_nm${NM}_w${W_STR}_h${H_STR}.txt"
echo "=============================================="
