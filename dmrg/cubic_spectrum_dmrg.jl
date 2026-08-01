# Cubic DMRG — single symmetry sector
#
# Usage:
#   julia cubic_spectrum_dmrg.jl [Z] [P2] [n_states] [nm] [w] [h] [--resume]
#
# Arguments:
#   Z           : Total parity sector (0 or 1), default 0
#   P2          : N2 parity sector (0 or 1), default 0
#   n_states    : Number of states to compute, default 8
#   nm          : System size, default 12
#   w           : Cubic anisotropy coupling, default 0.6
#   h           : Polarization parameter, default h_opt(nm, w) (see ../hopt.jl)
#   --resume    : Resume from a previous checkpoint if interrupted
#
# Prerequisites:
#   Run cubic_spectrum_dmrg_build_mpo.jl first to build the MPO cache.
#   Or use run_dmrg.sh which handles both steps automatically.
#
# Output:
#   dmrg_Z<Z>_P2<P2>_nm<nm>_w<w>_h<h>.dat  (binary)
#   dmrg_Z<Z>_P2<P2>_nm<nm>_w<w>_h<h>.txt  (human-readable)
#
# Checkpointing:
#   A checkpoint is saved after every sweep so the run can be resumed with
#   --resume if interrupted. Checkpoints are deleted on successful completion.
#
# Note: ITensors supports at most 4 quantum numbers. We use:
#   Ne, Lz2, Z (total parity), P2 (N2 parity)

using FuzzifiED
using ITensors, ITensorMPS
using LinearAlgebra
using Serialization
using Printf

include(joinpath(@__DIR__, "..", "hopt.jl"))

FuzzifiED.ElementType = Float64

# Threading Configuration
# ------------------------------------------------------------
# Set BLAS threads to use all available threads (from JULIA_NUM_THREADS)
LinearAlgebra.BLAS.set_num_threads(1)

# Enable BlockSparse threading (critical for QN conservation parallelism)
ITensors.enable_threaded_blocksparse(true)

# Disable Strided threading to avoid oversubscription (let BLAS handle inner parallelism)
if isdefined(ITensors, :Strided)
    ITensors.Strided.set_num_threads(1)
end

println("="^60)
println("Threading Configuration:")
println("  Threads.nthreads() = $(Threads.nthreads())")
println("  BLAS.get_num_threads() = $(LinearAlgebra.BLAS.get_num_threads())")
println("  ITensors.using_threaded_blocksparse() = $(ITensors.using_threaded_blocksparse())")
println("="^60)

# Helper function to compute mid-chain von Neumann entropy
# IMPORTANT: Works on a COPY to avoid modifying the MPS during DMRG
function midchain_entropy(psi::MPS)
    psi_copy = copy(psi)  # Make a copy to avoid corrupting DMRG state
    N = length(psi_copy)
    mid = div(N, 2)
    orthogonalize!(psi_copy, mid)
    _, S, _ = svd(psi_copy[mid], (linkind(psi_copy, mid-1), siteind(psi_copy, mid)))
    SvN = 0.0
    for n in 1:dim(S, 1)
        p = S[n, n]^2
        if p > 1e-14
            SvN -= p * log(p)
        end
    end
    return SvN
end

# 1. Parse Arguments
resume_mode = "--resume" in ARGS
numeric_args = filter(a -> !startswith(a, "-"), ARGS)
target_Z  = length(numeric_args) > 0 ? parse(Int,     numeric_args[1]) : 0
target_P2 = length(numeric_args) > 1 ? parse(Int,     numeric_args[2]) : 0
n_states  = length(numeric_args) > 2 ? parse(Int,     numeric_args[3]) : 1
nm        = length(numeric_args) > 3 ? parse(Int,     numeric_args[4]) : 12
w         = length(numeric_args) > 4 ? parse(Float64, numeric_args[5]) : 0.6
h         = length(numeric_args) > 5 ? parse(Float64, numeric_args[6]) : default_h(nm, w)
println("="^60)
println("Cubic DMRG - Sector (Z=$target_Z, P2=$target_P2)")
println("System size: nm=$nm | States to compute: $n_states")
if resume_mode
    println("Resume mode: ENABLED")
end
println("="^60)

# 2. Parameters (w and h come from the command line; see section 1)
nf = 4
no = nm * nf

println("Model parameters: h=$(round(h, digits=4)), w=$w")

# 3. Load cached MPOs and sites (built by cubic_spectrum_dmrg_build_mpo.jl)
h_str = replace(@sprintf("%.3f", h), "." => "p")
w_str = replace(string(w), "." => "p")
mpo_cache_file = "mpo_cache_nm$(nm)_w$(w_str)_h$(h_str).jls"

if !isfile(mpo_cache_file)
    error("""
    MPO cache not found: $mpo_cache_file

    Run the MPO build script first:
      julia --project=. cubic_spectrum_dmrg_build_mpo.jl $nm $w $h

    Or use run_dmrg.sh which handles both steps automatically.
    """)
end

println("Loading cached MPOs from $mpo_cache_file...")
t_load = time()
mpo_cache = deserialize(mpo_cache_file)

# Guard against a stale cache built with different parameters
if mpo_cache.nm != nm || !isapprox(mpo_cache.w, w; atol = 1e-9) || !isapprox(mpo_cache.h, h; atol = 1e-9)
    error("""
    MPO cache parameter mismatch in $mpo_cache_file
      cache:     nm=$(mpo_cache.nm), w=$(mpo_cache.w), h=$(mpo_cache.h)
      requested: nm=$nm, w=$w, h=$h
    Delete the cache file and rebuild it.
    """)
end

sites = mpo_cache.sites
hmt = mpo_cache.hmt
l2_mpo = mpo_cache.l2_mpo
c2_mpo = mpo_cache.c2_mpo
println("  MPOs loaded ($(round(time() - t_load, digits=1))s)")

# 6. Initialize MPS for specific sector
state = ["0" for i=1:no]
for m in 1:nm
    state[(m-1)*nf + 4] = "1" # Fill with scalars (Z=0, P2=0)
end

if target_Z == 1 && target_P2 == 1
    state[(1-1)*nf + 4] = "0"; state[(1-1)*nf + 2] = "1" # Swap scalar for flavor 2 (Z=1, P2=1)
elseif target_Z == 1 && target_P2 == 0
    state[(1-1)*nf + 4] = "0"; state[(1-1)*nf + 1] = "1" # Swap scalar for flavor 1 (Z=1, P2=0)
elseif target_Z == 0 && target_P2 == 1
    state[(1-1)*nf + 4] = "0"; state[(1-1)*nf + 2] = "1" # Site 1: Flv 2 (Z=1, P2=1)
    state[(2-1)*nf + 4] = "0"; state[(2-1)*nf + 1] = "1" # Site 2: Flv 1 (Z=1, P2=0) -> Total Z=0, P2=1
end
psi_init = MPS(sites, state)

# 7. DMRG Parameters (scaled by system size)
# Larger systems need higher bond dimensions and more sweeps for convergence
# Convergence tolerances for early stopping (per paper methodology)
function get_dmrg_params(nm)
    if nm <= 12
        max_bond = 1500
        nsweeps = 20
        energy_tol = 1e-6
    elseif nm <= 14
        max_bond = 2000
        nsweeps = 26
        energy_tol = 1e-6
    elseif nm <= 16
        max_bond = 2500
        nsweeps = 30
        energy_tol = 2e-6
    elseif nm <= 20
        max_bond = 3000
        nsweeps = 36
        energy_tol = 5e-6
    else  # nm >= 22
        max_bond = 3500
        nsweeps = 40
        energy_tol = 1e-5
    end

    # Build maxdim schedule: ramp up more gradually, then hold
    ramp = [20, 50, 100, 200, 400, round(Int, max_bond/2)]
    hold = fill(max_bond, nsweeps - length(ramp))
    maxdim = vcat(ramp, hold)

    # Noise schedule: constant 1E-5 for ergodicity in phase space (per paper)
    # "The amplitude of the noise was generally kept at 10^-5"
    noise = fill(1E-5, nsweeps)

    # Cutoff: slightly relax for very large systems to manage memory
    cutoff = nm <= 14 ? 1E-8 : 1E-7

    return nsweeps, maxdim, cutoff, noise, energy_tol
end

nsweeps, maxdim, cutoff, noise, energy_tol = get_dmrg_params(nm)

println("DMRG params: nsweeps=$nsweeps, max_bond=$(maximum(maxdim)), cutoff=$cutoff, energy_tol=$energy_tol")

# Checkpoint functions (using Julia Serialization)
# Note: w and h are part of the name so that runs at different couplings or
# different polarisation parameters cannot pick up each other's checkpoints.
function checkpoint_filename(Z, P2, nm, state_n)
    return "checkpoint_Z$(Z)_P2$(P2)_nm$(nm)_w$(w_str)_h$(h_str)_state$(state_n).jls"
end

function save_checkpoint(filename, psi, sweep, energy, prev_energies, prev_entropies, converged_count)
    # Write to temp file first, then rename (atomic operation)
    tmpfile = filename * ".tmp"
    serialize(tmpfile, (psi=psi, sweep=sweep, energy=energy, prev_energies=prev_energies, prev_entropies=prev_entropies, converged_count=converged_count))
    mv(tmpfile, filename; force=true)
end

function load_checkpoint(filename)
    data = deserialize(filename)
    # Handle both old format (without entropies) and new format
    prev_entropies = hasproperty(data, :prev_entropies) ? data.prev_entropies : Float64[]
    converged_count = hasproperty(data, :converged_count) ? data.converged_count : 0
    return data.psi, data.sweep, data.energy, data.prev_energies, prev_entropies, converged_count
end

# Save/load completed state results for multi-state resume
function completed_state_filename(Z, P2, nm, state_n)
    return "completed_Z$(Z)_P2$(P2)_nm$(nm)_w$(w_str)_h$(h_str)_state$(state_n).jls"
end

function save_completed_state(Z, P2, nm, state_n, psi, E, L2, C2)
    filename = completed_state_filename(Z, P2, nm, state_n)
    tmpfile = filename * ".tmp"
    serialize(tmpfile, (psi=psi, E=E, L2=L2, C2=C2))
    mv(tmpfile, filename; force=true)
end

function load_completed_state(Z, P2, nm, state_n)
    filename = completed_state_filename(Z, P2, nm, state_n)
    if isfile(filename)
        data = deserialize(filename)
        return data.psi, data.E, data.L2, data.C2, true
    end
    return nothing, nothing, nothing, nothing, false
end

function cleanup_all_checkpoints(Z, P2, nm, n_states)
    # Clean up sweep checkpoints
    for n in 1:n_states
        ckpt = checkpoint_filename(Z, P2, nm, n)
        if isfile(ckpt); rm(ckpt; force=true); end
    end
    # Clean up completed state files
    for n in 1:n_states
        cfile = completed_state_filename(Z, P2, nm, n)
        if isfile(cfile); rm(cfile; force=true); end
    end
end

# Custom observer for checkpointing with dual convergence (saves at end of each full sweep)
# Convergence criterion: both energy AND mid-chain entropy changes below threshold (per paper)
mutable struct CheckpointObserver <: AbstractObserver
    filename::String
    prev_energies::Vector{Float64}
    prev_entropies::Vector{Float64}  # Track mid-chain entropies
    last_saved_sweep::Int
    # Early stopping parameters
    min_sweeps::Int           # Minimum sweeps before early stop (to allow bond dim ramp-up)
    energy_tol::Float64       # Energy change threshold for convergence
    entropy_tol::Float64      # Entropy change threshold for convergence
    converged_count::Int      # Count of consecutive converged sweeps
    required_converged::Int   # How many consecutive converged sweeps needed
end

# Constructor with defaults
function CheckpointObserver(filename, prev_energies, last_sweep; min_sweeps=10, energy_tol=1e-7, entropy_tol=1e-7, required_converged=2)
    CheckpointObserver(filename, prev_energies, Float64[], last_sweep, min_sweeps, energy_tol, entropy_tol, 0, required_converged)
end

function ITensorMPS.measure!(o::CheckpointObserver; psi, sweep, half_sweep, energy, kwargs...)
    # Only save at end of full sweep (half_sweep == 2 means right-to-left pass done)
    if half_sweep == 2 && sweep > o.last_saved_sweep
        push!(o.prev_energies, energy)
        
        # Compute and track mid-chain entropy
        S_mid = midchain_entropy(psi)
        push!(o.prev_entropies, S_mid)
        
        # Check for BOTH energy AND entropy convergence (per paper methodology)
        if length(o.prev_energies) >= 2 && length(o.prev_entropies) >= 2
            dE = abs(o.prev_energies[end] - o.prev_energies[end-1])
            dS = abs(o.prev_entropies[end] - o.prev_entropies[end-1])
            
            # Both must be converged
            if dE < o.energy_tol && dS < o.entropy_tol
                o.converged_count += 1
            else
                o.converged_count = 0
            end
        end
        
        save_checkpoint(o.filename, psi, sweep, energy, o.prev_energies, o.prev_entropies, o.converged_count)
        o.last_saved_sweep = sweep
        #println("  [Checkpoint saved: sweep $sweep, E=$(round(energy, digits=6)), S=$(round(S_mid, digits=4))]")
    end
end

function ITensorMPS.checkdone!(o::CheckpointObserver; sweep, kwargs...)
    # Early stopping: if we've done enough sweeps and BOTH energy AND entropy have converged
    if sweep >= o.min_sweeps && o.converged_count >= o.required_converged
        println("    [Early stop: energy AND entropy converged after $sweep sweeps]")
        return true
    end
    return false
end

# Subspace re-diagonalization: form small H, S, L², C² matrices in the
# DMRG basis and solve the generalized eigenvalue problem to remove the
# systematic energy bias introduced by the penalty projector.
function subspace_rediagonalize(states::Vector{MPS}, hmt, l2_mpo, c2_mpo)
    n = length(states)
    if n <= 1
        return nothing
    end

    println("\n" * "="^60)
    println("Subspace re-diagonalization ($n states)")
    println("="^60)

    # Build full n×n matrices (off-diagonals needed for proper transformation)
    H_sub  = zeros(Float64, n, n)
    S_sub  = zeros(Float64, n, n)
    L2_sub = zeros(Float64, n, n)
    C2_sub = zeros(Float64, n, n)

    for i in 1:n
        for j in i:n
            H_sub[i, j]  = real(inner(states[i]', hmt,    states[j]))
            S_sub[i, j]  = real(inner(states[i],          states[j]))
            L2_sub[i, j] = real(inner(states[i]', l2_mpo, states[j]))
            C2_sub[i, j] = real(inner(states[i]', c2_mpo, states[j]))

            if i != j  # Hermitian
                H_sub[j, i]  = H_sub[i, j]
                S_sub[j, i]  = S_sub[i, j]
                L2_sub[j, i] = L2_sub[i, j]
                C2_sub[j, i] = C2_sub[i, j]
            end
        end
    end

    println("  S_sub diagonal (overlaps): ", [round(S_sub[i,i], digits=8) for i in 1:n])
    println("  S_sub off-diag (max |S_ij|, i≠j): ",
            round(maximum(abs(S_sub[i,j]) for i in 1:n for j in 1:n if i != j), digits=8))
    κ = cond(S_sub)
    println("  cond(S_sub) = ", round(κ, sigdigits=4))
    if κ > 1e10
        println("  WARNING: S_sub is ill-conditioned; re-diag results may be unreliable.")
    end

    # Solve generalised eigenvalue problem H c = E S c
    vals, vecs = eigen(Hermitian(H_sub), Hermitian(S_sub))

    # Transform observables: O_k = c_k' * O_sub * c_k
    results_rediag = []
    for k in 1:n
        c = vecs[:, k]
        E_k  = vals[k]
        L2_k = real(c' * L2_sub * c)
        C2_k = real(c' * C2_sub * c)
        push!(results_rediag, (state=k, E=E_k, L2=L2_k, C2=C2_k))
        println("  State $k: E = $(round(E_k, digits=8))  L² = $(round(L2_k, digits=4))  C² = $(round(C2_k, digits=4))")
    end

    return (
        results_rediag = results_rediag,
        H_sub  = H_sub,
        S_sub  = S_sub,
        L2_sub = L2_sub,
        C2_sub = C2_sub,
        eigvals = vals,
        eigvecs = vecs
    )
end

# Storage for results
results = []
prev_states = MPS[]

# Initial link dimension for randomMPS (scales with system size)
init_linkdim = nm <= 12 ? 10 : (nm <= 18 ? 20 : 40)

for n in 1:n_states
    println("\nComputing state $n / $n_states...")

    # First check if this state was already completed in a previous run
    loaded_psi, loaded_E, loaded_L2, loaded_C2, already_done = load_completed_state(target_Z, target_P2, nm, n)
    if already_done
        println("  State $n already completed (loaded from checkpoint)")
        println("  E = $loaded_E | L^2 = $loaded_L2 | C^2 = $loaded_C2")
        push!(results, (state=n, E=loaded_E, L2=loaded_L2, C2=loaded_C2))
        push!(prev_states, loaded_psi)
        continue
    end

    ckpt_file = checkpoint_filename(target_Z, target_P2, nm, n)
    start_sweep = 1
    psi_start = nothing
    prev_energies = Float64[]
    prev_entropies = Float64[]
    converged_count = 0

    # Check for checkpoint to resume from
    if resume_mode && isfile(ckpt_file)
        println("  Found checkpoint: $ckpt_file")
        psi_start, last_sweep, last_energy, prev_energies, prev_entropies, converged_count = load_checkpoint(ckpt_file)
        start_sweep = last_sweep + 1
        println("  Resuming from sweep $start_sweep (last E=$last_energy, converged_count=$converged_count)")

        if start_sweep > nsweeps
            println("  Already completed all sweeps, skipping to results...")
            E = last_energy
            psi = psi_start
            @goto compute_observables
        end
    else
        # Fresh start
        if n == 1
            psi_start = psi_init
        else
            psi_start = randomMPS(sites, state; linkdims=init_linkdim)
        end
    end

    # Adjust DMRG parameters for remaining sweeps
    remaining_sweeps = nsweeps - start_sweep + 1
    maxdim_remaining = maxdim[start_sweep:end]
    noise_remaining = noise[start_sweep:end]

    # Create checkpoint observer with restored state
    observer = CheckpointObserver(ckpt_file, prev_energies, start_sweep - 1; energy_tol=energy_tol)
    observer.prev_entropies = prev_entropies
    observer.converged_count = converged_count

    # Run DMRG
    if isempty(prev_states)
        # Ground state: use two-argument form
        E, psi = dmrg(hmt, psi_start;
                      nsweeps=remaining_sweeps,
                      maxdim=maxdim_remaining,
                      cutoff=[cutoff],
                      noise=noise_remaining,
                      observer=observer)
    else
        # Excited states: use three-argument form with penalty
        E, psi = dmrg(hmt, prev_states, psi_start;
                      nsweeps=remaining_sweeps,
                      maxdim=maxdim_remaining,
                      cutoff=[cutoff],
                      noise=noise_remaining,
                      weight=100.0,
                      observer=observer)

    end

    @label compute_observables
    L2 = real(inner(psi', l2_mpo, psi))
    C2 = real(inner(psi', c2_mpo, psi))

    println("State $n: E = $E | L^2 = $L2 | C^2 = $C2")
    push!(results, (state=n, E=E, L2=L2, C2=C2))
    push!(prev_states, psi)

    # Save completed state for multi-state resume (keep until all states done)
    save_completed_state(target_Z, target_P2, nm, n, psi, E, L2, C2)
    println("  [Completed state $n saved for resume]")
    
    # Remove sweep checkpoint (no longer needed, we have the completed state)
    if isfile(ckpt_file)
        rm(ckpt_file; force=true)
    end
end

# Subspace re-diagonalization (only meaningful for n_states > 1)
rediag_data = subspace_rediagonalize(prev_states, hmt, l2_mpo, c2_mpo)

# All states completed - clean up all checkpoint files
cleanup_all_checkpoints(target_Z, target_P2, nm, n_states)
println("[All checkpoints cleaned up]")

# 8. Save Results (h_str, w_str defined earlier for MPO cache)
outfile = "dmrg_Z$(target_Z)_P2$(target_P2)_nm$(nm)_w$(w_str)_h$(h_str).dat"
println("\n" * "="^60)
println("Saving results to $outfile...")

output_data = (
    Z = target_Z,
    P2 = target_P2,
    nm = nm,
    h = h,
    w = w,
    n_states = n_states,
    results = results
)

if !isnothing(rediag_data)
    output_data = merge(output_data, (
        results_rediag = rediag_data.results_rediag,
        H_sub  = rediag_data.H_sub,
        S_sub  = rediag_data.S_sub,
        L2_sub = rediag_data.L2_sub,
        C2_sub = rediag_data.C2_sub
    ))
    println("  [Re-diagonalized results included]")
end

serialize(outfile, output_data)

# Also write human-readable summary
txtfile = "dmrg_Z$(target_Z)_P2$(target_P2)_nm$(nm)_w$(w_str)_h$(h_str).txt"
open(txtfile, "w") do io
    println(io, "Cubic DMRG Results")
    println(io, "="^40)
    println(io, "Sector: Z=$target_Z, P2=$target_P2")
    println(io, "System size: nm=$nm")
    println(io, "Parameters: h=$(round(h, digits=4)), w=$w")
    println(io, "="^40)

    println(io, "\nDMRG Results:")
    println(io, "State | Energy | L^2 | C^2")
    for r in results
        println(io, "$(r.state) | $(r.E) | $(r.L2) | $(r.C2)")
    end

    if !isnothing(rediag_data)
        println(io, "\nSubspace Re-diagonalized Results:")
        println(io, "State | Energy | L^2 | C^2")
        for r in rediag_data.results_rediag
            println(io, "$(r.state) | $(r.E) | $(r.L2) | $(r.C2)")
        end

        println(io, "\nOverlap matrix S_sub:")
        for i in 1:n_states
            println(io, "  ", [round(rediag_data.S_sub[i,j], digits=8) for j in 1:n_states])
        end
    end
end
println("Saved human-readable summary to $txtfile")
println("="^60)
