# Build and cache MPOs for Cubic DMRG
# Usage: julia cubic_spectrum_dmrg_build_mpo.jl [nm] [w] [h]
#
# Arguments:
#   nm : system size (default: 12)
#   w  : cubic anisotropy coupling (default: 0.6)
#   h  : polarization parameter (default: h_opt(nm, w), see ../hopt.jl)
#
# This script builds the Hamiltonian and observable MPOs and caches them
# to disk. Run this BEFORE launching the DMRG sector jobs. If the cache
# already exists it is left untouched.
#
# Output:
#   mpo_cache_nm<nm>_w<w>_h<h>.jls

using FuzzifiED
using ITensors, ITensorMPS
using LinearAlgebra
using Serialization
using Printf

include(joinpath(@__DIR__, "..", "hopt.jl"))

FuzzifiED.ElementType = Float64

# Parse arguments
nm = length(ARGS) >= 1 ? parse(Int,     ARGS[1]) : 12
w  = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.6
h  = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : default_h(nm, w)

println("="^60)
println("Building MPO Cache for nm=$nm")
println("="^60)

# Parameters
nf = 4
no = nm * nf

println("Model parameters: h=$(round(h, digits=4)), w=$w")

# Cache filename; must match the one the DMRG driver looks for
h_str = replace(@sprintf("%.3f", h), "." => "p")
w_str = replace(string(w), "." => "p")
cache_file = "mpo_cache_nm$(nm)_w$(w_str)_h$(h_str).jls"

# Skip the (expensive) build if the cache is already present
if isfile(cache_file)
    println("\nMPO cache already exists: $cache_file")
    println("Nothing to do.")
    exit(0)
end

# Matrices for Hamiltonian
mat_0 = diagm([0, 0, 0, 1.0])
mat_V = [
    [ 0 0 0  1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0 -1 0 ] / sqrt(2),
    [ 0 0 0 -1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0  1 0 ] / sqrt(2) * im,
    [ 0 0 0  0 ; 0 0 0 1 ; 0 0 0  0 ; 0 1  0 0 ]
]
mat_Px = [ [0.5 0 -0.5 0]; [0 0 0 0]; [-0.5 0 0.5 0]; [0 0 0 0] ]
mat_Py = [ [0.5 0 0.5 0]; [0 0 0 0]; [0.5 0 0.5 0]; [0 0 0 0] ]
mat_Pz = diagm([0.0, 1.0, 0.0, 0.0])
mats_cubic = [mat_Px, mat_Py, mat_Pz]

mat_A = [
    [ 0  1 0 0 ; 1 0  1 0 ; 0 1  0 0 ; 0 0 0 0 ] / sqrt(2),
    [ 0 -1 0 0 ; 1 0 -1 0 ; 0 1  0 0 ; 0 0 0 0 ] / sqrt(2) * im,
    [ 1  0 0 0 ; 0 0  0 0 ; 0 0 -1 0 ; 0 0 0 0 ]
]

# Site definition (same QN structure used by all sectors)
println("\nBuilding site indices...")
sites = GetSites([
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, nf),
    GetFlavQNDiag(nm, nf, Dict([f => 1 for f = 1:3]), 0, 2), # Total Parity Z
    GetFlavQNDiag(nm, nf, Dict([2 => 1]), 1, 2)             # N2 Parity P2
])

# Build MPOs
println("\nBuilding Hamiltonian MPO...")
t_start = time()
terms_o3 = (GetDenIntTerms(nm, nf, [6.5, 1.0])
            - 1.4 * 2 * GetDenIntTerms(nm, nf, [6.5, 1.0], mat_V)
            - 2 * h * GetPolTerms(nm, nf, mat_0))
terms_cubic = w * GetDenIntTerms(nm, nf, [6.5, 1.0], mats_cubic)
hmt = MPO(OpSum(SimplifyTerms(terms_o3 + terms_cubic)), sites)
println("  Done ($(round(time() - t_start, digits=1))s)")

println("Building L² MPO...")
t_l2 = time()
l2_mpo = MPO(OpSum(GetL2Terms(nm, nf)), sites)
println("  Done ($(round(time() - t_l2, digits=1))s)")

println("Building C² MPO...")
t_c2 = time()
c2_mpo = MPO(OpSum(4 * GetC2Terms(nm, nf, mat_A)), sites)
println("  Done ($(round(time() - t_c2, digits=1))s)")

total_time = time() - t_start
println("\nTotal MPO build time: $(round(total_time, digits=1))s")

# Save cache (cache_file was computed above)
println("\nSaving MPO cache to $cache_file...")
serialize(cache_file, (
    nm = nm,
    nf = nf,
    h = h,
    w = w,
    sites = sites,
    hmt = hmt,
    l2_mpo = l2_mpo,
    c2_mpo = c2_mpo
))

# Report file size
file_size = filesize(cache_file) / 1024^2
println("  Cache file size: $(round(file_size, digits=1)) MB")

println("\n" * "="^60)
println("MPO cache ready: $cache_file")
println("="^60)
