# Exact Diagonalization of the Cubic Model on the Fuzzy Sphere
#
# Loops over all 16 symmetry sectors (Z × P2 × R × X), diagonalises each,
# and writes the combined spectrum to disk.
#
# Usage:
#   julia cubic_spectrum_ed.jl [nm] [w] [h]
#
# Arguments:
#   nm  : system size (default: 8)
#   w   : cubic anisotropy coupling (default: 0.6)
#   h   : polarization parameter (default: h_opt(nm, w), see ../hopt.jl)
#
# Output:
#   spectrum_nm<nm>_w<w>_h<h>.txt   (scaled energies, quantum numbers, irrep labels)
#
# Example:
#   julia cubic_spectrum_ed.jl 12 0.2

using FuzzifiED
using LinearAlgebra
using Serialization
using Printf

include(joinpath(@__DIR__, "..", "hopt.jl"))

FuzzifiED.ElementType = Float64

# ── Parse arguments ────────────────────────────────────────────────────────────
nm = length(ARGS) >= 1 ? parse(Int,     ARGS[1]) : 8
w  = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.6

# Default h is the optimal value for this (nm, w); see ../hopt.jl
h = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : default_h(nm, w)

println("="^60)
println("Cubic ED Spectrum")
println("  nm = $nm  |  w = $w  |  h = $(round(h, digits=4))")
println("="^60)

# ── Matrices ───────────────────────────────────────────────────────────────────
mat_0 = diagm([0, 0, 0, 1.0])
mat_A = [
    [ 0  1 0 0 ; 1 0  1 0 ; 0 1  0 0 ; 0 0 0 0 ] / √2,
    [ 0 -1 0 0 ; 1 0 -1 0 ; 0 1  0 0 ; 0 0 0 0 ] / √2 * im,
    [ 1  0 0 0 ; 0 0  0 0 ; 0 0 -1 0 ; 0 0 0 0 ]
]
mat_V = [
    [ 0 0 0  1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0 -1 0 ] / √2,
    [ 0 0 0 -1 ; 0 0 0 0 ; 0 0 0 -1 ; 1 0  1 0 ] / √2 * im,
    [ 0 0 0  0 ; 0 0 0 1 ; 0 0 0  0 ; 0 1  0 0 ]
]
mat_Px = [ 0.5 0 -0.5 0 ; 0 0 0 0 ; -0.5 0 0.5 0 ; 0 0 0 0 ]
mat_Py = [ 0.5 0  0.5 0 ; 0 0 0 0 ;  0.5 0 0.5 0 ; 0 0 0 0 ]
mat_Pz = diagm([0.0, 1.0, 0.0, 0.0])
mats_cubic = [mat_Px, mat_Py, mat_Pz]

nf = 4
no = nm * nf

# ── Quantum number setup ───────────────────────────────────────────────────────
qnd = [
    GetNeQNDiag(no),
    GetLz2QNDiag(nm, nf),
    GetFlavQNDiag(nm, nf, Dict([f => 1 for f = 1 : nf - 1]), 0, 2),
    GetFlavQNDiag(nm, nf, Dict([2 => 1]), 0, 2)
]
qnf = [
    GetRotyQNOffd(nm, nf),
    GetFlavPermQNOffd(nm, nf, Dict([1 => 3, 3 => 1]))
]

# Hamiltonian terms (sector-independent)
tms_o3    = GetDenIntTerms(nm, nf, [6.5, 1.0]) - 1.4 * 2 * GetDenIntTerms(nm, nf, [6.5, 1.0], mat_V) - 2 * h * GetPolTerms(nm, nf, mat_0)
tms_cubic = w * GetDenIntTerms(nm, nf, [6.5, 1.0], mats_cubic)
tms_hmt   = SimplifyTerms(tms_o3 + tms_cubic)

tms_l2 = GetL2Terms(nm, nf)
tms_c2 = 4 * GetC2Terms(nm, nf, mat_A)

# ── Struct for a single energy level ──────────────────────────────────────────
struct EnergyLevel
    E::Float64; L2::Float64; C2::Float64
    Z::Int; P2::Int; R::Int; X::Int
end

# ── Loop over all 16 sectors ───────────────────────────────────────────────────
all_levels = EnergyLevel[]

sectors = [(Z=Z, P2=P2, R=R, X=X)
           for Z in [0,1] for P2 in [0,1] for R in [1,-1] for X in [1,-1]]

for (idx, sec) in enumerate(sectors)
    Z, P2, R, X = sec.Z, sec.P2, sec.R, sec.X
    print("Sector $idx/16  (Z=$Z P2=$P2 R=$R X=$X): ")

    cfs = Confs(no, [nm, 0, Z, P2], qnd)
    bs  = Basis(cfs, [R, X], qnf)
    print("dim=$(bs.dim)  ")

    hmt  = Operator(bs, tms_hmt)
    enrg, st = GetEigensystem(OpMat(hmt), 8)

    l2  = Operator(bs, tms_l2)
    c2  = Operator(bs, tms_c2)
    l2m = OpMat(l2)
    c2m = OpMat(c2)

    for i in eachindex(enrg)
        push!(all_levels, EnergyLevel(
            real(enrg[i]),
            real(st[:,i]' * l2m * st[:,i]),
            real(st[:,i]' * c2m * st[:,i]),
            Z, P2, R, X
        ))
    end
    println("done")
end

sort!(all_levels, by = x -> x.E)
println()
println("Total states: $(length(all_levels))")

# ── Scale energies ─────────────────────────────────────────────────────────────
E_vac = all_levels[1].E
phis  = filter(x -> x.Z == 1 && abs(x.C2 - 2) < 0.5, all_levels)
E_phi = isempty(phis) ? all_levels[2].E : phis[1].E

# ── Cluster into degenerate multiplets ────────────────────────────────────────
function build_clusters(levels; tol=1e-5)
    clusters = Vector{EnergyLevel}[]
    cur = [levels[1]]
    for i in 2:length(levels)
        if abs(levels[i].E - levels[i-1].E) < tol
            push!(cur, levels[i])
        else
            push!(clusters, cur)
            cur = [levels[i]]
        end
    end
    push!(clusters, cur)
    return clusters
end

clusters = build_clusters(all_levels)

function label_cluster(c)
    rep  = c[1]
    deg  = length(c)
    par  = rep.Z == 0 ? "g" : "u"
    lbl  = deg == 1 ? "A1$par/A2$par" :
           deg == 2 ? "E$par"         :
           deg == 3 ? "T1$par/T2$par" :
           deg == 4 ? "deg(4)"        :
           "deg($deg)"
    if deg == 1 && abs(rep.C2) < 0.5;       lbl = "A1$par"; end
    if deg == 3 && abs(rep.C2 - 2) < 0.5;   lbl = "T1$par"; end
    scaled = round(0.5189 * (rep.E - E_vac) / (E_phi - E_vac), digits=6)
    return scaled, rep.L2, rep.C2, deg, lbl
end

unique_states = [label_cluster(c) for c in clusters]

# ── Write output ───────────────────────────────────────────────────────────────
w_str  = replace(string(w), "." => "p")
h_str  = replace(@sprintf("%.3f", h), "." => "p")
outfile = "spectrum_nm$(nm)_w$(w_str)_h$(h_str).txt"

open(outfile, "w") do io
    println(io, "# Cubic ED Spectrum  nm=$nm  w=$w  h=$(round(h, digits=6))")
    println(io, "# ScaledE normalised: E_vac=0, E_phi=0.5189")
    println(io, "#")
    println(io, "# ScaledE  | L2      | C2      | Degen | Irrep")
    println(io, "#" * "-"^55)
    for (se, l2, c2, deg, lbl) in unique_states
        @printf(io, "  %8.5f  | %6.3f  | %6.3f  | %5d | %s\n", se, l2, c2, deg, lbl)
    end
end

# ── Print summary ──────────────────────────────────────────────────────────────
println("="^60)
println("Spectrum summary (unique irreps, sorted by energy):")
println("  ScaledE   |   L2   |   C2   | Degen | Irrep")
println("  " * "-"^52)
for (se, l2, c2, deg, lbl) in unique_states[1:min(12, end)]
    @printf("  %8.5f  | %6.3f | %6.3f | %5d | %s\n", se, l2, c2, deg, lbl)
end
println()
println("Output written to: $outfile")
println("="^60)
