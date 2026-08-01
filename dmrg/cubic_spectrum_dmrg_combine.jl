# Combine DMRG results from all 4 sector jobs
#
# Usage:
#   julia cubic_spectrum_dmrg_combine.jl [nm] [w] [h]
#
# Arguments:
#   nm: system size (default: 12)
#   w:  cubic coupling (default: 0.6)
#   h:  polarization parameter (default: h_opt(nm, w), see ../hopt.jl)
#
# Expects input files: dmrg_Z<Z>_P2<P2>_nm<nm>_w<w>_h<h>.dat for Z,P2 in {0,1}
#
# Output:
#   dmrg_combined_nm<nm>_w<w>_h<h>.txt

using Serialization
using Printf

include(joinpath(@__DIR__, "..", "hopt.jl"))

nm = length(ARGS) >= 1 ? parse(Int,     ARGS[1]) : 12
w  = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.6
h  = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : default_h(nm, w)

# Format h and w for filename matching
h_str = replace(@sprintf("%.3f", h), "." => "p")
w_str = replace(string(w), "." => "p")

println("="^60)
println("Combining DMRG results for nm=$nm, h=$(round(h, digits=2)), w=$w")
println("="^60)

# Collect all results
all_results = []
missing_sectors = []

for Z in [0, 1], P2 in [0, 1]
    filename = "dmrg_Z$(Z)_P2$(P2)_nm$(nm)_w$(w_str)_h$(h_str).dat"
    if !isfile(filename)
        push!(missing_sectors, (Z=Z, P2=P2))
        println("WARNING: Missing $filename")
        continue
    end

    data = deserialize(filename)

    # Prefer re-diagonalized results when available
    has_rediag = hasproperty(data, :results_rediag)
    results_to_use = has_rediag ? data.results_rediag : data.results

    for r in results_to_use
        push!(all_results, (Z=Z, P2=P2, state=r.state, E=r.E, L2=r.L2, C2=r.C2))
    end

    source_label = has_rediag ? "re-diag" : "original"
    println("Loaded (Z=$Z, P2=$P2): $(length(results_to_use)) states [$source_label]")
end

if !isempty(missing_sectors)
    println("\nWARNING: Missing $(length(missing_sectors)) sector file(s)")
    println("Missing: $missing_sectors")
end

if isempty(all_results)
    println("ERROR: No results found!")
    exit(1)
end

# Sort by energy
sort!(all_results, by=x -> x.E)

# Find vacuum energy for scaling
E_vac = all_results[1].E

# Find phi energy (Z=1, C2~2) for scaling
phi_states = filter(x -> x.Z == 1 && abs(x.C2 - 2) < 0.5, all_results)
E_phi = isempty(phi_states) ? all_results[2].E : phi_states[1].E

# Write combined output
outfile = "dmrg_combined_nm$(nm)_w$(w_str)_h$(h_str).txt"
println("\nWriting combined results to $outfile...")

open(outfile, "w") do io
    println(io, "Cubic DMRG Combined Results")
    println(io, "="^60)
    println(io, "System size: nm=$nm")
    println(io, "Parameters: h=$(round(h, digits=4)), w=$w")
    println(io, "="^60)
    println(io)
    println(io, "All States (sorted by energy):")
    println(io, "ScaledE | Energy | L^2 | C^2 | Z | P2 | State#")
    println(io, "-"^70)

    for r in all_results
        scaled_E = 0.5189 * (r.E - E_vac) / (E_phi - E_vac)
        println(io, "$(round(scaled_E, digits=4)) | $(round(r.E, digits=6)) | $(round(r.L2, digits=4)) | $(round(r.C2, digits=4)) | $(r.Z) | $(r.P2) | $(r.state)")
    end

    println(io)
    println(io, "Reference energies:")
    println(io, "  E_vacuum = $E_vac")
    println(io, "  E_phi = $E_phi")
end

# Print summary to console
println("\n" * "="^60)
println("RESULTS SUMMARY (sorted by energy)")
println("="^60)
println("ScaledE | Energy | L^2 | C^2 | Z | P2")
println("-"^60)
for r in all_results[1:min(10, length(all_results))]
    scaled_E = 0.5189 * (r.E - E_vac) / (E_phi - E_vac)
    println("$(round(scaled_E, digits=4)) | $(round(r.E, digits=4)) | $(round(r.L2, digits=2)) | $(round(r.C2, digits=2)) | $(r.Z) | $(r.P2)")
end
if length(all_results) > 10
    println("... ($(length(all_results) - 10) more states)")
end
println("="^60)
println("\nOutput written to: $outfile")
