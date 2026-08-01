# Optimal polarisation parameter h_opt(nm, w)
#
# Shared by the ED and DMRG drivers so that all of them agree on the value of h
# for a given (nm, w). Include with
#
#   include(joinpath(@__DIR__, "..", "hopt.jl"))
#
# The table holds the values actually used to produce the published spectra:
# nm = 12..18 come from the DMRG optimisation of g_S(h) = 0, and nm = 20, 22 are
# extrapolated with h_opt(nm) = h_c + A / nm^(b/2). The fit parameters below are
# a fallback for system sizes that are not tabulated; note that they reproduce
# the tabulated values only to about +/- 0.002, so the table takes precedence.

const HOPT_TABLE = Dict(
    0.2 => Dict(12 => 15.394, 14 => 15.314, 16 => 15.260, 18 => 15.218, 20 => 15.187, 22 => 15.162),
    0.4 => Dict(12 => 15.358, 14 => 15.276, 16 => 15.220, 18 => 15.177, 20 => 15.146, 22 => 15.121),
    0.6 => Dict(12 => 15.323, 14 => 15.239, 16 => 15.181, 18 => 15.140, 20 => 15.110, 22 => 15.087),
    0.8 => Dict(12 => 15.288, 14 => 15.202, 16 => 15.142, 18 => 15.097, 20 => 15.063, 22 => 15.036),
)

# (h_c, A, b) of the finite-size fit h_opt(nm) = h_c + A / nm^(b/2)
const HOPT_FIT = Dict(
    0.2 => (14.985, 12.69, 2.77),
    0.4 => (14.957, 15.40, 2.94),
    0.6 => (14.955, 24.69, 3.39),
    0.8 => (14.855, 15.40, 2.87),
)

"""
    default_h(nm, w)

Optimal polarisation parameter for system size `nm` at cubic coupling `w`.

Returns the tabulated value when one exists, otherwise evaluates the
finite-size fit. Errors if `w` is not one of the couplings studied, in which
case h must be supplied explicitly on the command line.
"""
function default_h(nm::Integer, w::Real)
    wk = round(float(w), digits = 6)
    if haskey(HOPT_TABLE, wk) && haskey(HOPT_TABLE[wk], nm)
        return HOPT_TABLE[wk][nm]
    elseif haskey(HOPT_FIT, wk)
        h_c, A, b = HOPT_FIT[wk]
        return h_c + A * nm^(-b / 2)
    else
        error("""
        No default h is available for w = $w.
        Couplings with tabulated fits: $(sort(collect(keys(HOPT_FIT)))).
        Pass h explicitly as a command-line argument.
        """)
    end
end
