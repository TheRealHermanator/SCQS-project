using DelimitedFiles
using Plots
using TNRKit
using TensorKit
include("6VertexModel.jl")

"""
    compute_TRG_data(ac_val, bc_val, ntruncdim=16, nmaxiter=25)

Calculate the Tensor Renormalization Group (TRG) data based on the values of (a/c) and (b/c).

# Arguments
- `ac_val::Float64`: The ratio a/c.
- `bc_val::Float64`: The ratio b/c.
- `ntruncdim::Int`: The truncation dimension for the TRG algorithm (default: 16).
- `nmaxiter::Int`: The maximum number of iterations for the TRG algorithm (default: 25).

# Returns
- `data::Vector`: The TRG data computed for the given parameters.
"""
function compute_TRG_data(ac_val, bc_val, ntruncdim=16, nmaxiter=25)
    T = six_vertex_tensor(ac_val, bc_val, 1)  # sixvertex(a=ac_val, b=bc_val, c=1) # six_vertex_tensor(ac_val, bc_val, 1)
    scheme = BTRG(T)
    data = run!(scheme, truncdim(ntruncdim), maxiter(nmaxiter); verbosity=1)
    return data
end

"""
    compute_lnz(ac_val, bc_val, ntruncdim=16, nmaxiter=25)

Calculate the natural logarithm of the partition function Z (ln(Z)) for given values of (a/c) and (b/c).

# Arguments
- `ac_val::Float64`: The ratio a/c.
- `bc_val::Float64`: The ratio b/c.
- `ntruncdim::Int`: The truncation dimension for the TRG algorithm (default: 16).
- `nmaxiter::Int`: The maximum number of iterations for the TRG algorithm (default: 25).

# Returns
- `lnz::Float64`: The natural logarithm of the partition function Z.
"""
function compute_lnz(ac_val, bc_val, ntruncdim=16, nmaxiter=25)
    # Get data
    data = compute_TRG_data(ac_val, bc_val, ntruncdim, nmaxiter)

    # Calculate ln(Z)
    lnz = 0.0
    for (i, d) in enumerate(data)
        lnz += log(d) * 2.0^(1 - i)
    end
    return lnz
end

"""
    get_lnZ_grid(ac_vals, bc_vals, ntruncdim=16, nmaxiter=25)

Generate a grid of ln(Z) values for different combinations of (a/c) and (b/c).

# Arguments
- `ac_vals::AbstractVector`: A vector of a/c values.
- `bc_vals::AbstractVector`: A vector of b/c values.
- `ntruncdim::Int`: The truncation dimension for the TRG algorithm (default: 16).
- `nmaxiter::Int`: The maximum number of iterations for the TRG algorithm (default: 25).

# Returns
- `grid::Matrix`: A 2D matrix where each element corresponds to ln(Z) for a specific (a/c, b/c) pair.
"""
function get_lnZ_grid(ac_vals, bc_vals, ntruncdim=16, nmaxiter=25)
    # Define grid
    grid = zeros(length(ac_vals), length(bc_vals))

    # Fill in the grid
    for (i, ac) in enumerate(ac_vals)
        for (j, bc) in enumerate(bc_vals)
            println("Major iteration: $i\nMinor iteration: $j")
            lnZ = compute_lnz(ac, bc, ntruncdim, nmaxiter)
            grid[i, j] = lnZ
        end
    end
    return grid
end

# Constants for the computation
ntruncdim = 16  # Truncation dimension for the TRG algorithm
nmaxiter = 20   # Maximum number of iterations for the TRG algorithm
len = 200         # Number of points in the range for a/c and b/c

# Define the range of a/c and b/c values
ac_vals = range(0.01, 1.99; length=len)
bc_vals = range(0.01, 1.99; length=len)

# Compute the grid of ln(Z) values
lnz_grid = get_lnZ_grid(ac_vals, bc_vals, ntruncdim, nmaxiter)

# Store the computed data to CSV files
writedlm("6VertexModel_Jarid/Data/ac_vals_$len.csv", ac_vals)
writedlm("6VertexModel_Jarid/Data/bc_vals_$len.csv", bc_vals)
writedlm("6VertexModel_Jarid/Data/lnZ_$len.csv", lnz_grid)

# Generate a heatmap of ln(Z) values
heatmap(bc_vals, ac_vals, lnz_grid;
        xlabel="b/c", ylabel="a/c", title="ln(Z) Energy per Site",
        colorbar_title="ln(Z)", aspect_ratio=:equal)