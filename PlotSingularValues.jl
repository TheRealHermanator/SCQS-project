using DelimitedFiles
using Plots
using TNRKit
using TensorKit
include("6VertexModel.jl")
include("ModifiedTRGflow.jl")


# Variables
nmaxiter = 25



function compute_TRG_data_sing(ac_val, bc_val; ntruncdim=16, nmaxiter=nmaxiter)
    T = six_vertex_tensor(ac_val, bc_val, 1)  # sixvertex(a=ac_val, b=bc_val, c=1) # six_vertex_tensor(ac_val, bc_val, 1)
    scheme = BTRG(T)
    data, sing = run_sing!(scheme, truncdim(ntruncdim), maxiter(nmaxiter); verbosity=1)
    return data, sing
end

data, sing = compute_TRG_data_sing(1.5, 1.5; ntruncdim=16, nmaxiter=25)


# Prepare data for scatter plot
x = Int[]
y = Float64[]

for i in 1:nmaxiter-1
    svals = values(sing[i])[1]  # Extract the singular values for the i-th iteration
    for val in svals
        push!(x, i)              # Store the iteration index
        push!(y, val)            # Store the singular value
    end
end

# Scatter plot
scatter(x, y,
    xlabel = "Iteration",
    ylabel = "Singular value",
    title = "Singular values over iterations (log scale)",
    yscale = :log10,
    legend = false)
