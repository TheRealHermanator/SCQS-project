using DelimitedFiles
using Plots
using TNRKit
using TensorKit
include("6VertexModel.jl")
include("ModifiedTRGflow.jl")


# Variables
nmaxiter = 30
a = 1.5
b = 1.5

function extract_svals(sing)
    x = Int[]
    y = Float64[]
    for i in 1:nmaxiter-1
        svals = sing[i]  # Extract the singular values for the i-th iteration
        maxval = maximum(svals)
        for val in svals
            push!(x, i)              # Store the iteration index
            push!(y, val / maxval)            # Store the singular value
        end
    end
    return x, y
end

#The following code is repeated for the three different schemes, I was too lazy to find out how to vectorize it
function compute_TRG_data_sing(ac_val, bc_val; ntruncdim=16, nmaxiter=nmaxiter)
    T = six_vertex_tensor(ac_val, bc_val, 1)  # sixvertex(a=ac_val, b=bc_val, c=1) # six_vertex_tensor(ac_val, bc_val, 1)
    scheme = TRG(T)
    data, sing = run_sing!(scheme, truncdim(ntruncdim), maxiter(nmaxiter); verbosity=0)
    return data, sing
end

data1, sing1 = compute_TRG_data_sing(a, b; ntruncdim=16, nmaxiter=nmaxiter)
x1, y1 = extract_svals(sing1)
# Scatter plot
plt = scatter(x1, y1,
    xlabel = "Iteration",
    ylabel = "Singular value",
    title = "Singular Values Throughout TRG Scheme",
    yscale = :log10,
    legend = false,
    )
ylims!(5e-3, 1.5e+0)
savefig(plt, "SVs TRG")

function compute_HOTRG_data_sing(ac_val, bc_val; ntruncdim=16, nmaxiter=nmaxiter)
    T = six_vertex_tensor(ac_val, bc_val, 1)  # sixvertex(a=ac_val, b=bc_val, c=1) # six_vertex_tensor(ac_val, bc_val, 1)
    scheme = HOTRG(T)
    data, sing = run_sing!(scheme, truncdim(ntruncdim), maxiter(nmaxiter); verbosity=0)
    return data, sing
end


_, sing2 = compute_HOTRG_data_sing(a, b; ntruncdim=16, nmaxiter=nmaxiter)
x2, y2 = extract_svals(sing2)
# Scatter plot
plt = scatter(x2, y2,
    xlabel = "Iteration",
    ylabel = "Singular value",
    title = "Singular Values Throughout HOTRG Scheme",
    yscale = :log10,
    legend = false,
    )
ylims!(5e-3, 1.5e+0)
savefig(plt, "SVs HOTRG")

function compute_BTRG_data_sing(ac_val, bc_val; ntruncdim=16, nmaxiter=nmaxiter)
    T = six_vertex_tensor(ac_val, bc_val, 1)  # sixvertex(a=ac_val, b=bc_val, c=1) # six_vertex_tensor(ac_val, bc_val, 1)
    scheme = BTRG(T)
    data, sing = run_sing!(scheme, truncdim(ntruncdim), maxiter(nmaxiter); verbosity=0)
    return data, sing
end

_, sing3 = compute_BTRG_data_sing(a, b; ntruncdim=16, nmaxiter=nmaxiter)
x3, y3 = extract_svals(sing3)
# Scatter plot
plt = scatter(x3, y3,
    xlabel = "Iteration",
    ylabel = "Singular value",
    title = "Singular Values Throughout BTRG Scheme",
    yscale = :log10,
    legend = false,
    )
ylims!(5e-3, 1.5e+0)
savefig(plt, "SVs BTRG")