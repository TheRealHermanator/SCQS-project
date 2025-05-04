using DelimitedFiles
using Plots
using TNRKit
using TensorKit
include("6VertexModel.jl")
include("ModifiedTRGflow.jl")

# Variables
nmaxiter = 25
select_iter = 21
as = LinRange(0.1, 1.0, 10)
b = 1.5 #such that we have a phase transition for a = 0.5 (delta = 1)
x = Float64[]
y = Float64[]

function compute_BTRG_data_sing(ac_val, bc_val; ntruncdim=16, nmaxiter=nmaxiter)
    T = six_vertex_tensor(ac_val, bc_val, 1)  # sixvertex(a=ac_val, b=bc_val, c=1) # six_vertex_tensor(ac_val, bc_val, 1)
    scheme = BTRG(T)
    data, sing = run_sing!(scheme, truncdim(ntruncdim), maxiter(nmaxiter); verbosity=1)
    return data, sing
end
# Iterate over different a's
for a in as
    data, sing = compute_BTRG_data_sing(a, b; ntruncdim=16, nmaxiter=nmaxiter)
    svals = sing[select_iter]
    maxval = maximum(svals)
    for val in svals
        push!(x, a)              # Store the iteration index
        push!(y, val / maxval)            # Store the singular value
    end
end

# Scatter plot
plt = scatter(x, y,
    xlabel = "a",
    ylabel = "Singular Values",
    title = "Singular Values In Different Phases (Iteration: $select_iter)",
    yscale = :log10,
    legend = false,
    )
line = LinRange(5e-3, 1.1e+0, 1000)
plot!(0.5 * ones(1000), line, color = :red)
annotate!(0.3, 1.5, "∆ > 1")
annotate!(0.5, 1.5, "∆ = 1")
annotate!(0.8, 1.5, "∆ < 1")
ylims!(5e-3, 2e+0)
savefig(plt, "SVs delta Iter$select_iter")
