using DelimitedFiles
using Plots
using TNRKit
using TensorKit
include("6VertexModel.jl")
include("ModifiedTRGflow.jl")

# Variables
nmaxiter = 20
as_1 = LinRange(0.1, 0.5, 5)
as_2 = LinRange(0.6, 1.0, 5)
b = 1.5 #such that we have a phase transition for a = 0.5 (delta = 1)
x = Float64[]
y = Float64[]

function compute_TRG_data_sing(ac_val, bc_val; ntruncdim=16, nmaxiter=nmaxiter)
    T = six_vertex_tensor(ac_val, bc_val, 1)  # sixvertex(a=ac_val, b=bc_val, c=1) # six_vertex_tensor(ac_val, bc_val, 1)
    scheme = BTRG(T)
    data, sing = run_sing!(scheme, truncdim(ntruncdim), maxiter(nmaxiter); verbosity=1)
    return data, sing
end


# Scatter plot
plt = scatter(
    xlabel = "\$a\$",
    ylabel = "Singular Values",
    title = "Singular Values In Different Phases at RG step $nmaxiter",
    yscale = :log10,
    legend = false,
    )

# Plotting the singular values for different a values
# First the ones with delta > 1 are plotted since it's easier to see the degeneracies
for a in as_1
    data, sing = compute_TRG_data_sing(a, b; ntruncdim=16, nmaxiter=nmaxiter)
    svals = sing[nmaxiter]
    maxval = maximum(svals)
    deg = 0
    for val in svals
        if abs(val-maxval)<1e-5
            deg +=1
        else
            scatter!((a, val/maxval), markercolor=:grey)
        end
        
    end
    if deg>1
        scatter!((a, 1), markercolor=:blue, markershape=:star5)
        annotate!((a-0.05, 1-0.2), "$deg")
    else
        scatter!((a, 1), markercolor=:blue)
    end
end
# Then the ones with delta < 1, we don't have to worry about degeneracies here
for a in as_2
    data, sing = compute_TRG_data_sing(a, b; ntruncdim=16, nmaxiter=nmaxiter)
    svals = sing[nmaxiter]
    maxval = maximum(svals)

    for val in svals
       scatter!((a, val/maxval), markercolor =:grey) 
    end
end
line = LinRange(5e-3, 1.1e+0, 1000)
plot!(0.5 * ones(1000), line, color = :red)
annotate!(0.3, 1.5, "\$∆ > 1\$")
annotate!(0.5, 1.5, "\$∆ = 1\$")
annotate!(0.8, 1.5, "\$∆ < 1\$")
ylims!(5e-3, 2e+0)
savefig(plt, "SVs delta iter $nmaxiter")