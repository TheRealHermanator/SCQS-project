using DelimitedFiles
using Plots
using TNRKit
using TensorKit
include("6VertexModel.jl")
include("ModifiedTRGflow.jl")


# Variables
nmaxiter = 20 #since the SV's become rubbish after 25th iteration
as = [0.25, 0.25, 1, 1.5]
bs = [0.25, 1.5, 1, 0.25]
phases = ["III", "I", "IV", "II"]

function compute_TRG_data_cft(ac_val, bc_val; ntruncdim=16, nmaxiter=nmaxiter)
    T = six_vertex_tensor(ac_val, bc_val, 1)  
    scheme = BTRG(T, finalize=finalize_cftdata!)
    data = run!(scheme, truncdim(ntruncdim), maxiter(nmaxiter); verbosity=1)
    return data
end
for i in range(1,4)
    a = as[i]
    b = bs[i]
    fil = phases[i]
    delta = (a^2 + b^2 - 1)/(2*a*b)
    data = compute_TRG_data_cft(a, b, ntruncdim = 16, nmaxiter = nmaxiter)
    xs = Int[]
    ys = Float64[]
    for (j, vec) in enumerate(data)
        for val in vec
            push!(xs, j)
            push!(ys, exp(-2π*val))
        end
    end
    
    # Plot as scatter
    scatter(xs, ys, xlabel="Renormalisation Step", ylabel="λn(L)/λ0(L)", title="Eigenvalues for Δ=$delta in phase $fil")
    savefig("CFT_$fil.png")      
end
