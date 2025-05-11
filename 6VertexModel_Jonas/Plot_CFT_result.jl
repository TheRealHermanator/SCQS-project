using DelimitedFiles
using Plots
using TNRKit
using TensorKit
using LaTeXStrings
include("6VertexModel.jl")
include("ModifiedTRGflow.jl")
# Variables
nmaxiter = 20 #since the SV's become rubbish after 25th iteration
as = [0.25, 0.25, 1, 1.5]
bs = [0.25, 1.5, 1, 0.25]
phases = ["III", "I", "IV", "II"]

#Simple function that performs the CFT calculations
"""
    compute_TRG_data_cft(ac_val, bc_val; ntruncdim=16, nmaxiter=nmaxiter)

Compute the TRG (Tensor Renormalization Group) data for given `ac_val` and `bc_val` parameters.
Uses the six-vertex tensor and BTRG scheme with CFT data finalization.

# Arguments
- `ac_val`: Parameter a/c for the six-vertex model.
- `bc_val`: Parameter b/c for the six-vertex model.
- `ntruncdim`: (optional) Truncation dimension for the TRG (default: 16).
- `nmaxiter`: (optional) Maximum number of iterations (default: `nmaxiter`).

# Returns
- `data`: The result of the TRG flow, typically a vector of eigenvalue vectors per iteration.
"""
function compute_TRG_data_cft(ac_val, bc_val; ntruncdim=16, nmaxiter=nmaxiter)
    T = six_vertex_tensor(ac_val, bc_val, 1)  
    scheme = BTRG(T, finalize=finalize_cftdata!)
    data = run!(scheme, truncdim(ntruncdim), maxiter(nmaxiter); verbosity=1)
    return data
end

"""
    degeneracy(vec)

Check a vector of eigenvalues for degeneracies and classify minima and maxima.

# Arguments
- `vec`: Vector of eigenvalues (or related quantities).

# Returns
- For `length(vec) == 0`: `nothing`
- For `length(vec) == 1` or `2`: Tuple of (min, max, min_deg, max_deg, rest)
- For `length(vec) > 2`: Tuple of (min, max, min_deg, max_deg, rest), where `rest` contains non-degenerate values.
"""
function degeneracy(vec)
    l = length(vec)
    if l == 0
        return nothing  # or throw error / custom behavior
    elseif l == 1
        λ = exp(-2π * vec[1])
        if λ >0.5  #then it's a maximum at that iteration
            return -1, λ, 0, 1, Float64[]#
        else    #then it's a minimum at that iteration
            return λ, -1, 1, 0, Float64[]
        end
    elseif l == 2
        vec_sort = sort(vec)
        min_val, max_val = vec_sort[1], vec_sort[2]
        l_min, l_max = exp(-2π * min_val), exp(-2π * max_val)

        if abs(l_min - l_max) < 1e-8
            # They're effectively the same value
            if l_min >0.5 #then it's a maximum at that iteration
                return -1, l_min, 0, 2, Float64[] #-1 to scatter outside the view
            else #then it's a minimum at that iteration
                return l_min, -1, 2, 0, Float64[]
            end
        else #they are different
            return l_min, l_max, 1, 1, Float64[]
        end
    else #for length >2 we can just work as follows
        vec_sort = sort(vec)
        min_val, max_val = vec_sort[1], vec_sort[end]
        l_max, l_min = exp(-2π * min_val), exp(-2π * max_val)

        deg_min, deg_max = 1, 1  # Start at 1 to count min/max themselves
        rest = Float64[]

        for i in 2:l-1  # skip vec_sort[1] and vec_sort[end]
            lambda_i = exp(-2π * vec_sort[i])
            if abs(lambda_i - l_min) < 1e-3 #then it's one of the degenerate minima
                deg_min += 1
            elseif abs(lambda_i - l_max) <  1e-3 #then it's one of the degenerate maxima
                deg_max += 1
            else #it's a separate point and can be plotted as such
                push!(rest, lambda_i)
            end
        end

        return l_min, l_max, deg_min, deg_max, rest
    end
end


 """
 The code below generates scatter plots for the transfer matrix eigenvalues
 for different phases of the six-vertex model. 
 """
 
for i in range(1,4)
    a = as[i]
    b = bs[i]
    fil = phases[i]
    delta = (a^2 + b^2 - 1)/(2*a*b)
    data = compute_TRG_data_cft(a, b, ntruncdim = 16, nmaxiter = nmaxiter)
    plt = scatter(xlabel="Renormalisation Step", 
    ylabel="\$\\frac{λ_n}{λ_0}\$", 
    title="Transfer matrix eigenvalues for \$Δ\$=$(delta) in phase $(fil) ", 
    legend=false)
    ylims!(-0.5, 1.5)
    for (j, vec) in enumerate(data)
        # Check for degeneracies
        # and classify min/max  
        l_min, l_max, deg_min, deg_max, rest = degeneracy( vec)
        if deg_min>1
            scatter!(plt,(j, l_min), markershape=:star5, markercolor=:red)
            annotate!(j, l_min+0.1, "$deg_min")
        else
            scatter!(plt, (j, l_min), markercolor=:red)
        end
        if deg_max>1
            scatter!(plt, (j, l_max),markershape=:star5, markercolor=:blue)
            annotate!(j, l_max -0.1, "$deg_max", markercolor=:blue)
        else
            scatter!(plt, (j, l_max), markercolor=:blue)
        end
        for val in rest
            scatter!(plt, (j, val), markercolor=:black)
        end
    end
    savefig("CFT_$fil.png")      
end
