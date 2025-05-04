using DelimitedFiles
using Plots


function get_gradient(x_vals::Vector{<:Real}, y_vals::Vector{<:Real}, grid::AbstractMatrix{<:Real})
    nx, ny = length(x_vals), length(y_vals)
    grad = zeros(nx-1, ny-1)

    dx = x_vals[2] - x_vals[1]
    dy = y_vals[2] - y_vals[1]

    for i in 1:nx-1
        for j in 1:ny-1
            dfdx = (grid[i+1, j] - grid[i, j]) / dx
            dfdy = (grid[i, j+1] - grid[i, j]) / dy
            grad[i, j] = sqrt(dfdx^2 + dfdy^2)
        end
    end

    # Midpoints for plotting
    mid = v -> @. 0.5 * (v[1:end-1] + v[2:end])
    return mid(x_vals), mid(y_vals), grad
end


len = 101
typeplot = 2 # 0: ln(Z), 1: 1st derivative, 2: 2nd derivative

ac_vals = readdlm("Data/ac_vals_$len.csv")[:]
bc_vals = readdlm("Data/bc_vals_$len.csv")[:]
lnz_grid = readdlm("Data/lnZ_$len.csv")

# 1st derivative
ac_vals_1diff, bc_vals_1diff, lnz_grid_1diff = get_gradient(ac_vals, bc_vals, lnz_grid)

# 2nd derivative
ac_vals_2diff, bc_vals_2diff, lnz_grid_2diff = get_gradient(ac_vals_1diff, bc_vals_1diff, lnz_grid_1diff)


if typeplot == 0
    plt = heatmap(ac_vals, bc_vals, lnz_grid;
        xlabel="a/c", ylabel="b/c", title="ln(Z) Energy per Site",
        colorbar_title="ln(Z)", aspect_ratio=:equal)
    display(plt)
    savefig(plt, "Plots/Phase_len$(len)_dif$typeplot.png")
elseif typeplot == 1
    plt = heatmap(ac_vals_1diff, bc_vals_1diff, lnz_grid_1diff;
        xlabel="a/c", ylabel="b/c", title="1st Derivative ln(Z) Energy per Site",
        colorbar_title="|∇ln(Z)|", aspect_ratio=:equal)
    display(plt)
    savefig(plt, "Plots/Phase_len$(len)_dif$typeplot.png")
elseif typeplot == 2
    plt = heatmap(ac_vals_2diff, bc_vals_2diff, lnz_grid_2diff;
        xlabel="a/c", ylabel="b/c", title="2nd Derivative ln(Z) Energy per Site",
        colorbar_title="|∇²ln(Z)|", aspect_ratio=:equal)
    display(plt)
    savefig(plt, "Plots/Phase_len$(len)_dif$typeplot.png")
end