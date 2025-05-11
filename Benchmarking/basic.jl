T_array = range(0.5, 5.0, length=20)
using QuadGK, TNRKit, TensorKit
function onsager_free_energy(T; J=1.0)
    if T <= 0
        return NaN  # or error
    end
    β = 1 / T
    a = sinh(2β * J)
    c = cosh(2β * J)
    
    integrand(kx, ky) = log(abs(c^2 - a * (cos(kx) + cos(ky))))

    # Nested integral
    integral, _ = QuadGK.quadgk(kx -> QuadGK.quadgk(ky -> integrand(kx, ky), -π, π)[1], -π, π)

    f = -(log(2) + integral / (8π^2)) / β
    return f
end

f_array = [onsager_free_energy(T) for T in T_array]
using Plots
scatter(T_array, f_array, xlabel="T", ylabel="Free Energy", title="Onsager Free Energy", legend=false)
vline!([2.269])
savefig("test.png")
a = onsager_free_energy(1/ising_βc)
b = f_onsager
print(a,b)