using QuadGK, TensorKit, TNRKit, Plots
using LaTeXStrings
function onsager_free_energy(beta; J=1.0)
    if beta <= 0
        return NaN  # or error
    end
    a = sinh(2*beta * J)
    c = cosh(2*beta * J)
    
    integrand(kx, ky) = log(abs(c^2 - a * (cos(kx) + cos(ky))))

    # Nested integral
    integral, _ = QuadGK.quadgk(kx -> QuadGK.quadgk(ky -> integrand(kx, ky), -π, π)[1], -π, π)

    f = -(log(2) + integral / (8π^2)) / beta
    return f
end
BTRGdata = []
TRGdata = []
HOTRGdata = []
beta = range(0.4, 0.5, length=50)
#BTRG
for T in beta
    ising = classical_ising_symmetric(T)
    scheme = BTRG(ising)
    data = run!(scheme, truncdim(16), maxiter(22), verbosity=0)
    global lnz = 0
    for (i, d) in enumerate(data)
        global lnz += log(d) * 2.0^(1 - i)
    end
    global f_ising = 0
    global f_ising += (lnz) * -1 / T
    f_onsager = onsager_free_energy(T; J=1.0)
    addage1 = abs((f_ising - f_onsager)/ f_onsager)
    push!(BTRGdata, addage1)
end
scatter(beta, BTRGdata, label="BTRG", m=4, markershape=:star5, title=L"Accuracy around $\beta_c$", yscale=:log10)
#TRG
for T in beta
    ising = classical_ising_symmetric(T)
    scheme = TRG(ising)
    data = run!(scheme, truncdim(16), maxiter(22), verbosity=0)
    global lnz = 0
    for (i, d) in enumerate(data)
        global lnz += log(d) * 2.0^(1 - i)
    end
    global f_ising = 0
    global f_ising += (lnz) * -1 / T
    f_onsager = onsager_free_energy(T; J=1.0)
    addage2 = abs((f_ising - f_onsager)/ f_onsager)
    push!(TRGdata, addage2)
end
scatter!(beta, TRGdata, markershape=:rtriangle, label="TRG", m=4)
#HOTRG
for T in beta
    ising = classical_ising_symmetric(T)
    scheme = HOTRG(ising)
    data = run!(scheme, truncdim(16), maxiter(22), verbosity=0)
    global lnz = 0
    for (i, d) in enumerate(data)
        global lnz += log(d) * 4.0^(1 - i)
    end
    global f_ising = 0
    global f_ising += (lnz) * -1 / T
    f_onsager = onsager_free_energy(T; J=1.0)
    addage3 = abs((f_ising - f_onsager)/ f_onsager)
    push!(HOTRGdata, addage3)
end
scatter!(beta, HOTRGdata, label="HOTRG", markershape=:circle, m=4)
vline!([1/2.269], label=L"Critical $\beta$")
xlabel!(L"$\beta$")
ylabel!("Relative error")
savefig("benchinT.png")