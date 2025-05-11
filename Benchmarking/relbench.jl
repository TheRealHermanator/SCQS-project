using TensorKit, Plots, TNRKit
x=1:30
y = addage = abs((f_onsager - f_onsager) / f_onsager)

TRGdata = []
BTRGdata = []
HOTRGdata = []
#First we do the benchmarking for TRG
for n in 1:30
    data = Any[]
    T = classical_ising_symmetric(ising_βc)
    scheme = TRG(T)
    data = run!(scheme, truncdim(16), maxiter(n), verbosity=0)
    global lnz = 0
    for (i, d) in enumerate(data)
        global lnz += log(d) * 2.0^(1 - i)
    end
    global f_ising = 0
    global f_ising += (lnz) * -1 / ising_βc
    addage1 = abs((f_ising - f_onsager) / f_onsager)
    push!(TRGdata, addage1)
end
scatter(x, TRGdata, label="TRG", markershape=:rtriangle,m=8, title="Convergence rate", yscale=:log10)
xlabel!("Number of iterations")
ylabel!("Relative error")
#Now we do the same for BTRG
for n in 1:30
    data = Any[]
    T = classical_ising_symmetric(ising_βc)
    scheme = BTRG(T)
    data = run!(scheme, truncdim(16), maxiter(n), verbosity=0)
    global lnz = 0
    for (i, d) in enumerate(data)
        global lnz += log(d) * 2.0^(1 - i)
    end
    global f_ising = 0
    global f_ising += (lnz) * -1 / ising_βc
    addage2 = abs((f_ising - f_onsager) / f_onsager)
    push!(BTRGdata, addage2)
end
scatter!(x, BTRGdata, label="BTRG",markershape=:star5, m = (5, :transparent))
#And for HOTRG

for n in 1:15
    data = Any[]
    T = classical_ising_symmetric(ising_βc)
    scheme = HOTRG(T)
    data = run!(scheme, truncdim(16), maxiter(n), verbosity=0)
    data2 = run!(scheme, truncdim(16), maxiter(2), verbosity=0)
    global lnz = 0
    for (i, d) in enumerate(data)
        global lnz += log(d) * 4.0^(1 - i)
    end
    global f_ising = 0
    global f_ising += (lnz) * -1 / ising_βc
    addage3 = abs((f_ising - f_onsager) / f_onsager)
    push!(HOTRGdata, addage3)
end
x3 = 1:2:30
print(x3, HOTRGdata)
scatter!(x3, HOTRGdata,markershape=:circle, label="HOTRG")


savefig("newfig.png")  