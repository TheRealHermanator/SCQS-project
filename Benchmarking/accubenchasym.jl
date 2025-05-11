using TensorKit, Plots, TNRKit
x=1:30
y = f_onsager
plot(x, fill(y, length(x)), label="Onsager", title="Convergence rate with asymmetric Ising")
xlabel!("Number of max iterations")
ylabel!("Free energy")
TRGdata = []
BTRGdata = []
HOTRGdata = []
#First we do the benchmarking for TRG
for n in 1:30
    data = Any[]
    T = classical_ising(ising_βc)
    scheme = TRG(T)
    data = run!(scheme, truncdim(16), maxiter(n), verbosity=0)
    global lnz = 0
    for (i, d) in enumerate(data)
        global lnz += log(d) * 2.0^(1 - i)
    end
    global f_ising = 0
    global f_ising += (lnz) * -1 / ising_βc
    push!(TRGdata, f_ising)
end
scatter!(x, TRGdata, label="TRG", m=8)
#Now we do the same for BTRG
for n in 1:30
    data = Any[]
    T = classical_ising(ising_βc)
    scheme = BTRG(T)
    data = run!(scheme, truncdim(16), maxiter(n), verbosity=0)
    global lnz = 0
    for (i, d) in enumerate(data)
        global lnz += log(d) * 2.0^(1 - i)
    end
    global f_ising = 0
    global f_ising += (lnz) * -1 / ising_βc
    push!(BTRGdata, f_ising)
end
scatter!(x, BTRGdata, label="BTRG", m = (5, :transparent))
#And for HOTRG
for n in 1:30
    data = Any[]
    T = classical_ising(ising_βc)
    scheme = HOTRG(T)
    data = run!(scheme, truncdim(16), maxiter(n), verbosity=0)
    data2 = run!(scheme, truncdim(16), maxiter(2), verbosity=0)
    global lnz = 0
    for (i, d) in enumerate(data)
        global lnz += log(d) * 2.0^(1 - i)
    end
    global f_ising = 0
    global f_ising += (lnz) * -1 / ising_βc
    push!(HOTRGdata, f_ising)
end
scatter!(x, HOTRGdata, label="HOTRG")

savefig("accubenchasym.png")  