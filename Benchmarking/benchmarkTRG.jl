using TNRKit, TensorKit
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 30.0
BenchmarkTools.DEFAULT_PARAMETERS.samples = 1
function runTRG(symmetric=true)
    if symmetric 
        T = classical_ising_symmetric(ising_βc) # partition function of classical Ising model at the critical point
    else
        T = classical_ising(ising_βc) # partition function of classical Ising model at the critical point
    end
    scheme = TRG(T) # Bond-weighted TRG (excellent choice)
    data = run!(scheme, truncdim(16), maxiter(25), verbosity=0) # max bond-dimension of 16, for 25 iterations
    return data
end
symmetric = false
@btime runTRG(symmetric)