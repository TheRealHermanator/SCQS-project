Welcome to our repository for the Tensor Network Renormalization (TNR) project, developed as part of the elective course Strongly Correlated Quantum Systems. The project explores both algorithmic benchmarking and physical results using TNR techniques.

# 📌 Main objectives

1)	Benchmark three TNR schemes: TRG, BTRG, and HOTRG, using the classical Ising model.
2)	Compute the phase diagram of the six-vertex model using tensor network methods.


# 📂 Repository Structure

🔬 **Benchmarking**
Contains code for benchmarking TRG, BTRG, and HOTRG with respect to:

-	Efficiency (runtime, memory usage)
- Convergence rate
-	Behavior near critical points
-	Impact of truncation dimension

🧊 **Partition function phase diagram**

Code to compute and visualize the partition function and its derivatives over a parameter grid for the six-vertex model. The result reveals second-order phase transitions via gradients and Laplacians.

🌀 **Partition function different phases**

Tools to analyze different phases of the six-vertex model through:
•	Singular values for different phases 
•	Eigenvalues of the transfer matrix, highlighting gapped vs. critical phases

# 🧰 Dependencies

This repository makes extensive use of:
-	[TNRKit.jl](https://github.com/VictorVanthilt/TNRKit.jl)
-	[TensorKit.jl](https://github.com/Jutho/TensorKit.jl)


A star is appreciated! :)
