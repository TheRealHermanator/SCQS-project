using TensorKit


function six_vertex_tensor(a::Real, b::Real, c::Real)
"""
Create a TRG family compatible 6-vertex model tensor with weights a, b, c.
Returns a TensorMap of shape (left, down ← up, right).
"""
    V = ℂ^2
    T = zeros(Float64, V ⊗ V ← V ⊗ V)

    # Arrow convention: 0 = out, 1 = in
    # Index convention (clockwise): (u, r, d, l)
    # Image source: wikipedia on Ice type model
    configs = [
        ((0, 0, 1, 1), a),  # Vertex 1
        ((1, 1, 0, 0), a),  # Vertex 2
        ((1, 0, 0, 1), b),  # Vertex 3
        ((0, 1, 1, 0), b),  # Vertex 4
        ((0, 1, 0, 1), c),  # Vertex 5
        ((1, 0, 1, 0), c),  # Vertex 6
    ]

    for ((u, r, d, l), weight) in configs
        # Julia starts with index 1, so we add 1. Stupid, I know. Do like Python bro
        T[l + 1, d + 1, u + 1, r + 1] = weight
    end

    return T
end

six_vertex_tensor(1, 2, 3)