#### initialization

function initialize_β(y::Vector{<:Real}, X::Array{<:Real})

    # removing the offset
    X̃ = X .- mean(X, dims=1)
    ỹ = y .- mean(y)

    β = X̃\ỹ

    return β

end

function initialize_U(Uᵢ::Vector{<:Real}, W::SparseMatrixCSC{Int64,Int64}, κ::Int64, S::GridPointStructure, S̄::GridPointStructure)

    m = size(W, 1)

    U = Array{Float64}(undef, m)
    U[S.V] = Uᵢ

    Waa = W[S̄.V,S̄.V]
    Wab = W[S̄.V,S.V]

    b = -κ*Wab*Uᵢ
    Q = κ*Waa
    U[S̄.V] = gmrfsample(b,Q)

    return U

end

function initialize_rf(data::DataStructure)

#     Y = data.Y
#     X₁ᵢ = data.X₁ᵢ
#     X₂ᵢ = data.X₂ᵢ
#     G = data.G
#     S = data.S
#     S̄ = data.S̄

    # Initial values
    κ₁ = 100
    κ₂ = 80

    n = length(data.Y)

    μ₀ = Array{Float64}(undef, n)
    ϕ₀ = Array{Float64}(undef, n)
    ξ₀ = 0.0

    for i=1:n
        fd = Extremes.gumbelfitpwm(data.Y[i])
        μ₀[i] = fd.θ̂[1]
        ϕ₀[i] = fd.θ̂[2]
    end

    β₁ = initialize_β(log.(μ₀), data.X₁ᵢ)
    Uᵢ = log.(μ₀) - data.X₁ᵢ*β₁
    U = initialize_U(Uᵢ, data.G.W, κ₁, data.S, data.S̄)

    β₂ = initialize_β(ϕ₀, data.X₂ᵢ)
    Vᵢ = ϕ₀ - data.X₂ᵢ*β₂
    V = initialize_U(Vᵢ, data.G.W, κ₂, data.S, data.S̄)

#     ini = Dict(:μ₀ => μ₀, :σ₀ => exp.(ϕ₀), :ξ₀ => ξ₀,
#         :U => U, :β₁ => β₁, :κ₁ => κ₁,
#         :V => V, :β₂ => β₂, :κ₂ => κ₂)

#     return ini
    
    return μ₀, exp.(ϕ₀), ξ₀, U, β₁, κ₁, V, β₂, κ₂

end