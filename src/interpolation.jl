"""
	interpolation()

The latent fields U and V are interpolated for each MCMC iteration to the grid cells containing no observation using the Gaussian Markov random field.

"""
function interpolation(C::Chains, data::DataStructure)

    S = data.S
    S̄ = data.S̄
    G = data.G

    n = length(data.Y)

    u = dropdims(C[:,["u[$i]" for i=1:n],1].value, dims=3)
    κ₁ = vec(C[:,"κ₁",1].value)

    v = dropdims(C[:,["v[$i]" for i=1:n],1].value, dims=3)
    κ₂ = vec(C[:,"κ₂",1].value)

    β₁ = vec(C[:,"β₁",1].value)
    β₂ = vec(C[:,"β₂",1].value)
    ξ = vec(C[:,"ξ",1].value)

    niter = length(κ₁)


    Waa = G.W[S̄.V,S̄.V]
    Wab = G.W[S̄.V,S.V]

    # number of grid cells
    m = prod(G.gridSize)

    U = Array{Float64}(undef, m, niter)
    V = Array{Float64}(undef, m, niter)

    F₀ = cholesky(Waa)

    @showprogress for j=1:niter

        U[S.V,j] = u[j,:]     # on connait U aux points des stations
        V[S.V,j] = v[j,:]

        b = -κ₁[j] * Wab * U[ S.V, j ]       # terme relié à la moyenne \mu
        Q = κ₁[j] * Waa                      # precision matrix (?)
        U[ S̄.V , j] = gmrfsample(b, Q, F₀)

        b = -κ₂[j] * Wab * V[ S.V, j]
        Q = κ₂[j] * Waa
        V[ S̄.V , j] = gmrfsample(b, Q, F₀)

    end

    parmnames = vcat(["u[$i]" for i=1:m], ["v[$i]" for i=1:m], "κ₁", "β₁", "κ₂", "β₂","ξ")
    res = vcat(U, V, κ₁', β₁', κ₂', β₂', ξ')
    completeChain = MambaLite.Chains(collect(res'), names=parmnames)

    return completeChain

end

function gmrfsample(b::Vector{<:Real},Q::AbstractArray{<:Real})

    F = cholesky(Q)

    μ = F\b

    z = randn(length(b))

    # Pivoting is on by default in SuiteSparse (https://julialinearalgebra.github.io/SuiteSparse.jl/latest/cholmod.html)
    v = F.UP\z

    y = μ + v

    return y

end

function gmrfsample(b::Vector{<:Real},Q::AbstractArray{<:Real}, F₀::SuiteSparse.CHOLMOD.Factor{Float64})

    F = cholesky!(F₀, Q)

    μ = F\b

    z = randn(length(b))

    # Pivoting is on by default in SuiteSparse (https://julialinearalgebra.github.io/SuiteSparse.jl/latest/cholmod.html)
    v = F.UP\z

    y = μ + v

    return y

end