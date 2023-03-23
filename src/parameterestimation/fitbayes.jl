function mcmc(data::DataStructure; niter::Int=5000, warmup::Int=2000, stepsize::Array{Float64}=[.1, .15, .02], thin::Int64=5, returnUV::Bool=false)
    
    ini = initialize_rf(data)
    
    return mcmc(data, ini, niter=niter, warmup=warmup, stepsize=stepsize, thin=thin, returnUV=returnUV)
    
end


function mcmc(data::DataStructure, initial_values::Tuple{Vector{Float64}, Vector{Float64}, Float64, Vector{Float64}, Vector{Float64}, Real, Vector{Float64}, Vector{Float64}, Real}; niter::Int=5000, warmup::Int=2000, stepsize::Array{Float64}=[.1, .15, .02], thin::Int64=5, returnUV::Bool=false)
    
    n = length(data.Y)
    p₁ = size(data.X₁ᵢ, 2)
    p₂ = size(data.X₂ᵢ, 2)
    
    # Pre-initialize outputs 
    u = Array{Float64}(undef, n, niter)
    β₁ = Array{Float64}(undef, p₁, niter)
    κ₁ = Array{Float64}(undef, niter)

    v = Array{Float64}(undef, n, niter)
    β₂ = Array{Float64}(undef, p₂, niter)
    κ₂ = Array{Float64}(undef, niter)

    ξ = Array{Float64}(undef, niter)
    
    # Initial values
    μ₀, σ₀, ξ[1], U, β₁[:,1], κ₁[1], V,β₂[:,1], κ₂[1] = initial_values
    u[:,1] = U[data.S.V]
    v[:,1] = V[data.S.V]
    
    F₁ = GMRF.iGMRF(data.G, 1, κ₁[1])
    F₂ = GMRF.iGMRF(data.G, 1, κ₂[1])

    δ = randn(n)

    @showprogress for j=2:niter

        # Generate a candidate for {Uᵢ : i ∈ S}
        Ũ = u[:,j-1] + stepsize[1]*randn!(δ)

        # Computing the corresponding candidates for μ at the grid cells containing the observations
        μ̃ = exp.(data.X₁ᵢ*β₁[:,j-1] + Ũ);

        # Evaluate the log-likelihood ratio at the data level between the candidates and the present state
        logL = datalevel_loglike.(data.Y, μ̃, σ₀, ξ[j-1]) - datalevel_loglike.(data.Y, μ₀, σ₀, ξ[j-1])

        # Updating the latent field U
        latentgmrf_update!(U, F₁, data.S, data.S̄, Ũ, logL)

        u[:,j] = U[data.S.V]

        x = data.X₁ᵢ*β₁[:, j-1]+u[:,j]
        μ₀ = exp.(x)
        β₁[:, j] = data.X₁ᵢ \ ( x .- mean(U[data.S.V]) )  

        # Sampling the new precision parameter
        κ₁[j] = latentfieldprecision_sample(F₁, U)

        # Updating the iGMRF object with the new precision
        F₁ = GMRF.iGMRF(data.G, 1, κ₁[j])


        # Generate a candidate for {Vᵢ : i ∈ S}
        Ṽ = v[:,j-1] + stepsize[2]*randn!(δ)

        # Computing the corresponding candidates for σ at the grid cells containing the observations
        σ̃ = exp.(data.X₂ᵢ*β₂[:, j-1] + Ṽ)

        # Evaluate the log-likelihood ratio at the data level between the candidates and the present state
        logL = datalevel_loglike.(data.Y, μ₀, σ̃, ξ[j-1]) - datalevel_loglike.(data.Y, μ₀, σ₀, ξ[j-1])

        # Updating the latent field V
        latentgmrf_update!(V, F₂, data.S, data.S̄, Ṽ, logL)

        v[:,j] = V[data.S.V]

        x = data.X₂ᵢ*β₂[:,j-1]+v[:,j]
        σ₀ = exp.(x)
        β₂[:, j] = data.X₂ᵢ \ ( x .- mean(V[data.S.V]) )
        

        # Sampling the new precision parameter
        κ₂[j] = latentfieldprecision_sample(F₂, V)

        # Updating the iGMRF object with the new precision
        F₂ = GMRF.iGMRF(data.G, 1, κ₂[j])

        ξ̃ = ξ[j-1] + stepsize[3]*randn()
        logL = sum(datalevel_loglike.(data.Y, μ₀, σ₀, ξ̃) - datalevel_loglike.(data.Y, μ₀, σ₀, ξ[j-1]))
        if logL > log(rand())
            ξ[j] = ξ̃
        else
            ξ[j] = ξ[j-1]
        end

    end

    parmnames = vcat(["u[$i]" for i=1:n], ["v[$i]" for i=1:n], "κ₁", ["β₁[$i]" for i=1:p₁], "κ₂", ["β₂[$i]" for i=1:p₂], "ξ")
    res = vcat(u, v, κ₁', β₁, κ₂', β₂, ξ')
    C = MambaLite.Chains(collect(res'), names=parmnames)
    C = C[warmup+1:thin:end,:,:]
    
    if returnUV
        return C, U, V
    else
        return C
    end
    
end


function mcmc_valid(data::DataStructure, initial_values::Tuple{Vector{Float64}, Vector{Float64}, Float64, Vector{Float64}, Vector{Float64}, Real, Vector{Float64}, Vector{Float64}, Real}, gridpoint::Int64, X₁::Matrix{Float64}, X₂::Matrix{Float64}; niter::Int=5000, warmup::Int=2000, stepsize::Array{Float64}=[.1, .15, .02], thin::Int64=5) # u, v, xi
    
    n = length(data.Y)
    p₁ = size(data.X₁ᵢ, 2)
    p₂ = size(data.X₂ᵢ, 2)
    
    # Pre-initialize outputs 
    u = Array{Float64}(undef, n, niter)
    β₁ = Array{Float64}(undef, p₁, niter)
    κ₁ = Array{Float64}(undef, niter)

    v = Array{Float64}(undef, n, niter)
    β₂ = Array{Float64}(undef, p₂, niter)
    κ₂ = Array{Float64}(undef, niter)

    ξ = Array{Float64}(undef, niter)
    μᵥ = Array{Float64}(undef, niter)
    σᵥ = Array{Float64}(undef, niter)
    
    # Initial values
    μ₀, σ₀, ξ[1], U, β₁[:,1], κ₁[1], V,β₂[:,1], κ₂[1] = initial_values
    u[:,1] = U[data.S.V]
    v[:,1] = V[data.S.V]
    
    μᵥ[1] = exp((X₁*β₁[:,1])[1] + U[gridpoint])
    σᵥ[1] = exp((X₂*β₂[:,1])[1] + V[gridpoint])
    
    F₁ = GMRF.iGMRF(data.G, 1, κ₁[1])
    F₂ = GMRF.iGMRF(data.G, 1, κ₂[1])

    δ = randn(n)

    @showprogress for j=2:niter

        # Generate a candidate for {Uᵢ : i ∈ S}
        Ũ = u[:,j-1] + stepsize[1]*randn!(δ)

        # Computing the corresponding candidates for μ at the grid cells containing the observations
        μ̃ = exp.(data.X₁ᵢ*β₁[:,j-1] + Ũ);

        # Evaluate the log-likelihood ratio at the data level between the candidates and the present state
        logL = datalevel_loglike.(data.Y, μ̃, σ₀, ξ[j-1]) - datalevel_loglike.(data.Y, μ₀, σ₀, ξ[j-1])

        # Updating the latent field U
        latentgmrf_update!(U, F₁, data.S, data.S̄, Ũ, logL)

        u[:,j] = U[data.S.V]

        x = data.X₁ᵢ*β₁[:, j-1]+u[:,j]
        μ₀ = exp.(x)
        β₁[:, j] = data.X₁ᵢ \ ( x .- mean(U[data.S.V]) )  

        # Sampling the new precision parameter
        κ₁[j] = latentfieldprecision_sample(F₁, U)

        # Updating the iGMRF object with the new precision
        F₁ = GMRF.iGMRF(data.G, 1, κ₁[j])


        # Generate a candidate for {Vᵢ : i ∈ S}
        Ṽ = v[:,j-1] + stepsize[2]*randn!(δ)

        # Computing the corresponding candidates for σ at the grid cells containing the observations
        σ̃ = exp.(data.X₂ᵢ*β₂[:, j-1] + Ṽ)

        # Evaluate the log-likelihood ratio at the data level between the candidates and the present state
        logL = datalevel_loglike.(data.Y, μ₀, σ̃, ξ[j-1]) - datalevel_loglike.(data.Y, μ₀, σ₀, ξ[j-1])

        # Updating the latent field V
        latentgmrf_update!(V, F₂, data.S, data.S̄, Ṽ, logL)

        v[:,j] = V[data.S.V]

        x = data.X₂ᵢ*β₂[:,j-1]+v[:,j]
        σ₀ = exp.(x)
        β₂[:, j] = data.X₂ᵢ \ ( x .- mean(V[data.S.V]) )
        

        # Sampling the new precision parameter
        κ₂[j] = latentfieldprecision_sample(F₂, V)

        # Updating the iGMRF object with the new precision
        F₂ = GMRF.iGMRF(data.G, 1, κ₂[j])

        ξ̃ = ξ[j-1] + stepsize[3]*randn()
        logL = sum(datalevel_loglike.(data.Y, μ₀, σ₀, ξ̃) - datalevel_loglike.(data.Y, μ₀, σ₀, ξ[j-1]))
        if logL > log(rand())
            ξ[j] = ξ̃
        else
            ξ[j] = ξ[j-1]
        end
        
        μᵥ[j] = exp((X₁*β₁[:,j])[1] + U[gridpoint])
        σᵥ[j] = exp((X₂*β₂[:,j])[1] + V[gridpoint])

    end

    parmnames = vcat(["u[$i]" for i=1:n], ["v[$i]" for i=1:n], "κ₁", ["β₁[$i]" for i=1:p₁], "κ₂", ["β₂[$i]" for i=1:p₂], "ξ")
    res = vcat(u, v, κ₁', β₁, κ₂', β₂, ξ')
    C = MambaLite.Chains(collect(res'), names=parmnames)
    C = C[warmup+1:thin:end,:,:]
    
    parmnames2 = vcat("μ", "σ", "ξ")
    res2 = vcat(μᵥ', σᵥ', ξ')
    C2 = MambaLite.Chains(collect(res2'), names=parmnames2)
    C2 = C2[warmup+1:thin:end,:,:]
    
    return C, C2
    
end


function latentgmrf_update!(x::Vector{<:Real}, F::GMRF.iGMRF, S::GridPointStructure, S̄::GridPointStructure, x̃::Vector{<:Real}, logL::Vector{<:Real})


    for j in eachindex(S.CondIndSubset)

        V = S.CondIndSubset[j]

        ind = S.CondIndIndex[j]

        latentcondgmrf_update!(x, F, V, x̃[ind], logL[ind])

    end

    for j in eachindex(S̄.CondIndSubset)

        V = S̄.CondIndSubset[j]

        latentcondgmrf_update!(x, F, V)

    end
end

function latentcondgmrf_update!(x::Vector{<:Real}, F::GMRF.iGMRF, V::Vector{<:Real})

    #@timeit "fullcond" pd = GMRF.fullconditionals(F,x)[V]
    pd = InterpolationIDF.fullconditionals(F,x)[V]

    x̃ = rand.(pd)

    setindex!(x,x̃,V)

end

function latentcondgmrf_update!(x::Vector{<:Real}, F::GMRF.iGMRF, V::Vector{<:Real}, x̃::Vector{<:Real}, logL::Vector{<:Real})

    u = rand(length(V))

    #@timeit "fullcond" pd = GMRF.fullconditionals(F,x)[V]
    pd = InterpolationIDF.fullconditionals(F,x)[V]
    
    lf = logpdf.(pd,x̃) - logpdf.(pd,x[V])

    lr = logL + lf

    ind = lr .> log.(u)

    setindex!(x,x̃[ind],V[ind])

end

function fullconditionals(F::iGMRF, y::Vector{<:Real})::Vector{NormalCanon}

    κ = F.κ

    W̄ = F.G.W̄
    W = F.G.W

    #@timeit "Q" Q = κ * Array(diag(F.G.W))
    Q = κ * diag(W).nzval
    
    h = -κ*(W̄*y)

    pd = NormalCanon.(h,Q)

    return pd

end

function latentfieldprecision_sample(F::GMRF.iGMRF, x::Vector{<:Real})

    m = prod(F.G.gridSize)
    k = F.rankDeficiency
    W = F.G.W

    b = dot(x, W*x)

    # Informative prior Gamma(1,1/100)
    pd = Gamma( (m - k)/2 + 1 , 1/(b/2 + 1/100) )

    # Improper prior
    # pd = Gamma( (m - k)/2 , 2/b)

    κ = rand(pd)

    return κ

end

"""
    iGMRFupdate(F::GMRF.iGMRF,κ::Real)

Return an iGMRF object similar to F but with the precision parameter equals to κ.

### Arguments
- `F` : An iGMRF object. See the GMRF.jl package for more details.
- `κ` : The new iGMRF precision parameter. It should be a positive real number.

### Details

The function uses the GMRF.jl package.

### Examples

\```
 julia> F = iGMRFupdate(F,10)
\```

"""
function iGMRFupdate(F::GMRF.iGMRF,κ::Real)

    @assert κ>0

    F̃ = GMRF.iGMRF(F.G, F.rankDeficiency, κ, F.W, F.W̄)

    return F̃

end


"""
    regressioncoefficient_sample(F::GMRF.iGMRF, Covariate::Vector{<:Real}, U::Vector{<:Real}, β::Real)

Sample the regression coefficients of the latent iGMRF mean using the complete conditional distribution.

### Arguments
- `F`: An iGMRF object. See the GMRF.jl package for more details.
- `Covariate`: The covariate values at every vertices.
- `x`: The current state of the iGMRF field.

### Details

The function uses the GMRF.jl package.

### Examples

\```
 julia> κ = latentfieldprecision_sample(F, x)
\```

"""
function regressioncoefficient_sample(F::GMRF.iGMRF, Covariate::Vector{<:Real}, x::Vector{<:Real})

    κ = F.κ
    W = F.G.W

#     x = Covariate*β + U

    h = κ*Covariate'*(W*x)

    pd = NormalCanon(h,κ*Qₓ)

    β̃ = rand(pd)

    return β̃

end

# diagnostic tools

function accrate(θ::Array{<:Real})

    d = θ[:,2:end] - θ[:,1:end-1]
    rate = 1 .- mean(d .≈ 0)

    return rate

end