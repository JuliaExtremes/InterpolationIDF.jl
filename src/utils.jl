"""
	create_datastructure(G::GMRF.GridStructure, station_list::DataFrame, m₁::Int64, m₂::Int64, covariate::Vector{Float64})

"""
function create_datastructure(G::GMRF.GridStructure, station_list::DataFrame, m₁::Int64, m₂::Int64, covariate::Vector{Float64})
    
    # Station data
    Y = station_list.Data
    n = length(Y)
    m = m₁*m₂

    V = station_list[:,:GridCell]
    condIndIndex = [ findall(in(intersect(G.condIndSubset[j],V)),V) for j=1:length(G.condIndSubset)  ]
    condIndSubset = [ V[condIndIndex[j]] for j=1:length(condIndIndex)  ]
    S = GridPointStructure(V, condIndSubset, condIndIndex)

    V̄ = setdiff(1:m,V)
    condIndIndex = [ findall(in(intersect(G.condIndSubset[j],V̄)),V̄) for j=1:length(G.condIndSubset) ]
    condIndSubset = [ V̄[condIndIndex[j]] for j=1:length(condIndIndex) ]
    S̄ = GridPointStructure(V̄, condIndSubset, condIndIndex)

    # Covariates 
    X₁ = log.(covariate)
    X₁ᵢ = X₁[S.V]

    X₂ = log.(covariate)
    X₂ᵢ = X₂[S.V]

    return DataStructure(Y, X₁ᵢ, X₂ᵢ, G, S, S̄)
    
end

"""
	create_datastructure(G::GMRF.GridStructure, station_list::DataFrame, m₁::Int64, m₂::Int64, spatial_cov::Vector{Float64}, local_cov::Vector{Float64})

"""
function create_datastructure(G::GMRF.GridStructure, station_list::DataFrame, m₁::Int64, m₂::Int64, spatial_cov::Vector{Float64}, local_cov::Vector{Float64})
    
    # Station data
    Y = station_list.Data
    n = length(Y)
    m = m₁*m₂

    V = station_list[:,:GridCell]
    condIndIndex = [ findall(in(intersect(G.condIndSubset[j],V)),V) for j=1:length(G.condIndSubset)  ]
    condIndSubset = [ V[condIndIndex[j]] for j=1:length(condIndIndex)  ]
    S = GridPointStructure(V, condIndSubset, condIndIndex)

    V̄ = setdiff(1:m,V)
    condIndIndex = [ findall(in(intersect(G.condIndSubset[j],V̄)),V̄) for j=1:length(G.condIndSubset) ]
    condIndSubset = [ V̄[condIndIndex[j]] for j=1:length(condIndIndex) ]
    S̄ = GridPointStructure(V̄, condIndSubset, condIndIndex)

    # Covariates 
    X₁_spatial = spatial_cov
    X₁_local = local_cov
    X₁ᵢ = hcat(X₁_spatial[S.V, :], X₁_local)
    
    X₂_spatial = spatial_cov
    X₂_local = local_cov
    X₂ᵢ = hcat(X₂_spatial[S.V, :], X₂_local)

    return DataStructure(Y, X₁ᵢ, X₂ᵢ, G, S, S̄)
    
end

"""
    idf_load(stationID::AbstractString, path::AbstractString)

Loads an IDF CSV file.

"""
function idf_load(stationID::AbstractString, path::AbstractString)

    filename = path*stationID*".csv"

    df = CSV.read(filename, DataFrame)
    rename!(df,:Année => :Year)

    df2 = stack(df, Not(:Year); variable_name=:Duration, value_name=:Pcp)
    dropmissing!(df2,:Pcp)

    return df2

end

"""
    nnsearch(X::Matrix{<:Real}, points::Matrix{<:Real})

Nearest neighbor search.

"""
function nnsearch(X::Matrix{<:Real}, points::Matrix{<:Real})

    nPoints = size(points,2)
    ind = zeros(Int64,nPoints)

    for i=1:nPoints
        ind[i] = nnsearch(X,points[:,i])
    end

    return ind

end

function nnsearch(X::Matrix{<:Real}, point::Vector{<:Real})

    d = X .- point
    d² = dropdims(sum(d.^2,dims=1),dims=1)

    # Find the index of the minimum
    ind = argmin(d²)

    return ind

end

"""
    slicematrix()

"""
function slicematrix(A::AbstractMatrix{T}) where T
    
    n, m = size(A)
    
    B = Vector{T}[Vector{T}(undef, n) for _ in 1:m]
    
    for i in 1:m
        B[i] .= A[:, i]
    end
    
    return B
    
end

"""
    datalevel_loglike(Y::Vector{<:Real}, μ::Real, σ::Real, ξ::Real)

Compute the likelihood of the GEV parameters μ, σ and ξ as function of the data Y.

### Arguments
- `Y` : Vector of data.
- `μ` : GEV location parameter (real number).
- `σ` : GEV scale parameter (positive real number).
- `ξ` : GEV shape parameter (real number).

### Details

The function uses the Distributions.jl package.

### Examples

\```
 julia> datalevel_loglike(Y,0,1,.1)
\```

"""
function datalevel_loglike(Y::Vector{<:Real},μ::Real,σ::Real,ξ::Real)

    pd = GeneralizedExtremeValue(μ,σ,ξ)
    logL = loglikelihood(pd,Y)

    return logL

end

"""
    getGEV()

"""
function getGEV(C::Chains, data::DataStructure, X₁::Array{<:Real}, X₂::Array{<:Real})

    #  number of grid cells
    m = prod(data.G.gridSize)

    U = dropdims(C[:,["u[$i]" for i=1:m],1].value, dims=3)
    κ₁ = vec(C[:,"κ₁",1].value)

    V = dropdims(C[:,["v[$i]" for i=1:m],1].value, dims=3)
    κ₂ = vec(C[:,"κ₂",1].value)

    β₁ = dropdims(C[:,"β₁",1].value, dims=3)
    β₂ = dropdims(C[:,"β₂",1].value, dims=3)
    ξ = vec(C[:,"ξ",1].value)

    μ = exp.( β₁*X₁' + U )
    σ = exp.( β₂*X₂' + V )

    parmnames = vcat(["μ[$i]" for i=1:m], ["σ[$i]" for i=1:m], "ξ")
    res = vcat(μ', σ', ξ')
    
    gevChain = MambaLite.Chains(collect(res'), names=parmnames)

    return gevChain

end

function Base.write(name::AbstractString, c::MambaLite.AbstractChains)
  open(file -> serialize(file, c), name, "w")
end

function Base.read(name::AbstractString, ::Type{T}) where {T<:MambaLite.AbstractChains}
  c = open(deserialize, name, "r")
  isa(c, T) || throw(TypeError(:open, "read(\"$name\", $T)", T, c))
  c
end

function Base.write(name::AbstractString, c::Vector{Float64})
  open(file -> serialize(file, c), name, "w")
end

function Base.read(name::AbstractString, ::Type{T}) where {T<:Vector{Float64}}
  c = open(deserialize, name, "r")
  isa(c, T) || throw(TypeError(:open, "read(\"$name\", $T)", T, c))
  c
end
