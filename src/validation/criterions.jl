"""
    cvmcriterion()

Computes the Cramér-Von Mises criterion for one sample.
Assumes that there is no duplicate in x. 

"""
function cvmcriterion(pd::UnivariateDistribution,x̃::Array{<:Real,1})

    x = sort(x̃)
    n = length(x)

    T = 1/12/n + sum( ((2*i-1)/2/n - cdf(pd,x[i]) )^2 for i=1:n)

    ω² = T/n

    return ω²

end

"""
    cvmcriterion()

Computes the Cramér-Von Mises criterion for two samples.                              
Assumes that there is no duplicate in x and y.	

"""
function cvmcriterion(x̃::Array{<:Real,1}, ỹ::Array{<:Real,1})

    x = sort(x̃)
    y = sort!(ỹ)
    n = length(x)
    m = length(y)

    ind = [fill(1,n); fill(2,m)]

    perm = sortperm([x;y])
    ind = ind[perm]

    r = findall(ind .== 1)
    s = findall(ind .== 2)

    U = n*sum( (r[i]-i)^2 for i=1:n ) + m*sum( (s[j]-j)^2 for j=1:m )

    ω² = U/n^2/m^2 - (4*m*n-1)/6/m/n

    return ω²

end

# function get_ω̂(train::Chains, test::DataFrameRow)
    
#     train_μ = vec(dropdims(train[:,"μ[$(test.GridCell)]",1].value, dims=3));
#     train_σ = vec(dropdims(train[:,"σ[$(test.GridCell)]",1].value, dims=3));
#     train_ξ = vec(dropdims(train[:,"ξ",1].value, dims=3));
#     train_pd = GeneralizedExtremeValue.(train_μ, train_σ, train_ξ);

#     ω̂ = cvmcriterion.(train_pd, Ref(test.Data))
#     ω̂ = mean(ω̂) 
    
#     return ω̂
# end

function get_ω̂(train::Chains, test::DataFrameRow)
    
    train_μ = vec(dropdims(getindex(train, :, "μ", :).value, dims=3));
    train_σ = vec(dropdims(getindex(train, :, "σ", :).value, dims=3));
    train_ξ = vec(dropdims(getindex(train, :, "ξ", :).value, dims=3));
    
    train_pd = GeneralizedExtremeValue.(train_μ, train_σ, train_ξ);

    ω̂ = cvmcriterion.(train_pd, Ref(test.Data))
    ω̂ = mean(ω̂) 
    
    return ω̂
end

function get_ω̄(train::DataFrame, test::DataFrameRow)
    
    nearest_station = nearest_station_finder(train, test)
    nearest_station_data = nearest_station.Data
    pd_nearest_station = Extremes.getdistribution(Extremes.gumbelfitpwm(nearest_station_data))[1]

    ω̄ = cvmcriterion.(pd_nearest_station, Ref(test.Data))
    ω̄ = mean(ω̄) 
     
    return ω̄
end