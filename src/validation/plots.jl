function qqnaive(train::DataFrame, test::DataFrameRow)
    
    # Sample quantiles
    y_test, p_test = Extremes.ecdf(test.Data)

    # Predicted quantiles
    nearest_station = nearest_station_finder(train, test)
    nearest_station_data = nearest_station.Data
    pd_nearest_station = Extremes.getdistribution(Extremes.gumbelfitpwm(nearest_station_data))[1]
    q_nearest_station = quantile.(pd_nearest_station, p_test)

    qq_data = DataFrame(Model = q_nearest_station, Empirical = y_test);

#     plotTitle = "Naïve approach"

#     set_default_plot_size(8cm, 8cm)
#     p1 = Gadfly.plot(qq_data, y=:Model, x=:Empirical, Geom.point, Geom.abline(color="red", style=:dash),
#                 Guide.xlabel("Sample quantiles"), Guide.ylabel("Predicted quantiles"), Guide.title(plotTitle),
#                 Theme(discrete_highlight_color=c->nothing),
#                 Coord.cartesian(xmin=0, xmax=120, ymin=0, ymax=120));
    
#     return p1
    
    return qq_data
end

function qqspatial(train::Chains, test::DataFrameRow)
    
    # Sample quantiles
    y_test, p_test = Extremes.ecdf(test.Data)

    # Predicted quantiles
    train_μ = vec(dropdims(getindex(train, :, "μ", :).value, dims=3));
    train_σ = vec(dropdims(getindex(train, :, "σ", :).value, dims=3));
    train_ξ = vec(dropdims(getindex(train, :, "ξ", :).value, dims=3));
    train_pd = GeneralizedExtremeValue.(train_μ, train_σ, train_ξ);
    q_sim = Array{Float64}(undef, length(train_pd), length(y_test))
    i = 1
    for d in eachslice(train_pd, dims=1)
        q_sim[i,:] = quantile.(d, p_test)
        i +=1
    end
    q_train = vec(mean(q_sim, dims=1))

    qq_data = DataFrame(Model = q_train, Empirical = y_test);

#     plotTitle = "Spatial approach"

#     set_default_plot_size(8cm, 8cm)
#     p2 = Gadfly.plot(qq_data, y=:Model, x=:Empirical, Geom.point, Geom.abline(color="red", style=:dash),
#                 Guide.xlabel("Sample quantiles"), Guide.ylabel("Predicted quantiles"), Guide.title(plotTitle),
#                 Theme(discrete_highlight_color=c->nothing),
#                 Coord.cartesian(xmin=0, xmax=120, ymin=0, ymax=120));
    
#     return p2
    
    return qq_data
end

# function qqspatial(train::Chains, test::DataFrameRow)
    
#     # Sample quantiles
#     y_test, p_test = Extremes.ecdf(test.Data)

#     # Predicted quantiles
#     train_μ = vec(dropdims(train[:,"μ[$(test.GridCell)]",1].value, dims=3));
#     train_σ = vec(dropdims(train[:,"σ[$(test.GridCell)]",1].value, dims=3));
#     train_ξ = vec(dropdims(train[:,"ξ",1].value, dims=3));
#     train_pd = GeneralizedExtremeValue.(train_μ, train_σ, train_ξ);
#     q_sim = Array{Float64}(undef, length(train_pd), length(y_test))
#     i = 1
#     for d in eachslice(train_pd, dims=1)
#         q_sim[i,:] = quantile.(d, p_test)
#         i +=1
#     end
#     q_train = vec(mean(q_sim, dims=1))

#     qq_data = DataFrame(Model = q_train, Empirical = y_test);

# #     plotTitle = "Spatial approach"

# #     set_default_plot_size(8cm, 8cm)
# #     p2 = Gadfly.plot(qq_data, y=:Model, x=:Empirical, Geom.point, Geom.abline(color="red", style=:dash),
# #                 Guide.xlabel("Sample quantiles"), Guide.ylabel("Predicted quantiles"), Guide.title(plotTitle),
# #                 Theme(discrete_highlight_color=c->nothing),
# #                 Coord.cartesian(xmin=0, xmax=120, ymin=0, ymax=120));
    
# #     return p2
    
#     return qq_data
# end