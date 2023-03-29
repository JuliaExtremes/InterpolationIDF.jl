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

function densityplot(c::MambaLite.AbstractChains; legend::Bool=false,
                     trim::Tuple{Real, Real}=(0.025, 0.975), na...)
  nrows, nvars, nchains = size(c.value)
  plots = Array{Plot}(undef, nvars)
  pos = legend ? :right : :none
  for i in 1:nvars
    val = Array{Vector{Float64}}(undef, nchains)
    for j in 1:nchains
      qs = quantile(c.value[:, i, j], [trim[1], trim[2]])
      val[j] = c.value[qs[1] .<= c.value[:, i, j] .<= qs[2], i, j]
    end
    plots[i] = Gadfly.plot(x=[val...;], Geom.density(),
                    color=repeat(c.chains, inner=[length(c.range)]),
                    Scale.color_discrete(), Guide.colorkey(title="Chain"),
                    Guide.xlabel("Value", orientation=:horizontal),
                    Guide.ylabel("Density", orientation=:vertical),
                    Guide.title(c.names[i]), Theme(key_position=pos))
  end
  return plots
end

function trace_plot(c::MambaLite.AbstractChains)
    fig1 = Gadfly.plot(y = c.value, Geom.line,
    Guide.xlabel("Iteration"),
    Guide.ylabel("Value"),
    Guide.title(c.names[1]))

    fig2 = densityplot(c)[1]

    Gadfly.set_default_plot_size(24cm, 8cm)
    return hstack(fig1,fig2)
end
