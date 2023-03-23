function nearest_station_finder(train::DataFrame, test::DataFrameRow)
    
    stationLocation = collect([train[:,:Lat] train[:,:Lon]]')
    
    nearest_station_idx = nnsearch(stationLocation, [test.Lat, test.Lon])
    
    return train[nearest_station_idx, :]
end


function haversine_dist(p₁::Vector{Float64}, p₂::Vector{Float64})

    @assert length(p₁) == length(p₂) == 2

    earth_radius = 6371

    # Lat
    φ₁ = deg2rad(p₁[1])
    φ₂ = deg2rad(p₂[1])
    Δφ = φ₂ - φ₁

    # Lon
    Δλ = deg2rad(p₂[2] - p₁[2])

    # Haversine formula
    a = sin(Δφ/2)^2 + cos(φ₁)*cos(φ₂)*sin(Δλ/2)^2
    c = 2*atan(sqrt(a), sqrt(1-a))

    return c * earth_radius

end
