# struct DataStructure
#     Y::Vector{Vector{Float64}}
#     X₁ᵢ::Vector{Float64}
#     X₂ᵢ::Vector{Float64}
#     G::GMRF.GridStructure
#     S::GridPointStructure
#     S̄::GridPointStructure
# end

struct DataStructure
    Y::Vector{Vector{Float64}}
    X₁ᵢ::Array{Float64, 2}
    X₂ᵢ::Array{Float64, 2}
    G::GMRF.GridStructure
    S::GridPointStructure
    S̄::GridPointStructure
end

# function showDataStructure(io::IO, obj::DataStructure; prefix::String = "")
#     println(io, prefix, "DataStructure")
#     println(io, prefix, "Y :\t\t", typeof(obj.Y), "[", length(obj.Y), "]")
#     println(io, prefix, "X₁ᵢ :\t\t", typeof(obj.X₁ᵢ), "[", length(obj.X₁ᵢ), "]")
#     println(io, prefix, "X₂ᵢ :\t\t", typeof(obj.X₂ᵢ), "[", length(obj.X₂ᵢ), "]")
#     println(io, prefix, "G :")
#     GMRF.showGridStructure(io, obj.G, prefix = prefix*"\t")
#     println(io)
#     println(io, prefix, "S :")
#     showGridPointStructure(io, obj.S, prefix = prefix*"\t")
#     println(io)
#     println(io, prefix, "S̄ :")
#     showGridPointStructure(io, obj.S̄, prefix = prefix*"\t")
# end

function showDataStructure(io::IO, obj::DataStructure; prefix::String = "")
    println(io, prefix, "DataStructure")
    println(io, prefix, "Y :\t\t", typeof(obj.Y), "[", length(obj.Y), "]")
    println(io, prefix, "X₁ᵢ :\t\t", typeof(obj.X₁ᵢ), "[", size(obj.X₁ᵢ), "]")
    println(io, prefix, "X₂ᵢ :\t\t", typeof(obj.X₂ᵢ), "[", size(obj.X₂ᵢ), "]")
    println(io, prefix, "G :")
    GMRF.showGridStructure(io, obj.G, prefix = prefix*"\t")
    println(io)
    println(io, prefix, "S :")
    showGridPointStructure(io, obj.S, prefix = prefix*"\t")
    println(io)
    println(io, prefix, "S̄ :")
    showGridPointStructure(io, obj.S̄, prefix = prefix*"\t")
end

function Base.show(io::IO, obj::DataStructure)

    showDataStructure(io, obj)

end