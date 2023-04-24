"""
	GridPointStructure(V::Vector{Int64},
		CondIndSubset::Vector{Vector{Int64}},
		CondIndIndex::Vector{Vector{Int64}})

Creates a GridPointStructure type where
* V 
* CondIndSubset
* CondIndIndex

"""
struct GridPointStructure
    V::Vector{Int64}
    CondIndSubset::Vector{Vector{Int64}}
    CondIndIndex::Vector{Vector{Int64}}
end

"""
    showGridPointStructure(io::IO, obj::GridPointStructure; prefix::String = "")

Displays a GridPointStructure with the prefix `prefix` before every line.

"""
function showGridPointStructure(io::IO, obj::GridPointStructure; prefix::String = "")
    println(io, prefix, "GridPointStructure")
    println(io, prefix, "V :\t\t", typeof(obj.V), "[", length(obj.V), "]")
    println(io, prefix, "CondIndSubset :\t\t", typeof(obj.CondIndSubset), "[", length(obj.CondIndSubset), "]")
    println(io, prefix, "CondIndIndex :\t\t", typeof(obj.CondIndIndex), "[", length(obj.CondIndIndex), "]")
end

function Base.show(io::IO, obj::GridPointStructure)

    showGridPointStructure(io, obj)

end
