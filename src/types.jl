"""
    get_types_with(parent::Type, field::Symbol)

A type creation utility meant for searching types in OMAS database. This function
returns a list of types that are fields at any level below the parent data type which
have a particular field present in it.

Example:

get_types_with(OMAS.edge_profiles, :grid_subset_index)

returns all edge_profiles types that have a subfield named grid_subset_index.
"""
function get_types_with(parent::Type, field::Symbol)
    if field ∈ fieldnames(parent)
        return [parent]
    end
    ret = Type[]
    for f ∈ fieldnames(parent)
        T = typeof(getfield(parent(), f))
        if T <: AbstractArray
            eT = eltype(T)
            if field ∈ fieldnames(eT)
                append!(ret, [eT])
            else
                append!(ret, get_types_with(eT, field))
            end
        else
            if field ∈ fieldnames(T)
                append!(ret, [T])
            else
                append!(ret, get_types_with(T, field))
            end
        end
    end
    return ret
end

edge_profiles__prop_on_subset =
    Union{get_types_with(OMAS.edge_profiles, :grid_subset_index)...}