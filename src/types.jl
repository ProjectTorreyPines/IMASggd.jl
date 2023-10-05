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