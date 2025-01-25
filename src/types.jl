export get_types_with

"""
    get_types_with(parent::Type, field::Symbol)

A type creation utility meant for searching types in IMAS database. This function
returns a list of types that are fields at any level below the parent data type which
have a particular field present in it.

Example:

```julia
get_types_with(IMASdd.edge_profiles, :grid_subset_index)
```

returns all edge_profiles types that have a subfield named grid_subset_index. Note that
this function can give errors sometimes and for sake of speed and stability it is
recommended to always create union of types with hardcodes type names.
"""
function get_types_with(parent::Type, field::Symbol)
    if field ∈ fieldnames(parent)
        return [parent]
    end
    ret = Type[]
    for f ∈ fieldnames(parent)
        T = typeof(getfield(parent(), f))
        if T <: AbstractArray || T <: Tuple
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
    Union{
        IMASdd.edge_profiles__ggd___a_field_parallel{Float64},
        IMASdd.edge_profiles__ggd___e_field{Float64},
        IMASdd.edge_profiles__ggd___electrons__density{Float64},
        IMASdd.edge_profiles__ggd___electrons__density_fast{Float64},
        IMASdd.edge_profiles__ggd___electrons__distribution_function{Float64},
        IMASdd.edge_profiles__ggd___electrons__pressure{Float64},
        IMASdd.edge_profiles__ggd___electrons__pressure_fast_parallel{Float64},
        IMASdd.edge_profiles__ggd___electrons__pressure_fast_perpendicular{Float64},
        IMASdd.edge_profiles__ggd___electrons__temperature{Float64},
        IMASdd.edge_profiles__ggd___electrons__velocity{Float64},
        IMASdd.edge_profiles__ggd___ion___density{Float64},
        IMASdd.edge_profiles__ggd___ion___density_fast{Float64},
        IMASdd.edge_profiles__ggd___ion___energy_density_kinetic{Float64},
        IMASdd.edge_profiles__ggd___ion___pressure{Float64},
        IMASdd.edge_profiles__ggd___ion___pressure_fast_parallel{Float64},
        IMASdd.edge_profiles__ggd___ion___pressure_fast_perpendicular{Float64},
        IMASdd.edge_profiles__ggd___ion___state___density{Float64},
        IMASdd.edge_profiles__ggd___ion___state___density_fast{Float64},
        IMASdd.edge_profiles__ggd___ion___state___distribution_function{Float64},
        IMASdd.edge_profiles__ggd___ion___state___energy_density_kinetic{Float64},
        IMASdd.edge_profiles__ggd___ion___state___ionisation_potential{Float64},
        IMASdd.edge_profiles__ggd___ion___state___pressure{Float64},
        IMASdd.edge_profiles__ggd___ion___state___pressure_fast_parallel{Float64},
        IMASdd.edge_profiles__ggd___ion___state___pressure_fast_perpendicular{Float64},
        IMASdd.edge_profiles__ggd___ion___state___temperature{Float64},
        IMASdd.edge_profiles__ggd___ion___state___velocity{Float64},
        IMASdd.edge_profiles__ggd___ion___state___velocity_diamagnetic{Float64},
        IMASdd.edge_profiles__ggd___ion___state___velocity_exb{Float64},
        IMASdd.edge_profiles__ggd___ion___state___z_average{Float64},
        IMASdd.edge_profiles__ggd___ion___state___z_square_average{Float64},
        IMASdd.edge_profiles__ggd___ion___temperature{Float64},
        IMASdd.edge_profiles__ggd___ion___velocity{Float64},
        IMASdd.edge_profiles__ggd___j_anomalous{Float64},
        IMASdd.edge_profiles__ggd___j_diamagnetic{Float64},
        IMASdd.edge_profiles__ggd___j_heat_viscosity{Float64},
        IMASdd.edge_profiles__ggd___j_inertial{Float64},
        IMASdd.edge_profiles__ggd___j_ion_neutral_friction{Float64},
        IMASdd.edge_profiles__ggd___j_parallel{Float64},
        IMASdd.edge_profiles__ggd___j_parallel_viscosity{Float64},
        IMASdd.edge_profiles__ggd___j_perpendicular_viscosity{Float64},
        IMASdd.edge_profiles__ggd___j_pfirsch_schlueter{Float64},
        IMASdd.edge_profiles__ggd___j_total{Float64},
        IMASdd.edge_profiles__ggd___n_i_total_over_n_e{Float64},
        IMASdd.edge_profiles__ggd___neutral___density{Float64},
        IMASdd.edge_profiles__ggd___neutral___density_fast{Float64},
        IMASdd.edge_profiles__ggd___neutral___energy_density_kinetic{Float64},
        IMASdd.edge_profiles__ggd___neutral___pressure{Float64},
        IMASdd.edge_profiles__ggd___neutral___pressure_fast_parallel{Float64},
        IMASdd.edge_profiles__ggd___neutral___pressure_fast_perpendicular{Float64},
        IMASdd.edge_profiles__ggd___neutral___state___density{Float64},
        IMASdd.edge_profiles__ggd___neutral___state___density_fast{Float64},
        IMASdd.edge_profiles__ggd___neutral___state___distribution_function{Float64},
        IMASdd.edge_profiles__ggd___neutral___state___energy_density_kinetic{Float64},
        IMASdd.edge_profiles__ggd___neutral___state___pressure{Float64},
        IMASdd.edge_profiles__ggd___neutral___state___pressure_fast_parallel{Float64},
        IMASdd.edge_profiles__ggd___neutral___state___pressure_fast_perpendicular{
            Float64,
        },
        IMASdd.edge_profiles__ggd___neutral___state___temperature{Float64},
        IMASdd.edge_profiles__ggd___neutral___state___velocity{Float64},
        IMASdd.edge_profiles__ggd___neutral___state___velocity_diamagnetic{Float64},
        IMASdd.edge_profiles__ggd___neutral___state___velocity_exb{Float64},
        IMASdd.edge_profiles__ggd___neutral___temperature{Float64},
        IMASdd.edge_profiles__ggd___neutral___velocity{Float64},
        IMASdd.edge_profiles__ggd___phi_potential{Float64},
        IMASdd.edge_profiles__ggd___pressure_parallel{Float64},
        IMASdd.edge_profiles__ggd___pressure_perpendicular{Float64},
        IMASdd.edge_profiles__ggd___pressure_thermal{Float64},
        IMASdd.edge_profiles__ggd___t_i_average{Float64},
        IMASdd.edge_profiles__ggd___zeff{Float64},
        IMASdd.edge_profiles__ggd_fast___electrons__density{Float64},
        IMASdd.edge_profiles__ggd_fast___electrons__temperature{Float64},
        IMASdd.edge_profiles__ggd_fast___energy_thermal{Float64},
        IMASdd.edge_profiles__ggd_fast___ion___content{Float64},
        IMASdd.edge_profiles__ggd_fast___ion___density{Float64},
        IMASdd.edge_profiles__ggd_fast___ion___temperature{Float64},
        IMASdd.edge_profiles__statistics___quantity_2d___statistics_type{Float64},
    }
