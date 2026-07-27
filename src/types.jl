export get_types_with

"""
    get_types_with(
        parent::Type,
        field::Symbol;
        return_field_type::Bool=false,
    )::Array{Type}

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
function get_types_with(
    parent::Type,
    field::Symbol;
    return_field_type::Bool=false,
)::Array{Type}
    if field ∈ fieldnames(parent)
        ret = [parent]
    else
        ret = Type[]
        for f ∈ fieldnames(parent)
            if !startswith(string(f), "_")
                T = fieldtype(parent, f)
                if !(T <: AbstractArray{<:Number} || T <: Number)
                    if T <: AbstractArray || T <: Tuple
                        eT = eltype(T)
                        if eT <: Number
                            continue
                        end
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
            end
        end
    end
    if return_field_type
        fret = Type[]
        for T ∈ ret
            fT = fieldtype(T, field)
            if fT <: AbstractArray || fT <: Tuple
                append!(fret, [eltype(fT)])
            else
                append!(fret, [fT])
            end
        end
        return fret
    else
        return ret
    end
end

"""
    all__grid_ggd =
        Union{
            IMASdd.edge_profiles__grid_ggd{T},
            IMASdd.edge_sources__grid_ggd{T},
            IMASdd.edge_transport__grid_ggd{T},
            IMASdd.em_coupling__grid_ggd{T},
            IMASdd.ferritic__grid_ggd{T},
            IMASdd.mhd__grid_ggd{T},
            IMASdd.radiation__grid_ggd{T},
            IMASdd.runaway_electrons__grid_ggd{T},
            IMASdd.wall__description_ggd___grid_ggd{T},
        } where {T};

Union of all IMAS data dictionary `grid_ggd` types.
"""
all__grid_ggd =
    Union{
        IMASdd.edge_profiles__grid_ggd{T},
        IMASdd.edge_sources__grid_ggd{T},
        IMASdd.edge_transport__grid_ggd{T},
        IMASdd.em_coupling__grid_ggd{T},
        IMASdd.ferritic__grid_ggd{T},
        IMASdd.mhd__grid_ggd{T},
        IMASdd.radiation__grid_ggd{T},
        IMASdd.runaway_electrons__grid_ggd{T},
        IMASdd.wall__description_ggd___grid_ggd{T},
    } where {T};

"""
    all__space =
        Union{
            IMASdd.distribution_sources__source___ggd___grid__space{T},
            IMASdd.distributions__distribution___ggd___grid__space{T},
            IMASdd.edge_profiles__grid_ggd___space{T},
            IMASdd.edge_sources__grid_ggd___space{T},
            IMASdd.edge_transport__grid_ggd___space{T},
            IMASdd.em_coupling__grid_ggd___space{T},
            IMASdd.equilibrium__grids_ggd___grid___space{T},
            IMASdd.ferritic__grid_ggd__space{T},
            IMASdd.mhd__grid_ggd___space{T},
            IMASdd.radiation__grid_ggd___space{T},
            IMASdd.runaway_electrons__grid_ggd___space{T},
            IMASdd.tf__field_map___grid__space{T},
            IMASdd.transport_solver_numerics__boundary_conditions_ggd___grid__space{T},
            IMASdd.wall__description_ggd___grid_ggd___space{T},
            IMASdd.waves__coherent_wave___full_wave___grid__space{T},
        } where {T};

Union of all `space` types that are attributes of `grid_ggd` objects in IMAS.
"""
all__space =
    Union{
        IMASdd.distribution_sources__source___ggd___grid__space{T},
        IMASdd.distributions__distribution___ggd___grid__space{T},
        IMASdd.edge_profiles__grid_ggd___space{T},
        IMASdd.edge_sources__grid_ggd___space{T},
        IMASdd.edge_transport__grid_ggd___space{T},
        IMASdd.em_coupling__grid_ggd___space{T},
        IMASdd.equilibrium__grids_ggd___grid___space{T},
        IMASdd.ferritic__grid_ggd__space{T},
        IMASdd.mhd__grid_ggd___space{T},
        IMASdd.radiation__grid_ggd___space{T},
        IMASdd.runaway_electrons__grid_ggd___space{T},
        IMASdd.tf__field_map___grid__space{T},
        IMASdd.transport_solver_numerics__boundary_conditions_ggd___grid__space{T},
        IMASdd.wall__description_ggd___grid_ggd___space{T},
        IMASdd.waves__coherent_wave___full_wave___grid__space{T},
    } where {T};

"""
    all__grid_subset =
        Union{
            IMASdd.edge_profiles__grid_ggd___grid_subset{T},
            IMASdd.edge_sources__grid_ggd___grid_subset{T},
            IMASdd.edge_transport__grid_ggd___grid_subset{T},
            IMASdd.em_coupling__grid_ggd___grid_subset{T},
            IMASdd.ferritic__grid_ggd__grid_subset{T},
            IMASdd.mhd__grid_ggd___grid_subset{T},
            IMASdd.radiation__grid_ggd___grid_subset{T},
            IMASdd.runaway_electrons__grid_ggd___grid_subset{T},
            IMASdd.wall__description_ggd___grid_ggd___grid_subset{T},
        } where {T};

Union of all `grid_subset` types are attributes of `grid_ggd` objects in IMAS.
"""
all__grid_subset =
    Union{
        IMASdd.edge_profiles__grid_ggd___grid_subset{T},
        IMASdd.edge_sources__grid_ggd___grid_subset{T},
        IMASdd.edge_transport__grid_ggd___grid_subset{T},
        IMASdd.em_coupling__grid_ggd___grid_subset{T},
        IMASdd.ferritic__grid_ggd__grid_subset{T},
        IMASdd.mhd__grid_ggd___grid_subset{T},
        IMASdd.radiation__grid_ggd___grid_subset{T},
        IMASdd.runaway_electrons__grid_ggd___grid_subset{T},
        IMASdd.wall__description_ggd___grid_ggd___grid_subset{T},
    } where {T};

"""
    all__ggd =
        Union{
            IMASdd.distribution_sources__source___ggd{T},
            IMASdd.distributions__distribution___ggd{T},
            IMASdd.edge_profiles__ggd{T},
            IMASdd.edge_sources__source___ggd{T},
            IMASdd.edge_transport__model___ggd{T},
            IMASdd.equilibrium__time_slice___ggd{T},
            IMASdd.mhd__ggd{T},
            IMASdd.radiation__process___ggd{T},
            IMASdd.runaway_electrons__distribution__ggd{T},
            IMASdd.wall__description_ggd___ggd{T},
        } where {T};

Union of all Generalized Grid Description (ggd) objects in IMAS.
"""
all__ggd =
    Union{
        IMASdd.distribution_sources__source___ggd{T},
        IMASdd.distributions__distribution___ggd{T},
        IMASdd.edge_profiles__ggd{T},
        IMASdd.edge_sources__source___ggd{T},
        IMASdd.edge_transport__model___ggd{T},
        IMASdd.equilibrium__time_slice___ggd{T},
        IMASdd.mhd__ggd{T},
        IMASdd.radiation__process___ggd{T},
        IMASdd.runaway_electrons__distribution__ggd{T},
        IMASdd.wall__description_ggd___ggd{T},
    } where {T};

"""
    all__grid_subset_prop

A large union of all `ggd` properties that refer to a `grid_subset` for the dimensions
of the data.
"""
all__grid_subset_prop =
    Union{
        IMASdd.core_profiles__statistics___quantity_2d___statistics_type{T},
        IMASdd.distribution_sources__source___ggd___particles{T},
        IMASdd.distributions__distribution___ggd___expansion___grid_subset{T},
        IMASdd.distributions__distribution___ggd___expansion_fd3v___grid_subset{T},
        IMASdd.edge_profiles__ggd___a_field_parallel{T},
        IMASdd.edge_profiles__ggd___e_field{T},
        IMASdd.edge_profiles__ggd___electrons__density{T},
        IMASdd.edge_profiles__ggd___electrons__density_fast{T},
        IMASdd.edge_profiles__ggd___electrons__distribution_function{T},
        IMASdd.edge_profiles__ggd___electrons__pressure{T},
        IMASdd.edge_profiles__ggd___electrons__pressure_fast_parallel{T},
        IMASdd.edge_profiles__ggd___electrons__pressure_fast_perpendicular{T},
        IMASdd.edge_profiles__ggd___electrons__temperature{T},
        IMASdd.edge_profiles__ggd___electrons__velocity{T},
        IMASdd.edge_profiles__ggd___ion___density{T},
        IMASdd.edge_profiles__ggd___ion___density_fast{T},
        IMASdd.edge_profiles__ggd___ion___energy_density_kinetic{T},
        IMASdd.edge_profiles__ggd___ion___pressure{T},
        IMASdd.edge_profiles__ggd___ion___pressure_fast_parallel{T},
        IMASdd.edge_profiles__ggd___ion___pressure_fast_perpendicular{T},
        IMASdd.edge_profiles__ggd___ion___state___density{T},
        IMASdd.edge_profiles__ggd___ion___state___density_fast{T},
        IMASdd.edge_profiles__ggd___ion___state___distribution_function{T},
        IMASdd.edge_profiles__ggd___ion___state___energy_density_kinetic{T},
        IMASdd.edge_profiles__ggd___ion___state___ionisation_potential{T},
        IMASdd.edge_profiles__ggd___ion___state___pressure{T},
        IMASdd.edge_profiles__ggd___ion___state___pressure_fast_parallel{T},
        IMASdd.edge_profiles__ggd___ion___state___pressure_fast_perpendicular{T},
        IMASdd.edge_profiles__ggd___ion___state___temperature{T},
        IMASdd.edge_profiles__ggd___ion___state___velocity{T},
        IMASdd.edge_profiles__ggd___ion___state___velocity_diamagnetic{T},
        IMASdd.edge_profiles__ggd___ion___state___velocity_exb{T},
        IMASdd.edge_profiles__ggd___ion___state___z_average{T},
        IMASdd.edge_profiles__ggd___ion___state___z_square_average{T},
        IMASdd.edge_profiles__ggd___ion___temperature{T},
        IMASdd.edge_profiles__ggd___ion___velocity{T},
        IMASdd.edge_profiles__ggd___j_anomalous{T},
        IMASdd.edge_profiles__ggd___j_diamagnetic{T},
        IMASdd.edge_profiles__ggd___j_heat_viscosity{T},
        IMASdd.edge_profiles__ggd___j_inertial{T},
        IMASdd.edge_profiles__ggd___j_ion_neutral_friction{T},
        IMASdd.edge_profiles__ggd___j_parallel{T},
        IMASdd.edge_profiles__ggd___j_parallel_viscosity{T},
        IMASdd.edge_profiles__ggd___j_perpendicular_viscosity{T},
        IMASdd.edge_profiles__ggd___j_pfirsch_schlueter{T},
        IMASdd.edge_profiles__ggd___j_total{T},
        IMASdd.edge_profiles__ggd___n_i_total_over_n_e{T},
        IMASdd.edge_profiles__ggd___neutral___density{T},
        IMASdd.edge_profiles__ggd___neutral___density_fast{T},
        IMASdd.edge_profiles__ggd___neutral___energy_density_kinetic{T},
        IMASdd.edge_profiles__ggd___neutral___pressure{T},
        IMASdd.edge_profiles__ggd___neutral___pressure_fast_parallel{T},
        IMASdd.edge_profiles__ggd___neutral___pressure_fast_perpendicular{T},
        IMASdd.edge_profiles__ggd___neutral___state___density{T},
        IMASdd.edge_profiles__ggd___neutral___state___density_fast{T},
        IMASdd.edge_profiles__ggd___neutral___state___distribution_function{T},
        IMASdd.edge_profiles__ggd___neutral___state___energy_density_kinetic{T},
        IMASdd.edge_profiles__ggd___neutral___state___pressure{T},
        IMASdd.edge_profiles__ggd___neutral___state___pressure_fast_parallel{T},
        IMASdd.edge_profiles__ggd___neutral___state___pressure_fast_perpendicular{T},
        IMASdd.edge_profiles__ggd___neutral___state___temperature{T},
        IMASdd.edge_profiles__ggd___neutral___state___velocity{T},
        IMASdd.edge_profiles__ggd___neutral___state___velocity_diamagnetic{T},
        IMASdd.edge_profiles__ggd___neutral___state___velocity_exb{T},
        IMASdd.edge_profiles__ggd___neutral___temperature{T},
        IMASdd.edge_profiles__ggd___neutral___velocity{T},
        IMASdd.edge_profiles__ggd___phi_potential{T},
        IMASdd.edge_profiles__ggd___pressure_parallel{T},
        IMASdd.edge_profiles__ggd___pressure_perpendicular{T},
        IMASdd.edge_profiles__ggd___pressure_thermal{T},
        IMASdd.edge_profiles__ggd___t_i_average{T},
        IMASdd.edge_profiles__ggd___zeff{T},
        IMASdd.edge_profiles__ggd_fast___electrons__density{T},
        IMASdd.edge_profiles__ggd_fast___electrons__temperature{T},
        IMASdd.edge_profiles__ggd_fast___energy_thermal{T},
        IMASdd.edge_profiles__ggd_fast___ion___content{T},
        IMASdd.edge_profiles__ggd_fast___ion___density{T},
        IMASdd.edge_profiles__ggd_fast___ion___temperature{T},
        IMASdd.edge_profiles__statistics___quantity_2d___statistics_type{T},
        IMASdd.edge_sources__source___ggd___current{T},
        IMASdd.edge_sources__source___ggd___electrons__energy{T},
        IMASdd.edge_sources__source___ggd___electrons__particles{T},
        IMASdd.edge_sources__source___ggd___ion___energy{T},
        IMASdd.edge_sources__source___ggd___ion___momentum{T},
        IMASdd.edge_sources__source___ggd___ion___particles{T},
        IMASdd.edge_sources__source___ggd___ion___state___energy{T},
        IMASdd.edge_sources__source___ggd___ion___state___momentum{T},
        IMASdd.edge_sources__source___ggd___ion___state___particles{T},
        IMASdd.edge_sources__source___ggd___momentum{T},
        IMASdd.edge_sources__source___ggd___neutral___energy{T},
        IMASdd.edge_sources__source___ggd___neutral___momentum{T},
        IMASdd.edge_sources__source___ggd___neutral___particles{T},
        IMASdd.edge_sources__source___ggd___neutral___state___energy{T},
        IMASdd.edge_sources__source___ggd___neutral___state___momentum{T},
        IMASdd.edge_sources__source___ggd___neutral___state___particles{T},
        IMASdd.edge_sources__source___ggd___total_ion_energy{T},
        IMASdd.edge_sources__source___ggd_fast___ion___power{T},
        IMASdd.edge_transport__model___ggd___conductivity{T},
        IMASdd.edge_transport__model___ggd___electrons__energy__d{T},
        IMASdd.edge_transport__model___ggd___electrons__energy__d_pol{T},
        IMASdd.edge_transport__model___ggd___electrons__energy__d_radial{T},
        IMASdd.edge_transport__model___ggd___electrons__energy__flux{T},
        IMASdd.edge_transport__model___ggd___electrons__energy__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___electrons__energy__flux_pol{T},
        IMASdd.edge_transport__model___ggd___electrons__energy__flux_radial{T},
        IMASdd.edge_transport__model___ggd___electrons__energy__v{T},
        IMASdd.edge_transport__model___ggd___electrons__energy__v_pol{T},
        IMASdd.edge_transport__model___ggd___electrons__energy__v_radial{T},
        IMASdd.edge_transport__model___ggd___electrons__particles__d{T},
        IMASdd.edge_transport__model___ggd___electrons__particles__d_pol{T},
        IMASdd.edge_transport__model___ggd___electrons__particles__d_radial{T},
        IMASdd.edge_transport__model___ggd___electrons__particles__flux{T},
        IMASdd.edge_transport__model___ggd___electrons__particles__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___electrons__particles__flux_pol{T},
        IMASdd.edge_transport__model___ggd___electrons__particles__flux_radial{T},
        IMASdd.edge_transport__model___ggd___electrons__particles__v{T},
        IMASdd.edge_transport__model___ggd___electrons__particles__v_pol{T},
        IMASdd.edge_transport__model___ggd___electrons__particles__v_radial{T},
        IMASdd.edge_transport__model___ggd___ion___energy__d{T},
        IMASdd.edge_transport__model___ggd___ion___energy__d_pol{T},
        IMASdd.edge_transport__model___ggd___ion___energy__d_radial{T},
        IMASdd.edge_transport__model___ggd___ion___energy__flux{T},
        IMASdd.edge_transport__model___ggd___ion___energy__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___ion___energy__flux_pol{T},
        IMASdd.edge_transport__model___ggd___ion___energy__flux_radial{T},
        IMASdd.edge_transport__model___ggd___ion___energy__v{T},
        IMASdd.edge_transport__model___ggd___ion___energy__v_pol{T},
        IMASdd.edge_transport__model___ggd___ion___energy__v_radial{T},
        IMASdd.edge_transport__model___ggd___ion___momentum__d{T},
        IMASdd.edge_transport__model___ggd___ion___momentum__d_pol{T},
        IMASdd.edge_transport__model___ggd___ion___momentum__d_radial{T},
        IMASdd.edge_transport__model___ggd___ion___momentum__flux{T},
        IMASdd.edge_transport__model___ggd___ion___momentum__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___ion___momentum__flux_pol{T},
        IMASdd.edge_transport__model___ggd___ion___momentum__flux_radial{T},
        IMASdd.edge_transport__model___ggd___ion___momentum__v{T},
        IMASdd.edge_transport__model___ggd___ion___momentum__v_pol{T},
        IMASdd.edge_transport__model___ggd___ion___momentum__v_radial{T},
        IMASdd.edge_transport__model___ggd___ion___particles__d{T},
        IMASdd.edge_transport__model___ggd___ion___particles__d_pol{T},
        IMASdd.edge_transport__model___ggd___ion___particles__d_radial{T},
        IMASdd.edge_transport__model___ggd___ion___particles__flux{T},
        IMASdd.edge_transport__model___ggd___ion___particles__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___ion___particles__flux_pol{T},
        IMASdd.edge_transport__model___ggd___ion___particles__flux_radial{T},
        IMASdd.edge_transport__model___ggd___ion___particles__v{T},
        IMASdd.edge_transport__model___ggd___ion___particles__v_pol{T},
        IMASdd.edge_transport__model___ggd___ion___particles__v_radial{T},
        IMASdd.edge_transport__model___ggd___ion___state___energy__d{T},
        IMASdd.edge_transport__model___ggd___ion___state___energy__d_pol{T},
        IMASdd.edge_transport__model___ggd___ion___state___energy__d_radial{T},
        IMASdd.edge_transport__model___ggd___ion___state___energy__flux{T},
        IMASdd.edge_transport__model___ggd___ion___state___energy__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___ion___state___energy__flux_pol{T},
        IMASdd.edge_transport__model___ggd___ion___state___energy__flux_radial{T},
        IMASdd.edge_transport__model___ggd___ion___state___energy__v{T},
        IMASdd.edge_transport__model___ggd___ion___state___energy__v_pol{T},
        IMASdd.edge_transport__model___ggd___ion___state___energy__v_radial{T},
        IMASdd.edge_transport__model___ggd___ion___state___momentum__d{T},
        IMASdd.edge_transport__model___ggd___ion___state___momentum__d_pol{T},
        IMASdd.edge_transport__model___ggd___ion___state___momentum__d_radial{T},
        IMASdd.edge_transport__model___ggd___ion___state___momentum__flux{T},
        IMASdd.edge_transport__model___ggd___ion___state___momentum__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___ion___state___momentum__flux_pol{T},
        IMASdd.edge_transport__model___ggd___ion___state___momentum__flux_radial{T},
        IMASdd.edge_transport__model___ggd___ion___state___momentum__v{T},
        IMASdd.edge_transport__model___ggd___ion___state___momentum__v_pol{T},
        IMASdd.edge_transport__model___ggd___ion___state___momentum__v_radial{T},
        IMASdd.edge_transport__model___ggd___ion___state___particles__d{T},
        IMASdd.edge_transport__model___ggd___ion___state___particles__d_pol{T},
        IMASdd.edge_transport__model___ggd___ion___state___particles__d_radial{T},
        IMASdd.edge_transport__model___ggd___ion___state___particles__flux{T},
        IMASdd.edge_transport__model___ggd___ion___state___particles__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___ion___state___particles__flux_pol{T},
        IMASdd.edge_transport__model___ggd___ion___state___particles__flux_radial{T},
        IMASdd.edge_transport__model___ggd___ion___state___particles__v{T},
        IMASdd.edge_transport__model___ggd___ion___state___particles__v_pol{T},
        IMASdd.edge_transport__model___ggd___ion___state___particles__v_radial{T},
        IMASdd.edge_transport__model___ggd___momentum__d{T},
        IMASdd.edge_transport__model___ggd___momentum__d_pol{T},
        IMASdd.edge_transport__model___ggd___momentum__d_radial{T},
        IMASdd.edge_transport__model___ggd___momentum__flux{T},
        IMASdd.edge_transport__model___ggd___momentum__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___momentum__flux_pol{T},
        IMASdd.edge_transport__model___ggd___momentum__flux_radial{T},
        IMASdd.edge_transport__model___ggd___momentum__v{T},
        IMASdd.edge_transport__model___ggd___momentum__v_pol{T},
        IMASdd.edge_transport__model___ggd___momentum__v_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___energy__d{T},
        IMASdd.edge_transport__model___ggd___neutral___energy__d_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___energy__d_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___energy__flux{T},
        IMASdd.edge_transport__model___ggd___neutral___energy__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___neutral___energy__flux_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___energy__flux_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___energy__v{T},
        IMASdd.edge_transport__model___ggd___neutral___energy__v_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___energy__v_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___momentum__d{T},
        IMASdd.edge_transport__model___ggd___neutral___momentum__d_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___momentum__d_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___momentum__flux{T},
        IMASdd.edge_transport__model___ggd___neutral___momentum__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___neutral___momentum__flux_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___momentum__flux_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___momentum__v{T},
        IMASdd.edge_transport__model___ggd___neutral___momentum__v_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___momentum__v_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___particles__d{T},
        IMASdd.edge_transport__model___ggd___neutral___particles__d_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___particles__d_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___particles__flux{T},
        IMASdd.edge_transport__model___ggd___neutral___particles__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___neutral___particles__flux_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___particles__flux_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___particles__v{T},
        IMASdd.edge_transport__model___ggd___neutral___particles__v_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___particles__v_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___state___energy__d{T},
        IMASdd.edge_transport__model___ggd___neutral___state___energy__d_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___state___energy__d_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___state___energy__flux{T},
        IMASdd.edge_transport__model___ggd___neutral___state___energy__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___neutral___state___energy__flux_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___state___energy__flux_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___state___energy__v{T},
        IMASdd.edge_transport__model___ggd___neutral___state___energy__v_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___state___energy__v_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___state___momentum__d{T},
        IMASdd.edge_transport__model___ggd___neutral___state___momentum__d_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___state___momentum__d_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___state___momentum__flux{T},
        IMASdd.edge_transport__model___ggd___neutral___state___momentum__flux_limiter{
            T,
        },
        IMASdd.edge_transport__model___ggd___neutral___state___momentum__flux_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___state___momentum__flux_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___state___momentum__v{T},
        IMASdd.edge_transport__model___ggd___neutral___state___momentum__v_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___state___momentum__v_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___state___particles__d{T},
        IMASdd.edge_transport__model___ggd___neutral___state___particles__d_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___state___particles__d_radial{T},
        IMASdd.edge_transport__model___ggd___neutral___state___particles__flux{T},
        IMASdd.edge_transport__model___ggd___neutral___state___particles__flux_limiter{
            T,
        },
        IMASdd.edge_transport__model___ggd___neutral___state___particles__flux_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___state___particles__flux_radial{
            T,
        },
        IMASdd.edge_transport__model___ggd___neutral___state___particles__v{T},
        IMASdd.edge_transport__model___ggd___neutral___state___particles__v_pol{T},
        IMASdd.edge_transport__model___ggd___neutral___state___particles__v_radial{T},
        IMASdd.edge_transport__model___ggd___total_ion_energy__d{T},
        IMASdd.edge_transport__model___ggd___total_ion_energy__d_pol{T},
        IMASdd.edge_transport__model___ggd___total_ion_energy__d_radial{T},
        IMASdd.edge_transport__model___ggd___total_ion_energy__flux{T},
        IMASdd.edge_transport__model___ggd___total_ion_energy__flux_limiter{T},
        IMASdd.edge_transport__model___ggd___total_ion_energy__flux_pol{T},
        IMASdd.edge_transport__model___ggd___total_ion_energy__flux_radial{T},
        IMASdd.edge_transport__model___ggd___total_ion_energy__v{T},
        IMASdd.edge_transport__model___ggd___total_ion_energy__v_pol{T},
        IMASdd.edge_transport__model___ggd___total_ion_energy__v_radial{T},
        IMASdd.edge_transport__model___ggd_fast___electrons__particle_flux_integrated{
            T,
        },
        IMASdd.edge_transport__model___ggd_fast___electrons__power{T},
        IMASdd.edge_transport__model___ggd_fast___energy_flux_max{T},
        IMASdd.edge_transport__model___ggd_fast___ion___particle_flux_integrated{T},
        IMASdd.edge_transport__model___ggd_fast___neutral___particle_flux_integrated{T},
        IMASdd.edge_transport__model___ggd_fast___power{T},
        IMASdd.edge_transport__model___ggd_fast___power_ion_total{T},
        IMASdd.equilibrium__time_slice___ggd___b_field_r{T},
        IMASdd.equilibrium__time_slice___ggd___b_field_tor{T},
        IMASdd.equilibrium__time_slice___ggd___b_field_z{T},
        IMASdd.equilibrium__time_slice___ggd___j_parallel{T},
        IMASdd.equilibrium__time_slice___ggd___j_tor{T},
        IMASdd.equilibrium__time_slice___ggd___phi{T},
        IMASdd.equilibrium__time_slice___ggd___psi{T},
        IMASdd.equilibrium__time_slice___ggd___r{T},
        IMASdd.equilibrium__time_slice___ggd___theta{T},
        IMASdd.equilibrium__time_slice___ggd___z{T},
        IMASdd.mhd__ggd___a_field_r{T},
        IMASdd.mhd__ggd___a_field_tor{T},
        IMASdd.mhd__ggd___a_field_z{T},
        IMASdd.mhd__ggd___b_field_r{T},
        IMASdd.mhd__ggd___b_field_tor{T},
        IMASdd.mhd__ggd___b_field_z{T},
        IMASdd.mhd__ggd___electrons__temperature{T},
        IMASdd.mhd__ggd___j_r{T},
        IMASdd.mhd__ggd___j_tor{T},
        IMASdd.mhd__ggd___j_tor_r{T},
        IMASdd.mhd__ggd___j_z{T},
        IMASdd.mhd__ggd___mass_density{T},
        IMASdd.mhd__ggd___n_i_total{T},
        IMASdd.mhd__ggd___phi_potential{T},
        IMASdd.mhd__ggd___psi{T},
        IMASdd.mhd__ggd___t_i_average{T},
        IMASdd.mhd__ggd___velocity_parallel{T},
        IMASdd.mhd__ggd___velocity_parallel_over_b_field{T},
        IMASdd.mhd__ggd___velocity_r{T},
        IMASdd.mhd__ggd___velocity_tor{T},
        IMASdd.mhd__ggd___velocity_z{T},
        IMASdd.mhd__ggd___vorticity{T},
        IMASdd.mhd__ggd___vorticity_over_r{T},
        IMASdd.mhd__ggd___zeff{T},
        IMASdd.radiation__process___ggd___electrons__emissivity{T},
        IMASdd.radiation__process___ggd___ion___emissivity{T},
        IMASdd.radiation__process___ggd___ion___state___emissivity{T},
        IMASdd.radiation__process___ggd___neutral___emissivity{T},
        IMASdd.radiation__process___ggd___neutral___state___emissivity{T},
        IMASdd.runaway_electrons__distribution__ggd___expansion___grid_subset{T},
        IMASdd.runaway_electrons__distribution__ggd___expansion_fd3v___grid_subset{T},
        IMASdd.runaway_electrons__ggd_fluid___current_density{T},
        IMASdd.runaway_electrons__ggd_fluid___ddensity_dt_compton{T},
        IMASdd.runaway_electrons__ggd_fluid___ddensity_dt_dreicer{T},
        IMASdd.runaway_electrons__ggd_fluid___ddensity_dt_hot_tail{T},
        IMASdd.runaway_electrons__ggd_fluid___ddensity_dt_total{T},
        IMASdd.runaway_electrons__ggd_fluid___ddensity_dt_tritium{T},
        IMASdd.runaway_electrons__ggd_fluid___density{T},
        IMASdd.runaway_electrons__ggd_fluid___e_field_critical{T},
        IMASdd.runaway_electrons__ggd_fluid___e_field_dreicer{T},
        IMASdd.runaway_electrons__ggd_fluid___energy_density_kinetic{T},
        IMASdd.runaway_electrons__ggd_fluid___momentum_critical_avalanche{T},
        IMASdd.runaway_electrons__ggd_fluid___momentum_critical_hot_tail{T},
        IMASdd.runaway_electrons__ggd_fluid___pitch_angle{T},
        IMASdd.tf__field_map___a_field_r{T},
        IMASdd.tf__field_map___a_field_tor{T},
        IMASdd.tf__field_map___a_field_z{T},
        IMASdd.tf__field_map___b_field_r{T},
        IMASdd.tf__field_map___b_field_tor{T},
        IMASdd.tf__field_map___b_field_z{T},
        IMASdd.transport_solver_numerics__boundary_conditions_ggd___current{T},
        IMASdd.transport_solver_numerics__boundary_conditions_ggd___electrons__energy{
            T,
        },
        IMASdd.transport_solver_numerics__boundary_conditions_ggd___electrons__particles{
            T,
        },
        IMASdd.transport_solver_numerics__boundary_conditions_ggd___ion___energy{T},
        IMASdd.transport_solver_numerics__boundary_conditions_ggd___ion___particles{T},
        IMASdd.transport_solver_numerics__boundary_conditions_ggd___ion___state___energy{
            T,
        },
        IMASdd.transport_solver_numerics__boundary_conditions_ggd___ion___state___particles{
            T,
        },
        IMASdd.wall__description_ggd___component___type{T},
        IMASdd.wall__description_ggd___ggd___a_field{T},
        IMASdd.wall__description_ggd___ggd___e_field{T},
        IMASdd.wall__description_ggd___ggd___energy_fluxes__current__emitted{T},
        IMASdd.wall__description_ggd___ggd___energy_fluxes__current__incident{T},
        IMASdd.wall__description_ggd___ggd___energy_fluxes__kinetic__electrons__emitted{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__kinetic__electrons__incident{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__kinetic__ion___emitted{T},
        IMASdd.wall__description_ggd___ggd___energy_fluxes__kinetic__ion___incident{T},
        IMASdd.wall__description_ggd___ggd___energy_fluxes__kinetic__ion___state___emitted{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__kinetic__ion___state___incident{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__kinetic__neutral___emitted{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__kinetic__neutral___incident{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__kinetic__neutral___state___emitted{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__kinetic__neutral___state___incident{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__radiation__emitted{T},
        IMASdd.wall__description_ggd___ggd___energy_fluxes__radiation__incident{T},
        IMASdd.wall__description_ggd___ggd___energy_fluxes__recombination__ion___emitted{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__recombination__ion___incident{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__recombination__ion___state___emitted{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__recombination__ion___state___incident{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__recombination__neutral___emitted{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__recombination__neutral___incident{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__recombination__neutral___state___emitted{
            T,
        },
        IMASdd.wall__description_ggd___ggd___energy_fluxes__recombination__neutral___state___incident{
            T,
        },
        IMASdd.wall__description_ggd___ggd___j_total{T},
        IMASdd.wall__description_ggd___ggd___particle_fluxes__electrons__emitted{T},
        IMASdd.wall__description_ggd___ggd___particle_fluxes__electrons__incident{T},
        IMASdd.wall__description_ggd___ggd___particle_fluxes__ion___emitted{T},
        IMASdd.wall__description_ggd___ggd___particle_fluxes__ion___incident{T},
        IMASdd.wall__description_ggd___ggd___particle_fluxes__ion___state___emitted{T},
        IMASdd.wall__description_ggd___ggd___particle_fluxes__ion___state___incident{T},
        IMASdd.wall__description_ggd___ggd___particle_fluxes__neutral___emitted{T},
        IMASdd.wall__description_ggd___ggd___particle_fluxes__neutral___incident{T},
        IMASdd.wall__description_ggd___ggd___particle_fluxes__neutral___state___emitted{
            T,
        },
        IMASdd.wall__description_ggd___ggd___particle_fluxes__neutral___state___incident{
            T,
        },
        IMASdd.wall__description_ggd___ggd___phi_potential{T},
        IMASdd.wall__description_ggd___ggd___power_density{T},
        IMASdd.wall__description_ggd___ggd___psi{T},
        IMASdd.wall__description_ggd___ggd___recycling__ion___coefficient{T},
        IMASdd.wall__description_ggd___ggd___recycling__ion___state___coefficient{T},
        IMASdd.wall__description_ggd___ggd___recycling__neutral___coefficient{T},
        IMASdd.wall__description_ggd___ggd___recycling__neutral___state___coefficient{
            T,
        },
        IMASdd.wall__description_ggd___ggd___resistivity{T},
        IMASdd.wall__description_ggd___ggd___temperature{T},
        IMASdd.wall__description_ggd___ggd___v_biasing{T},
        IMASdd.wall__description_ggd___material___grid_subset{T},
        IMASdd.wall__description_ggd___thickness___grid_subset{T},
        IMASdd.waves__coherent_wave___full_wave___b_field__bi_normal{T},
        IMASdd.waves__coherent_wave___full_wave___b_field__normal{T},
        IMASdd.waves__coherent_wave___full_wave___b_field__parallel{T},
        IMASdd.waves__coherent_wave___full_wave___e_field__bi_normal{T},
        IMASdd.waves__coherent_wave___full_wave___e_field__minus{T},
        IMASdd.waves__coherent_wave___full_wave___e_field__normal{T},
        IMASdd.waves__coherent_wave___full_wave___e_field__parallel{T},
        IMASdd.waves__coherent_wave___full_wave___e_field__plus{T},
        IMASdd.waves__coherent_wave___full_wave___k_perpendicular{T},
    } where {T};

# NOTE: `Base.getproperty(::all__grid_ggd, ::Symbol)` (grid_ggd path-linking)
# moved to IMASdd (>= 8.6.0), where it belongs: it dispatches on IMASdd types.
