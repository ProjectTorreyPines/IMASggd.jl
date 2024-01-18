import NearestNeighbors: KDTree, knn
import StaticArrays: SVector
import Statistics: mean
import Interpolations: linear_interpolation

export interp
export get_kdtree
export get_TPS_mats

function get_kdtree(space::OMAS.edge_profiles__grid_ggd___space)
    grid_nodes = space.objects_per_dimension[1].object
    grid_faces = space.objects_per_dimension[3].object
    grid_faces = [cell for cell ∈ grid_faces if length(cell.nodes) == 4]
    grid_centers = [
        SVector{2}(mean([grid_nodes[node].geometry for node ∈ cell.nodes])) for
        cell ∈ grid_faces
    ]
    return KDTree(grid_centers; leafsize=10)
end

function get_kdtree(
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset,
)
    subset_centers = get_subset_centers(space, subset)
    return KDTree([SVector{2}(sc) for sc ∈ subset_centers]; leafsize=10)
end

"""
    interp(
    prop_values::Vector{T},
    kdtree::KDTree;
    use_nearest_n::Int=4,
    weighing::Function=(d) -> 1 / d,

) where {T <: Real}

Lowest level interpolation function. It takes a vector of property values and a KDTree
defined over a 2D space with the same number of nodes as the property values. It returns
a function that can be used to interpolate the property values at any point in the
space.
"""
function interp(
    prop_values::Vector{T},
    kdtree::KDTree;
    use_nearest_n::Int=4,
    weighing::Function=(d) -> 1 / d,
) where {T <: Real}
    function get_interp_val(x::Real, y::Real)
        nearest_indices, distances = knn(kdtree, Array([x, y]), use_nearest_n)
        values = [prop_values[ii] for ii ∈ nearest_indices]
        weights = weighing.(distances)
        if any(isinf.(weights))
            return values[distances.==0][1]
        end
        return sum(weights .* values) / sum(weights)
    end
    get_interp_val(xy::Tuple{Real, Real}) = get_interp_val(xy...)
    return get_interp_val
end

"""
    _G(x1::Tuple{U, U}, x2::Tuple{U, U}) where {U <: Real}

Helper function for the interpolation function. It calculates the Green's function for
minimizing bending energy of a surface.

http://www.geometrictools.com/Documentation/ThinPlateSplines.pdf Eq(28)
"""
function _G(x1::Tuple{U, U}, x2::Tuple{U, U}) where {U <: Real}
    r2 = sum((x1 .- x2) .^ 2) # r²
    # Note this uses log(r) / 8π = log(r²) / 16π
    return (r2 == 0) ? zero(U) : r2 * log(r2) * inv_16pi
end

"""
    _condition_y(y::Vector{T}) where {T <: Real}

Conditioning function on value vector before interpolation is performed. If the values
range from negative to positive range, then linear condition is done in which the mean
is subtracted and the range of the vector is used to normalize the values so they lie
between -1 and 1. If the values all have the same sign and they vary by more than 2
orders of magnitude, then log10 is taken after normalizing the values with their
minimum absolute value. This is done to avoid numerical issues with interpolation.
Return values are conditioned y and inverse conditioning function.
"""
function _condition_y(y::Vector{T}) where {T <: Real}
    ylims = extrema(y)
    do_log = (prod(ylims) > 0) && (ylims[2] / ylims[1] > 100)
    norm_by = do_log ? minimum(abs.(ylims)) * sign(ylims[1]) : ylims[2] - ylims[1]
    mean_y = mean(y)
    cy = do_log ? log10.(y ./ norm_by) : (y .- mean_y) ./ norm_by
    inv_cy = let norm_by = norm_by, do_log = do_log, mean_y = mean_y
        x -> do_log ? (10.0^x) * norm_by : (x * norm_by) + mean_y
    end
    return cy, inv_cy
end

"""
    get_TPS_mats(x::Vector{Tuple{U, U}}) where {U <: Real}

We can separate matrix operations of thin plate spline method to utilize
the calculation for itnerpolating over same space. Similar to the idea of
reusing KDTree like above.
"""
function get_TPS_mats(x::Vector{Tuple{U, U}}) where {U <: Real}
    # Setup matrices as defined in Eq (30)
    M = Matrix{U}(undef, length(x), length(x))
    for ii ∈ eachindex(x)
        for jj ∈ eachindex(x)
            M[ii, jj] = _G(x[ii], x[jj])
        end
    end
    N = Matrix{U}(undef, length(x), 3)
    for ii ∈ eachindex(x)
        N[ii, 1] = 1.0
        N[ii, 2] = x[ii][1]
        N[ii, 3] = x[ii][2]
    end
    Minv = M^(-1)
    return Minv, N, (N' * Minv * N)^(-1) * N' * Minv, x
end

function get_interp_val(r, z, x, a, b, inv_cy)
    tot = sum(a[k] * _G((r, z), x[k]) for k ∈ eachindex(a))
    tot += b[1] + r * b[2] + z * b[3]
    return inv_cy(tot)
end

"""
    interp(
    y::Vector{T},
    TPS_mats::Tuple{Matrix{U}, Matrix{U}, Matrix{U}, Vector{Tuple{U, U}}},

) where {T <: Real, U <: Real}

Lowest level function for Thin Plate Spline method
"""
function interp(
    y::Vector{T},
    TPS_mats::Tuple{Matrix{U}, Matrix{U}, Matrix{U}, Vector{Tuple{U, U}}},
) where {T <: Real, U <: Real}
    Minv, N, y2b, x = TPS_mats
    cy, inv_cy = _condition_y(y)
    # From Eq(31)
    b = y2b * cy
    a = Minv * (cy - N * b)
    return let x = x, a = a, b = b, inv_cy = inv_cy
        g(r, z) = get_interp_val(r, z, x, a, b, inv_cy)
        g(gp::Tuple{V, V}) where {V <: Real} = g(gp...)
    end
end

"""
    interp(y::Vector{T}, x::Vector{Tuple{U, U}}) where {T <: Real, U <: Real}

Thin plate smoothing interpolation function for a 2d space scalar function. The
algorithm has been adopted from:

http://www.geometrictools.com/Documentation/ThinPlateSplines.pdf

This is an implementation of Euler-Lagrange equation for minimizing bending energy of a
surface.
"""
function interp(y::Vector{T}, x::Vector{Tuple{U, U}}) where {T <: Real, U <: Real}
    length(x) == length(y) || error("Space (r, z) and values must have the same length")
    return interp(y, get_TPS_mats(x))
end

function get_TPS_mats(space::OMAS.edge_profiles__grid_ggd___space)
    nodes = [Tuple(node.geometry) for node ∈ space.objects_per_dimension[1].object]
    return get_TPS_mats(nodes)
end

"""
    interp(
    prop_values::Vector{T},
    space::OMAS.edge_profiles__grid_ggd___space

) where {T <: Real}

If the whole space is provided instead of a kdtree, calculate the kdtree for whole
space. Again, here it is assumed that the property values are porvided for each node
of the space.
"""
function interp(
    prop_values::Vector{T},
    space::OMAS.edge_profiles__grid_ggd___space,
) where {T <: Real}
    return interp(prop_values, get_TPS_mats(space))
end

function get_TPS_mats(
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset,
)
    return get_TPS_mats(get_subset_centers(space, subset))
end

"""
    interp(
    prop_values::Vector{Real},
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset

)

If a subset of the space is provided, calculate the kdtree for the subset. In this case
it is assumed that the property values are provided for each element of the subset.
"""
function interp(
    prop_values::Vector{T},
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset,
) where {T <: Real}
    return interp(prop_values, get_TPS_mats(space, subset))
end

"""
    interp(
    prop::edge_profiles__prop_on_subset,
    grid_ggd::OMAS.edge_profiles__grid_ggd,
    value_field::Symbol=:values

)

Example:
grid_ggd = dd.edge_profiles.grid_ggd[1]
get_electron_density = interp(dd.edge_profiles.ggd[1].electrons.density[1], grid_ggd)
get_e_field_par = interp(dd.edge_profiles.ggd[1].e_field[1], grid_ggd, :parallel)
"""
function interp(
    prop::edge_profiles__prop_on_subset,
    grid_ggd::OMAS.edge_profiles__grid_ggd,
    value_field::Symbol=:values,
)
    subset = get_grid_subset_with_index(grid_ggd, prop.grid_subset_index)
    space = grid_ggd.space[subset.element[1].object[1].space]
    return interp(getfield(prop, value_field), space, subset)
end

"""
    interp(
    prop_arr::Vector{T},
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset,
    value_field::Symbol=:values

) where {T <: edge_profiles__prop_on_subset}

Example:
sol = get_grid_subset_with_index(dd.edge_profiles.grid_ggd[1], 23)
get_electron_density = interp(dd.edge_profiles.ggd[1].electrons.density, space, sol)
"""
function interp(
    prop_arr::Vector{T},
    space::OMAS.edge_profiles__grid_ggd___space,
    subset::OMAS.edge_profiles__grid_ggd___grid_subset,
    value_field::Symbol=:values,
) where {T <: edge_profiles__prop_on_subset}
    prop = get_prop_with_grid_subset_index(prop_arr, subset.identifier.index)
    return interp(getfield(prop, value_field), space, subset)
end

"""
    interp(
    prop_arr::Vector{T},
    grid_ggd::OMAS.edge_profiles__grid_ggd,
    grid_subset_index::Int,
    value_field::Symbol=:values

) where {T <: edge_profiles__prop_on_subset}

Example:
get_n_e_sep = interp(dd.edge_profiles.ggd[1].electrons.density, grid_ggd, 16)
"""
function interp(
    prop_arr::Vector{T},
    grid_ggd::OMAS.edge_profiles__grid_ggd,
    grid_subset_index::Int,
    value_field::Symbol=:values,
) where {T <: edge_profiles__prop_on_subset}
    prop = get_prop_with_grid_subset_index(prop_arr, grid_subset_index)
    subset = get_grid_subset_with_index(grid_ggd, grid_subset_index)
    space = grid_ggd.space[subset.element[1].object[1].space]
    return interp(getfield(prop, value_field), space, subset)
end

function get_TPS_mats(grid_ggd::OMAS.edge_profiles__grid_ggd, grid_subset_index::Int)
    subset = get_grid_subset_with_index(grid_ggd, grid_subset_index)
    space = grid_ggd.space[subset.element[1].object[1].space]
    return get_TPS_mats(space, subset)
end

#! format off
"""
    interp(
    prop_arr::Vector{T},
    TPS_mats::Tuple{Matrix{U}, Matrix{U}, Matrix{U}, Vector{Tuple{U, U}}},
    grid_subset_index::Int,
    value_field::Symbol=:values,

) where {T <: edge_profiles__prop_on_subset}

Same use case as above but allows one to reuse previously calculated TPS matrices.

Example:
TPS_mat = get_TPS_mats(dd.edge_profiles.grid_ggd[1], 5)

for it ∈ eachindex(dd.edge_profiles.ggd)
    get_n_e = interp(dd.edge_profiles.ggd[it].electrons.density, TPS_mat_sep, 5)
    println("This time step has n_e at (0, 0) = ", get_n_e(, 0))
end

This will run faster has heavy matrix calculations will happen only once.
"""
#! format on
function interp(
    prop_arr::Vector{T},
    TPS_mats::Tuple{Matrix{U}, Matrix{U}, Matrix{U}, Vector{Tuple{U, U}}},
    grid_subset_index::Int,
    value_field::Val{V}=Val(:values),
) where {T <: edge_profiles__prop_on_subset, U <: Real, V}
    prop = get_prop_with_grid_subset_index(prop_arr, grid_subset_index)
    field = getfield(prop, V)
    return interp(field, TPS_mats)
end

const RHO_EXT_POS = [1.0001, 1.1, 5]
const RHO_EXT_NEG = [-5, -0.0001] # I guess this would be at a PF coil or something?

"""
    interp(eqt::OMAS.equilibrium__time_slice)

For a given equilibrium time slice, return a function that can be used to interpolate
from (r, z) space to rho (normalized toroidal flux coordinate)space.

Example:
rz2rho = interp(dd.equilibrium.time_slice[1])
rho = rz2rho.([(r, z) for r in 3:0.01:9, for z in -5:0.01:5])
"""
function interp(eqt::OMAS.equilibrium__time_slice)
    p1 = eqt.profiles_1d
    p2 = eqt.profiles_2d[1]
    gq = eqt.global_quantities
    psi_a = gq.psi_axis
    psi_b = gq.psi_boundary
    rhon_eq = p1.rho_tor_norm
    psi_eq = p1.psi
    psin_eq = (psi_eq .- psi_a) ./ (psi_b - psi_a)
    psirz = p2.psi
    psinrz = (psirz .- psi_a) ./ (psi_b - psi_a)
    r_eq = p2.grid.dim1
    z_eq = p2.grid.dim2
    # rho_N isn't defined on open flux surfaces, so it is extended by copying psi_N
    psin_eq_ext = copy(psin_eq)
    append!(psin_eq_ext, RHO_EXT_POS)
    rhon_eq_ext = copy(rhon_eq)
    append!(rhon_eq_ext, RHO_EXT_POS)
    prepend!(psin_eq_ext, RHO_EXT_NEG)
    prepend!(rhon_eq_ext, RHO_EXT_NEG)
    rz2psin = linear_interpolation((r_eq, z_eq), psinrz)
    psin2rhon = linear_interpolation(psin_eq_ext, rhon_eq_ext)
    g = let psin2rhon = psin2rhon, rz2psin = rz2psin
        (r, z) -> psin2rhon(rz2psin(r, z))
    end
    return g
end

"""
    interp(
    prop::Vector{T},
    prof::OMAS.core_profiles__profiles_1d,

) where {T <: Real}

Returns an inteprolation function for the core profile property values defined on
normalized toroidal flux coordinate rho.

Example:
core_profile_n_e = dd.core_profiles.profiles_1d[1].electrons.density
get_n_e = interp(core_profile_n_e, dd.core_profiles.profiles_1d[1])
get_n_e(1) # Returns electron density at rho = 1 (separatix)
"""
function interp(
    prop::Vector{T},
    prof::OMAS.core_profiles__profiles_1d,
) where {T <: Real}
    rho_prof = copy(prof.grid.rho_tor_norm)
    length(prop) == length(rho_prof) ||
        error("Property muast have same length as rho_tor_norm in core_profile")
    prepend!(rho_prof, RHO_EXT_NEG)
    append!(rho_prof, RHO_EXT_POS)
    p = copy(prop)
    prepend!(p, zeros(size(RHO_EXT_NEG)))
    append!(p, zeros(size(RHO_EXT_POS)))
    return linear_interpolation(rho_prof, p)
end

"""
    interp(
    prop::Vector{T},
    prof::OMAS.core_profiles__profiles_1d,
    rz2rho::Function,

)

Returns an inteprolation function in (R, Z) domain for the core profile property values
defined on normalized toroidal flux coordinate rho and with a provided function to
convert (R,Z) to rho.

Example:

rz2rho = interp(dd.equilibrium.time_slice[1])
core_profile_n_e = dd.core_profiles.profiles_1d[1].electrons.density
get_n_e = interp(core_profile_n_e, dd.core_profiles.profiles_1d[1], rz2rho)
get_n_e(5.0, 3.5) # Returns electron density at (R, Z) = (5.0, 3.5)
"""
function interp(
    prop::Vector{T},
    prof::OMAS.core_profiles__profiles_1d,
    rz2rho::Function,
) where {T <: Real}
    itp = interp(prop, prof)
    g = let itp = itp
        (r, z) -> itp(rz2rho(r, z))
    end
    return g
end

"""
    interp(
    prop::Vector{T},
    prof::OMAS.core_profiles__profiles_1d,
    eqt::OMAS.equilibrium__time_slice,

) where {T <: Real}

Returns an inteprolation function in (R, Z) domain for the core profile property values
defined on normalized toroidal flux coordinate rho and with a provided equilibrium time
slice to get (R, Z) to rho conversion.

Example:

eqt = dd.equilibrium.time_slice[1]
core_profile_n_e = dd.core_profiles.profiles_1d[1].electrons.density
get_n_e = interp(core_profile_n_e, dd.core_profiles.profiles_1d[1], eqt)
get_n_e(5.0, 3.5) # Returns electron density at (R, Z) = (5.0, 3.5)
"""
function interp(
    prop::Vector{T},
    prof::OMAS.core_profiles__profiles_1d,
    eqt::OMAS.equilibrium__time_slice,
) where {T <: Real}
    rz2rho = interp(eqt)
    return interp(prop, prof, rz2rho)
end
