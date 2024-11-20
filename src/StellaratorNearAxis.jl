module StellaratorNearAxis
import Base: *, /, -, +, ^, inv, exp, sincos, sin, cos, sinh, cosh, adjoint, similar
import Base: getindex, setindex!, zero, div
import Base: show
import LinearAlgebra: norm, cross, *, dot
using LinearAlgebra
using StaticArrays
using LoopVectorization
using OrdinaryDiffEq
using Requires
using Printf
using Memoization
using NLsolve
using JSON

include("./PowerSeries/PowerSeries.jl")
include("./DirectFS/DirectFS.jl")
include("./BiotSavart/BiotSavart.jl")


# function __init__()
#     @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" begin
        using CairoMakie
        include("./Plotting/SNA_CairoMakie.jl")
        export plot_curve!, flux_surface_plot!, flux_contours!
#     end
# end

export FluxPowerSeries, get_N, get_a, zero_FluxPowerSeries, flux_compose, fourier_matrix, 
    zero_SpatialPowerSeries, zero_SpectralPowerSeries, set_slice!, SpectralPowerSeries, 
    SpatialPowerSeries, view_slice, get_slice, composition_basis, compose, invert_coordinates, 
    s_deriv, theta_deriv, get_p0, Fnorm, distribute_p0, unsafe_distribute_p0, L2norm, 
    PowerSeriesRho, AbstractPowerSeries, ZeroPowerSeries, IdentityPowerSeries, vacuum_update_B, 
    vacuum_residual, grad, div, rho_deriv, get_IC_from_eq, vacuum_solve, get_flux_coordinates, 
    change_order, remove_zeros, evaluate, section, section_composition_basis, section_compose, 
    fit_SpectralPowerSeries, surface_integrate, rho_integrate, rho_antideriv, volume_integrate

export DirectNearAxisEquilibrium, InitialVacuumNearAxisEquilibrium, r0_to_r, to_arclength, 
    magnetic_trajectory

export get_positions_B

export Coil, load_coils, find_magnetic_axis, magnetic_map, axis_from_point, get_field_on_axis, 
    field_to_nae


end # module StellaratorNearAxis
