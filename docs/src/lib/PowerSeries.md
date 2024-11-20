# Documentation
This is the documentation for StellaratorNearAxis.jl structs and their operations.

## Contents
```@contents
Pages = [
    "PowerSeries.md"
]
```

## Index
```@index
Pages = [
    "PowerSeries.md"
]
```

## Structs
```@docs
    PowerSeriesRho
    ZeroPowerSeries
    IdentityPowerSeries
    FluxPowerSeries
    SpectralPowerSeries
    SpatialPowerSeries
```

## Constructors
```@docs
    zero_SpectralPowerSeries
    zero_SpatialPowerSeries
    zero_FluxPowerSeries
    similar(::SpatialPowerSeries{T}) where {T}
    SpatialPowerSeries(::SpectralPowerSeries{T}) where {T}
    SpectralPowerSeries(::SpatialPowerSeries{T}) where {T}
    fit_SpectralPowerSeries
```

## Operations
```@docs
    evaluate(::SpectralPowerSeries,::Union{Number, AbstractVector},::Union{Number, AbstractVector},::Union{Number, AbstractVector})
    -(::AbstractPowerSeries)
    +(::AbstractPowerSeries,::AbstractPowerSeries)
    -(::AbstractPowerSeries,::AbstractPowerSeries)
    *(::AbstractPowerSeries{T},::AbstractPowerSeries{T}) where {T}
    *(::Number, ::AbstractPowerSeries)
    inv(::AbstractPowerSeries)
    /(A::AbstractPowerSeries{T}, B::AbstractPowerSeries{T}) where {T}
    exp(::AbstractPowerSeries{T},::T) where {T}
    sinh(::AbstractPowerSeries{T},::T) where {T}
    cosh(::AbstractPowerSeries{T},::T) where {T}
    sincos(::AbstractPowerSeries,::Number)
    sin(::AbstractPowerSeries,::Number)
    cos(::AbstractPowerSeries,::Number)
    ^(::AbstractPowerSeries,::Number)
    rho_deriv
    theta_deriv
    s_deriv
    grad
    div
    surface_integrate
    rho_integrate
    rho_antideriv(::AbstractPowerSeries)
    volume_integrate
    LinearAlgebra.norm(::AbstractVector{<:AbstractPowerSeries})
    Fnorm
    L2norm(::AbstractPowerSeries, ::Number)
    flux_compose(::FluxPowerSeries{T}, ::SpatialPowerSeries{T}) where {T}
    composition_basis
    compose
    invert_coordinates
```

## Helper Functions
```@docs
    remove_zeros(::AbstractPowerSeries)
    change_order(::AbstractPowerSeries{T}, ::Integer) where {T}
    distribute_p0(::AbstractPowerSeries, ::Integer)
    unsafe_distribute_p0
```

## Coil Functions
```@docs
    Coil
    load_coils
    evaluate(::AbstractArray{T}, ::Coil) where {T}
    magnetic_trajectory
    find_magnetic_axis
    get_field_on_axis
    field_to_nae
```

## Direct Frenet-Serret Equilibrium Problem
```@docs
    DirectNearAxisEquilibrium
    InitialVacuumNearAxisEquilibrium
    vacuum_solve
    get_flux_coordinates
```

## Plotting
```@docs
    plot_curve!
    flux_surface_plot!
    flux_contours!
```