# Power Series
This is the documentation for power series structs and their operations.

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
    SpectralPowerSeries
    SpatialPowerSeries
```

## Constructors
```@docs
    zero_SpectralPowerSeries
    zero_SpatialPowerSeries
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
    rho_antideriv
    volume_integrate
    LinearAlgebra.norm(A::Vector{S}) where {S <: AbstractPowerSeries}
    Fnorm
```

## Helper Functions
```@docs
    remove_zeros
    
```