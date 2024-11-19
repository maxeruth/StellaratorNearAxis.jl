abstract type AbstractPowerSeries{T} end
abstract type AbstractPowerSeriesSlice{T} end

## Structs
"""
    PowerSeriesRho

A struct used for algebraically representing changes in order of power series, typically stored as `p0`.

# Example
```
A = zero_SpectralPowerSeries(1,1)
A[1].a[1] = 1.         # A is a constant SpectralPowerSeries
display(get_p0(A))     # The leading order of A is 0
ρ = PowerSeriesRho()
B = ρ^2*A
display(get_p0(B))     # The SpectralPowerSeries B is ρ^2
```
"""
struct PowerSeriesRho
    p0::Vector{Integer}

    function PowerSeriesRho(; p0::Integer=1)
        new([p0])
    end
end

struct ZeroPowerSeriesSlice;
    function ZeroPowerSeriesSlice()
        new()
    end
end

# Used to wrap a given order of the expansion with a ``slice'' object
# If view_slice is used to obtain the slice, there is no memory allocation,
# and editing the underlying data will edit the original array.
# If get_slice is used, the data is copied, and you can edit without fear!
struct SpectralPowerSeriesSlice{T} <: AbstractPowerSeriesSlice{T}
    a::AbstractArray{T}

    function SpectralPowerSeriesSlice(a::AbstractArray{T}) where {T}
        new{T}(a)
    end
end



# Used to wrap a given order of the expansion with a ``slice'' object
# If view_slice is used to obtain the slice, there is no memory allocation,
# and editing the underlying data will edit the original array.
# If get_slice is used, the data is copied, and you can edit without fear!
struct SpatialPowerSeriesSlice{T} <: AbstractPowerSeriesSlice{T}
    a::AbstractArray{T}

    function SpatialPowerSeriesSlice(a::AbstractArray{T}) where T
        new{T}(a)
    end
end

"""
    ZeroPowerSeries{T} <: AbstractPowerSeries{T}

A power series representing the zero element. 
Multiplication of an AbstractPowerSeries by a ZeroPowerSeries returns a ZeroPowerSeries, and
addition of an AbstractPowerSeries `A` with a ZeroPowerSeries returns `A`.

# Example
```
A = zero_SpectralPowerSeries(1,1)
A[1].a[1] = 1.
Z = ZeroPowerSeries()
display(A*Z) # Returns a ZeroPowerSeries
display(A+Z) # Returns A
```
"""
struct ZeroPowerSeries{T} <: AbstractPowerSeries{T}
    function ZeroPowerSeries(;T::DataType=Float64)
        new{T}()
    end
end

"""
    IdentityPowerSeries{T} <: AbstractPowerSeries{T}

A power series representing the identity element. 
Multiplication of an AbstractPowerSeries `A` by an IdentityPowerSeries returns `A`.

# Example
```
A = zero_SpectralPowerSeries(1,1)
A[1].a[1] = 2.
id = IdentityPowerSeries()
display(A*id) # Returns A
```
"""
struct IdentityPowerSeries{T} <: AbstractPowerSeries{T}
    function IdentityPowerSeries(; T::DataType=Float64)
        new{T}()
    end
end

"""
    FluxPowerSeries{T} <: AbstractPowerSeries{T}

Formal power series that is a function of the flux surface, i.e. `A = A₀ + A₁ (ρ²)¹ + A₂ (ρ²)² + ...`
"""
struct FluxPowerSeries{T} <: AbstractPowerSeries{T}
    a::AbstractVector{T}; # Coefficients
    p0::Vector{Integer};  # Extra factor of p0

    function FluxPowerSeries(a::AbstractVector{T}; p0::Integer=0) where T
        new{T}(a, [p0])
    end
end

"""
    SpectralPowerSeries{T} <: AbstractPowerSeries{T}

The basic spectral representation of near-axis series.
A `SpectralPowerSeries` is composed of `SpectralPowerSeriesSlice`s, each of which can be accessed by getindex via `A[n]`.
Each slice `n` contains an `M × n` array of Fourier coefficients, where 
`M` is the number of coefficients in the `s` direction and `n` is the number of coefficients in the `θ` direction. 
In the `s` direction, the coefficients are ordered as `1`, `sin(s)`, `cos(s)`, `sin(2s)`, ...
In the odd `n` case in the `θ` direction, the coefficients are ordered as `1`, `sin(2θ)`, `cos(2θ)`, `sin(4θ)`, ...
In the even `n` case in the `θ` direction, the coefficients are ordered as `cos(θ)`, `sin(θ)`, `cos(3θ)`, `sin(3θ)`, ...
Derivative routines are defined for SpectralPowerSeries; see [`rho_deriv`](@ref), [`theta_deriv`](@ref), [`s_deriv`](@ref), [`grad`](@ref), [`div`](@ref).

# Example
Create a power series representing `y = ρ sin(θ)`
```
M = 5                              # Number of modes in s
N = 3                              # Number of orders in rho
y = zero_SpectralPowerSeries(M, N) # Create a SpectralPowerSeries
y[2].a[1,2] = 1.                   # Set the ρ sin(θ) coefficient to 1
y
```
"""
struct SpectralPowerSeries{T} <: AbstractPowerSeries{T}
    a::AbstractVector{SpectralPowerSeriesSlice{T}}; # Fourier coefficients, layed out as
      # a = [[A0c0],  [A1c1,A1s1],  [A2c0,A2c2,A2s2],  [A3c1,A3s1,A3c3,A3s3], ...],
      # where each Ancm and Ansm are columns where, e.g.
      # A0c0(ϕ) = A0c0[1] + A0c0[2] cos(ϕ) + A0c0[3] sin(ϕ) + A0c0[4] cos(2ϕ) + ...
      # where there are 2M+1 coefficients
      # TODO: consider using a rfft instead

    p0::Vector{Integer};
    function SpectralPowerSeries(a::AbstractVector{SpectralPowerSeriesSlice{T}}; p0::Integer=0) where {T}
        new{T}(a, [p0])
    end
end



"""
    SpatialPowerSeries{T} <: AbstractPowerSeries{T}

The representation of near-axis series on collocation nodes.
A `SpatialPowerSeries` is composed of `SpatialPowerSeriesSlice`s, each of which can be accessed by getindex via `A[n]`.
Each slice `n` contains an `M × n` array of sampled values, where 
`M` is the number of samples in the `s` direction and `n` is the order in the `θ` direction. 
The `s` points are ordered as `2(j-1)π/M` for `1≤j≤M`, while in the `θ` points are ordered as `2π(k-1)/(2n-1)` for `1≤k≤n`.
Algebraic operations, such as `*`, `/`, `inv`, etc. are defined for SpatialPowerSeries

# Example
Create a power series representing `y = ρ sin(θ)`
```
Mc = 5                             # Number of points in s
N = 3                              # Number of orders in rho
y = zero_SpectralPowerSeries(1, N) # Create a SpectralPowerSeries
y[2].a[1,2] = 1.                   # Set the ρ sin(θ) coefficient to 1
y_spatial = SpatialPowerSeries(y, M=Mc)
y_spatial
```
"""
struct SpatialPowerSeries{T} <: AbstractPowerSeries{T}
    a::AbstractArray{SpatialPowerSeriesSlice{T}};
      # Values, layed out as
      # a = [[A0c0], [A1c1,A1s1], [A2c0,A2c2,A2s2], [A3c1,A3s1,A3c3,A3s3], ...],
      # where each Ancm and Ansm are columns of equisampled points
      # TODO: consider using a rfft instead

    p0::AbstractVector{Integer}; # Leading ϵ power

    function SpatialPowerSeries(a::AbstractArray{SpatialPowerSeriesSlice{T}};
                                p0::Integer=0) where T
        new{T}(a, [p0])
    end
end


## ZeroPowerSeriesSlice functions
function *(A::ZeroPowerSeriesSlice, B::AbstractPowerSeriesSlice)
    ZeroPowerSeriesSlice();
end
function *(A::AbstractPowerSeriesSlice, B::ZeroPowerSeriesSlice)
    ZeroPowerSeriesSlice();
end
function *(A::ZeroPowerSeriesSlice, B::ZeroPowerSeriesSlice)
    ZeroPowerSeriesSlice();
end
function +(A::ZeroPowerSeriesSlice, B::AbstractPowerSeriesSlice)
    B
end
function +(A::AbstractPowerSeriesSlice, B::ZeroPowerSeriesSlice)
    A
end
function +(A::ZeroPowerSeriesSlice, B::ZeroPowerSeriesSlice)
    B
end
function -(A::AbstractPowerSeriesSlice, B::ZeroPowerSeriesSlice)
    A
end

## default AbstractPowerSEriesSlice functions
Base.Broadcast.broadcastable(A::AbstractPowerSeries) = Ref(A); # Makes it so linear algebra works

function get_M(A::AbstractPowerSeriesSlice)
    size(A.a, 1);
end

function set_slice!(A::AbstractPowerSeriesSlice{T}, val::AbstractArray{T}) where {T}
    A.a[:, :] = val
end

function set_slice!(A::AbstractPowerSeriesSlice{T}, B::AbstractPowerSeriesSlice{T}) where {T}
    A.a[:, :] = B.a
end

function set_slice!(A::AbstractPowerSeriesSlice{T}, B::ZeroPowerSeriesSlice) where{T}
    A.a[:, :] .= zero(T)
end

## Include functions other PowerSeries
include("./ZeroPowerSeries.jl")
include("./IdentityPowerSeries.jl")
include("./Operations.jl")
include("./FluxPowerSeries.jl")
include("./SpatialPowerSeries.jl")
include("./SpectralPowerSeries.jl")

## PowerSeriesRho functions
function get_p0(ρ::PowerSeriesRho)
    ρ.p0[1]
end

function *(ρ::PowerSeriesRho, A::AbstractPowerSeries)
    B = deepcopy(A);
    set_p0!(B, get_p0(A) + get_p0(ρ))

    B
end

function *(A::AbstractPowerSeries, ρ::PowerSeriesRho)
    B = deepcopy(A);
    set_p0!(B, get_p0(A) + get_p0(ρ))

    B
end

function ^(ρ::PowerSeriesRho, n::Integer)
    PowerSeriesRho(;p0=get_p0(ρ)*n)
end

function inv(ρ::PowerSeriesRho)
    PowerSeriesRho(;p0=-get_p0(ρ))
end

function /(A::AbstractPowerSeries{T}, ρ::PowerSeriesRho) where {T}
    A*inv(ρ)
end

"""
    change_order(A::AbstractPowerSeries{T}, N::Integer) where {T}

Creates a new PowerSeries with the coefficients of `A` with 
order `N`. If `N > get_N(A)`, the extra coefficients are set to zero.
"""
function change_order(A::AbstractPowerSeries{T}, N::Integer) where {T}
    B = similar(A; N)
    M = min(N, get_N(A))
    for ii = 1:M
        B[ii] = A[ii]
    end

    B
end
