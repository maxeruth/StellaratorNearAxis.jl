abstract type AbstractPowerSeries{T} end
abstract type AbstractPowerSeriesSlice{T} end

## Structs

struct PowerSeries_ρ
    p0::Vector{Integer}

    function PowerSeries_ρ(; p0::Integer=1)
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


struct ZeroPowerSeries{T} <: AbstractPowerSeries{T}
    function ZeroPowerSeries(;T::DataType=Float64)
        new{T}()
    end
end

struct IdentityPowerSeries{T} <: AbstractPowerSeries{T}
    function IdentityPowerSeries(; T::DataType=Float64)
        new{T}()
    end
end

## Formal power series that is a function of the flux surface
# A = A₀ + A₁ (ϵ²ρ²)¹ + A₂ (ϵ² ρ²)² + ...
struct FluxPowerSeries{T} <: AbstractPowerSeries{T}
    a::AbstractVector{T}; # Coefficients
    p0::Vector{Integer};  # Extra factor of p0

    function FluxPowerSeries(a::AbstractVector{T}; p0::Integer=0) where T
        new{T}(a, [p0])
    end
end

# Formal analytic power series
# A = A0c0(ϕ)
#     + ϵ ρ (A1c1(ϕ) cos(θ0) + A1s1(ϕ) sin(ϑ0)) +
#     + (ϵ ρ)² (A2c0(ϕ) + A2c2(ϕ) cos(2 θ0) + A2s2(ϕ) sin(2ϑ0))
#     + ...
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



# A power series function struct
# Discretization is at the collocation nodes in ϕ and θ, and asymptotic in ρ
# Collocation in θ occurs only on the upper half circle due to even/odd symmetry
# Series takes the form A = (ϵρ)^p0 (A0c0 + ϵρ(A1c1 cos θ + A1s1 sin θ) + ...)
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

## PowerSeries_ρ functions
function get_p0(ρ::PowerSeries_ρ)
    ρ.p0[1]
end

function *(ρ::PowerSeries_ρ, A::AbstractPowerSeries)
    B = deepcopy(A);
    set_p0!(B, get_p0(A) + get_p0(ρ))

    B
end

function *(A::AbstractPowerSeries, ρ::PowerSeries_ρ)
    B = deepcopy(A);
    set_p0!(B, get_p0(A) + get_p0(ρ))

    B
end

function ^(ρ::PowerSeries_ρ, n::Integer)
    PowerSeries_ρ(;p0=get_p0(ρ)*n)
end

function inv(ρ::PowerSeries_ρ)
    PowerSeries_ρ(;p0=-get_p0(ρ))
end

function /(A::AbstractPowerSeries{T}, ρ::PowerSeries_ρ) where {T}
    A*inv(ρ)
end

function change_order(A::AbstractPowerSeries{T}, N::Integer) where {T}
    B = similar(A; N)
    M = min(N, get_N(A))
    for ii = 1:M
        B[ii] = A[ii]
    end

    B
end
