## FluxPowerSeries Functions
function get_N(A::FluxPowerSeries{T}) where {T}
    length(A.a)
end

function get_a(A::Number)
    A
end

function get_a(A::FluxPowerSeries{T}) where {T}
    A.a
end

function getindex(A::FluxPowerSeries{T}, inds...) where {T}
    A.a[inds]
end

# function getindex(A::FluxPowerSeries, jj::Integer)
#     (jj <= get_N(A)) ? A.a[jj] : ZeroPowerSeriesSlice()
# end

function getindex(A::FluxPowerSeries{T}, ind::Integer) where {T}
    (ind <= get_N(A)) ? A.a[ind] : ZeroPowerSeriesSlice()
end


function setindex!(A::FluxPowerSeries, X, ind)
    A.a[ind] = X;
end

function get_M(A::FluxPowerSeries); 1; end


function zero_FluxPowerSeries(T::DataType, N::Integer; p0::Integer=0)
    FluxPowerSeries(zeros(T,N); p0)
end

"""
    zero_FluxPowerSeries(N::Integer; p0::Integer=0)

Create a FluxPowerSeries of degree `N` and leading coefficient `p0` with zero entries.
"""
function zero_FluxPowerSeries(N::Integer; p0::Integer=0)
    zero_FluxPowerSeries(Float64, N; p0)
end

function similar(A::FluxPowerSeries{T}; N::Union{Integer, Nothing}=nothing,
                 p0::Union{Integer, Nothing}=nothing) where {T}
    N =  isnothing(N)  ? get_N(A)  : N;
    p0 = isnothing(p0) ? get_p0(A) : p0;

    FluxPowerSeries(zeros(T, N); p0)
end

# function distribute_p0(A::FluxPowerSeries{T}, p0::Integer) where{T}
#     p0A = get_p0(A);
#     @assert p0 <= p0A;
#     @assert mod(p0A - p0, 2) == 0

#     offset = (p0A - p0) ÷ 2;
#     N = get_N(A) + offset
#     a = zeros(T, N);
#     a[offset+1:end] = get_a(A);

#     FluxPowerSeries(a; p0)
# end

function distribute_p0(A::FluxPowerSeries); distribute_p0(A, 0); end

function evaluate(A::FluxPowerSeries{T}, ρs::AbstractVector) where {T}
    f = zeros(T, length(ρs))
    for ii = 1:get_N(A)
        f[:] += A[ii] * ρs.^(ii-1)
    end

    return f[:] .* ρs.^get_p0(A)
end

function evaluate(A::FluxPowerSeries{T}, ρ::Number) where {T}
    f = zero(T);
    for ii = 1:get_N(A)
#         f += A[ii] * ρ^(2(ii-1))
        f += A[ii] * ρ^(ii-1) 
    end

    return f * ρ^get_p0(A)
end

function Fnorm(A::FluxPowerSeries{T}) where {T}
    return sqrt(A.a'*A.a)
end

"""
    flux_compose(A::FluxPowerSeries{T}, B::SpatialPowerSeries{T}) where {T}

Composes a FluxPowerSeries `A` with a SpatialPowerSeries `B`, i.e. `A∘B`.
Useful, e.g., for composing the rotational transform with the flux. 
"""
function flux_compose(A::FluxPowerSeries{T}, B::SpatialPowerSeries{T}) where {T}
    N = get_N(A);
    p0 = get_p0(A);
    ρ = PowerSeriesRho();

    C = similar(B; p0=p0);
    C[1].a[:] .= A[1];

    Bii = IdentityPowerSeries();
    for ii = 3:2:N
        Bii = Bii*B
        # C = C + (ρ^(p0 + ii - 1)) * (A[ii]*Bii)
        C = C + (A[ii]*Bii)
    end

    C
end

function distribute_p0!(B::FluxPowerSeries, A::FluxPowerSeries, 
                        ii::Integer, offset::Integer)
    B[ii] = A[ii-offset]
end

function unsafe_distribute_p0!(B::FluxPowerSeries, A::FluxPowerSeries, 
                               ii::Integer, offset::Integer)
    B[ii] = A[ii-offset]
end

