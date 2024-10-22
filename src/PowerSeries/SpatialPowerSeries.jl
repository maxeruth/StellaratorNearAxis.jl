## SpatialPowerSeriesSlice Functions
function get_a(A::SpatialPowerSeriesSlice)
    A.a
end

function zero_SpatialPowerSeriesSlice(M::Integer, jj::Integer)
    SpatialPowerSeriesSlice(zeros(M, jj))
end

function zero_SpatialPowerSeriesSlice(T::DataType, M::Integer, jj::Integer)
    SpatialPowerSeriesSlice(zeros(T, M, jj))
end

function LinearAlgebra.norm(A::SpatialPowerSeriesSlice)
    norm(A.a)
end

function *(a::T, A::SpatialPowerSeriesSlice{T}) where {T}
    SpatialPowerSeriesSlice(a * A.a);
end

function +(A::SpatialPowerSeriesSlice{T}, B::SpatialPowerSeriesSlice{T}) where {T}
    SpatialPowerSeriesSlice(A.a + B.a)
end

function +(A::SpatialPowerSeriesSlice{T}, a::T) where {T}
    SpatialPowerSeriesSlice(A.a .+ a)
end

function +(a::T, A::SpatialPowerSeriesSlice{T}) where {T}
    +(A, a)
end

function -(A::SpatialPowerSeriesSlice{T}, B::SpatialPowerSeriesSlice{T}) where {T}
    SpatialPowerSeriesSlice(A.a - B.a)
end

function -(A::SpatialPowerSeriesSlice{T}, a::T) where {T}
    SpatialPowerSeriesSlice(A.a .- a)
end

function -(a::T, A::SpatialPowerSeriesSlice{T}) where {T}
    -(A, a)
end

function -(A::SpatialPowerSeriesSlice)
    SpatialPowerSeriesSlice(-A.a)
end

# Make an M × N fourier matrix F
function half_fourier_matrix(x::AbstractVector, N::Integer)
    F = zeros(length(x), N)

    odd = mod(N-1, 2);

    for ii = 1:2:N
        n = 1.0 * (ii-1+odd)
        F[:, ii] = cos.(n .* x);
    end

    for ii = 2:2:N
        n = 1.0 * (ii-odd)
        F[:, ii] = sin.(n .* x)
    end

    F
end

function half_fourier_points(M::Number)
    (0:M-1) .* (2π/(2M-1))
end


# Make an M × N fourier matrix F
function half_fourier_matrix(M::Integer, N::Integer)
    half_fourier_matrix(half_fourier_points(M), N::Integer)
end


# Matrix that takes a function on 2Nθ1+1 nodes and moves it to 2Nθ2+1 nodes
@memoize function collocation_raising_matrix(Nθ1::Integer, Nθ2::Integer)
    F1 = half_fourier_matrix(Nθ1, Nθ1);
    F2 = half_fourier_matrix(Nθ2, Nθ1);

    F2*inv(F1)
end


function *(A::SpatialPowerSeriesSlice{T}, B::SpatialPowerSeriesSlice{T}) where {T}
    Na = size(A.a, 2); Nb = size(B.a, 2); Nc = Na+Nb-1;
    Ra = collocation_raising_matrix(Na, Nc)
    Rb = collocation_raising_matrix(Nb, Nc)

    SpatialPowerSeriesSlice((A.a*Ra') .* (B.a*Rb'))
end

## SpatialPowerSeriesFunctions

function get_N(A::SpatialPowerSeries{T}) where {T}
    length(A.a)
end

function get_M(A::SpatialPowerSeries{T}) where {T}
    get_M(A[1])
end

function get_a(A::SpatialPowerSeries{T}) where {T}
    A.a
end

function get_p0(A::AbstractPowerSeries)
    A.p0[1]
end

function set_p0!(A::AbstractPowerSeries, p0::Integer)
    A.p0[1] = p0;
end


function get_slice(A::SpatialPowerSeries{T}, n::Integer) where {T}
    A.a[n]
end

function getindex(A::SpatialPowerSeries, jj::Integer)
    (jj <= get_N(A)) ? A.a[jj] : ZeroPowerSeriesSlice()
end

function setindex!(A::SpatialPowerSeries, X, ind::Integer)
    set_slice!(A.a[ind], X)
end


### Begin Added 
function *(A::SpatialPowerSeriesSlice{T}, a::T) where {T}
    *(a, A)
end

function times_j!(C::SpatialPowerSeries{T}, A::FluxPowerSeries{T},
                  B::SpatialPowerSeries{T}, jj::Integer) where {T}
#     println("C[jj].a ", size(C[jj].a))
    C[jj] = A[1]*B[jj]
    for kk = 2:jj
        Rb = collocation_raising_matrix(size(B[jj-kk+1].a, 2), size(C[jj].a, 2))
#         println("B[jj-kk+1].a ", B[jj-kk+1].a)
        C[jj] = C[jj] + A[kk]*SpatialPowerSeriesSlice(B[jj-kk+1].a*Rb')
    end
end

function times_j!(C::SpatialPowerSeries{T}, B::SpatialPowerSeries{T}, 
                  A::FluxPowerSeries{T}, jj::Integer) where {T}
    times_j!(C,A,B,jj)
end

# function *(A::FluxPowerSeries{T}, B::SpatialPowerSeries{T};
#            N::Integer=-1) where {T}
#     N = (N==-1) ? min(get_N(A), get_N(B)) : N;
#     M = get_M(A);
# #     println("A's M: ", M)
# #     println("B's M: ", get_M(B))
# #     @assert M == get_M(B)
#     p0 = get_p0(A) + get_p0(B)
#     C = similar(B; N, p0); # C needs to be Spatial Series
#     for jj = 1:N
#         times_j!(C, A, B, jj);
#     end
    
#     C
# end

# function *(A::SpatialPowerSeries{T}, B::FluxPowerSeries{T}; 
#            N::Integer=-1) where {T}
#     *(B, A, N=N)
# end

function +(A::SpatialPowerSeries, B::FluxPowerSeries)
    # This is made because we always want the output series type to be SpatialPowerSeries
    p0 = get_p0(A);
    p0B = get_p0(B);
    if p0 < p0B
        return +(A, distribute_p0(B, p0))
    elseif p0 > p0B
        return +(distribute_p0(A, p0B), B)
    end

    N = min(get_N(A), get_N(B))

    C = similar(A; N)
    for ii = 1:N
        C[ii] = A[ii] + B[ii];
    end

    C
end

function +(A::FluxPowerSeries, B::SpatialPowerSeries)
    +(B, A)
end

function -(A::SpatialPowerSeries, B::FluxPowerSeries)
    p0 = get_p0(A);
    p0B = get_p0(B);
    if p0 < p0B
        return -(A, distribute_p0(B, p0))
    elseif p0 > p0B
        return -(distribute_p0(A, p0B), B)
    end

    N = min(get_N(A), get_N(B))

    C = similar(A; N)
    for ii = 1:N
        C[ii] = A[ii] - B[ii];
    end

    C
end

function -(A::FluxPowerSeries, B::SpatialPowerSeries)
    -(B, A)
end

# This is inefficient. With some more logic, we wouldn't have to allocate B
function +(A::SpatialPowerSeries, a::Number; N::Integer=-1) 
    M = get_M(A);
    N = (N==-1) ? get_N(A) : N;
    
    B = zero_SpatialPowerSeries(M,1)
    B[1].a[:] .= a;

    +(A, B; N)
end

function +(a::Number, A::SpatialPowerSeries; N::Integer=-1)
    +(A, a; N)
end

function -(A::SpatialPowerSeries, a::Number; N::Integer=-1)
    +(A,-a;N)
end

function -(a::Number, A::SpatialPowerSeries; N::Integer=-1)
    +(-1. * A,a;N)
end

### End


"""
    
"""
# function surface_integrate

function zero_SpatialPowerSeries(M::Integer, N::Integer; p0::Integer = 0)
    a = [zero_SpatialPowerSeriesSlice(M, jj) for jj = 1:N];
    SpatialPowerSeries(a; p0)
end

function zero_SpatialPowerSeries(T::DataType, M::Integer, N::Integer; p0::Integer = 0)
    a = [zero_SpatialPowerSeriesSlice(T, M, jj) for jj = 1:N];
    SpatialPowerSeries(a; p0)
end


function similar(A::SpatialPowerSeries{T}; M::Union{Integer, Nothing}=nothing,
        N::Union{Integer, Nothing}=nothing, p0::Union{Integer, Nothing}=nothing) where {T}
#
    M =  isnothing(M)  ? get_M(A)  : M;
    N =  isnothing(N)  ? get_N(A)  : N;
    p0 = isnothing(p0) ? get_p0(A) : p0;

    zero_SpatialPowerSeries(T, M, N; p0)
end


function distribute_p0(A::SpatialPowerSeries{T}, p0::Integer) where {T}
    p0A = get_p0(A);
    @assert p0 <= p0A;
    @assert mod(p0A - p0, 2) == 0

    offset = p0A - p0;
    N = get_N(A) + offset

    B = similar(A; p0, N);
    for ii = offset+1:N
        F = collocation_raising_matrix(ii - offset, ii)
        B[ii] = get_a(A[ii-offset])*F'
    end

    B
end

function distribute_p0(A::SpatialPowerSeries); distribute_p0(A, 0); end


function surface_integrate!(B::FluxPowerSeries{T}, A::SpatialPowerSeries{T}, ii::Integer) where{T}
    @assert mod(ii, 2) == 1
    M = get_M(A);

    # This is inefficient but whatever
    Fθ_inv = inv(half_fourier_matrix(ii, ii))

    B.a[ii] = ((2π)^2 / M) * sum(get_a(A[ii])*Fθ_inv[1,:])
end


function Base.show(io::IO,  ::MIME"text/plain", A::SpatialPowerSeriesSlice{T}; showname=true, NM=11) where {T}
    M, jj = size(A.a);
    if showname; println(io, "SpatialPowerSeriesSlice{$(T)}"); end
    θs = half_fourier_points(jj)
    ss = (0:M-1) .* (2π/M)

    Njj = 7;

    top = @sprintf "n=%-3s     " (@sprintf "%d" jj-1)
    for ii = 1:jj
        if (jj <= Njj) || ((ii < Njj-2) || (ii >= jj-1) && (jj > Njj))
            s = @sprintf "θ=%.3f" θs[ii]
            top = top * (@sprintf "%-11s" s)
        elseif ii == Njj-1
            top = top * "…  "
        end
    end
    printstyled(io, "$top\n"; bold=true)
    # println("$top2")

    for ii = 1:M
        if (M <= NM) || ((ii < NM-1) || (ii >= M-1) && (M > NM))
            line = @sprintf "s=%.3f" ss[ii]
            line = @sprintf "%-10s" line
            printstyled(io, "$line"; bold=true)

            line = ""
            for kk = 1:jj
                if (jj <= Njj) || ((kk < Njj-2) || (kk >= jj-1) && (jj > Njj))
                    line = line * @sprintf "%-11s" (@sprintf "%.3e" A.a[ii, kk])
                elseif kk == Njj-1
                    line = line * "…  "
                end
            end

            println(io, "$line")
        elseif (ii==NM)
            line = @sprintf "%-10s" " :"
            printstyled(io, "$line"; bold=true)

            line = ""
            for kk = 1:jj
                if (jj <= Njj) || ((kk < Njj-2) || (kk >= jj-1) && (jj > Njj))
                    line = line * @sprintf "%-11s" (" :")
                elseif kk == Njj-1
                    line = line * "⋱  "
                end
            end
            println(io, "$line")
        end
    end
end

function Base.show(io::IO, M::MIME"text/plain",A::SpatialPowerSeries{T}) where {T}
    N = get_N(A);
    println(io, "SpatialPowerSeries{$(T)}");
    println(io, "  Order:        N=$(N)")
    println(io, "  No. s-points: M=$(get_M(A))")
    println(io, "  ρ offset:     p0=$(get_p0(A))")
    println(io, "")
    println(io, "Slices:")
    for ii = 1:min(N, 3)
    # for ii = 1:N
        show(io, M, A[ii]; showname=false, NM=7)
        println(io, "")
    end
    if N > 3
        println(io, " ⋮\n")
    end
end