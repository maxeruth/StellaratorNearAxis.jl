"""
    -(A::AbstractPowerSeries)

Get the negative power series of `A` (i.e., `A + (-A) == 0`).
"""
function -(A::AbstractPowerSeries)
    B = similar(A);
    for ii = 1:get_N(B)
        B[ii] = -A[ii]
    end

    B
end

"""
    *(a::Number, A::AbstractPowerSeries)

Multiply the power series `A` by a scalar number.
"""
function *(a::Number, A::AbstractPowerSeries; kwargs...)
    B = similar(A);
    for jj = 1:get_N(B)
        B[jj] = a*A[jj];
    end

    B
end

function *(A::AbstractPowerSeries, a::Number; kwargs...)
    return *(a, A)
end

function /(A::AbstractPowerSeries, a::Number)
    return *(A, 1/a)
end

"""
    +(A::AbstractPowerSeries, B::AbstractPowerSeries)

Add two AbstractPowerSeries. May not work for all inputs (e.g. a SpectralPowerSeries and a SpatialPowerSeries).
"""
function +(A::AbstractPowerSeries, B::AbstractPowerSeries; N::Integer=-1)
    p0 = get_p0(A);
    p0B = get_p0(B);
    if p0 < p0B
        return +(A, distribute_p0(B, p0))
    elseif p0 > p0B
        return +(distribute_p0(A, p0B), B)
    end
    
    N = (N==-1) ? min(get_N(A), get_N(B)) : N;

    C = similar(A; N)
    for ii = 1:N
        C[ii] = A[ii] + B[ii];
    end

    C
end

"""
    -(A::AbstractPowerSeries, B::AbstractPowerSeries)

Subtract two AbstractPowerSeries. May not work for all inputs (e.g. a SpectralPowerSeries and a SpatialPowerSeries).
"""
function -(A::AbstractPowerSeries, B::AbstractPowerSeries; N::Integer=-1)
    p0 = get_p0(A);
    p0B = get_p0(B);
    if p0 < p0B
        return -(A, distribute_p0(B, p0))
    elseif p0 > p0B
        return -(distribute_p0(A, p0B), B)
    end
    
    N = (N==-1) ? min(get_N(A), get_N(B)) : N;

    C = similar(A; N)
    for ii = 1:N
        C[ii] = A[ii] - B[ii];
    end

    C
end

# """
#     times_j!(C::AbstractPowerSeries{T}, A::AbstractPowerSeries{T},
#              B::AbstractPowerSeries{T}, jj::Integer) where {T}

# Set the `jj`th slice of `C` for multiplication (i.e.~`C[jj] = (A*B)[jj]).
# """
function times_j!(C::AbstractPowerSeries{T}, A::AbstractPowerSeries{T},
                  B::AbstractPowerSeries{T}, jj::Integer) where {T}
    #
    C[jj] = A[1]*B[jj]
    for kk = 2:jj
        C[jj] = C[jj] + A[kk]*B[jj-kk+1]
    end
end

"""
    *(A::AbstractPowerSeries{T}, B::AbstractPowerSeries{T};
           N::Integer=-1, careful::Bool=false, careful_tol::Number=1e-12) where {T}

Multiply two power series together to get `C = A*B` to order `N`.
If `careful==true`, the output order is determined by the leading
order by which `A` and `B` are nonzero to tolerance `careful_tol`.
In this case, the argument for `N` is ignored.
"""
function *(A::AbstractPowerSeries{T}, B::AbstractPowerSeries{T};
           N::Integer=-1, careful::Bool=false, careful_tol::Number=1e-12) where {T}
    
    # Get the correct size of the output
    if careful
        if N!=-1
            @warn ("Tried to set `N` and `careful` simultaneously in "*
            "`*(::AbstractPowerSeries,::AbstractPowerSeries)`. "*
            "Performing the careful multiplication.")
        end

        NA = get_N(A)
        NzeroA = 0
        for ii = 1:NA
            if norm(A[ii])<careful_tol
                NzeroA = ii
            else
                break
            end
        end

        NB = get_N(B)
        NzeroB = 0
        for ii = 1:NB
            if norm(B[ii])<careful_tol
                NzeroB = ii
            else
                break
            end
        end

        NA_eff = NA-NzeroA
        NB_eff = NB-NzeroB
        N_eff = min(NA_eff, NB_eff)

        p0_eff = NzeroA+NzeroB

        N = N_eff + p0_eff
    else
        N = (N==-1) ? min(get_N(A), get_N(B)) : N;
    end

    p0 = get_p0(A) + get_p0(B)
    
    # Create the correct type of output
    if typeof(A) <: SpatialPowerSeries
        M = get_M(A);
        if typeof(B) <: SpatialPowerSeries
            @assert M == get_M(B)
        end
        C = similar(A; N, p0);
    else        
        C = similar(B; N, p0)
    end
    
    # Multiply
    for jj = 1:N
        times_j!(C, A, B, jj);
    end

    C
end




"""
    /(A::AbstractPowerSeries{T}, B::AbstractPowerSeries{T};
           N::Integer=-1, careful::Bool=false, careful_tol::Number=1e-12) where {T}

Divide two power series to get `C = A/B` to order `N`. 
Currently just calls inv() and *(), but could potentially be faster with a better method.
"""
function /(A::AbstractPowerSeries{T}, B::AbstractPowerSeries{T};
           N::Integer=-1, careful::Bool=false, careful_tol::Number=1e-12) where {T}
    *(A,inv(B);N=N,careful,careful_tol)
end

# """
#     *(A::AbstractPowerSeries, B::AbstractArray{<:AbstractPowerSeries}; 
#            N::Integer=-1, careful::Bool=false, careful_tol::Number=1e-12)

# Treat multiplying by an AbstractPowerSeries as a scalar over arrays
# """
function *(A::AbstractPowerSeries, B::AbstractArray{<:AbstractPowerSeries}; 
           N::Integer=-1, careful::Bool=false, careful_tol::Number=1e-12)
    [*(A,Bi;N,careful,careful_tol) for Bi in B]
end

function *(B::AbstractArray{<:AbstractPowerSeries}, A::AbstractPowerSeries; 
            N::Integer=-1, careful::Bool=false, careful_tol::Number=1e-12)
    [*(A,Bi;N,careful,careful_tol) for Bi in B]
end


function *(A::AbstractMatrix{<:AbstractPowerSeries}, 
           B::AbstractVector{<:AbstractPowerSeries}; 
           N::Union{Integer,AbstractVector}=-1, 
           careful::Bool=false, careful_tol::Number=1e-12)
    p,q = size(A);
    @assert q==length(B)

    
    C = Vector{AbstractPowerSeries}(undef, p)
    C[:] = [ZeroPowerSeries() for ii = 1:p]

    Ns = ones(Integer, p)
    Ns[:] .= N

    for jj = 1:q
        for ii = 1:p
            C[ii] = +(C[ii], *(A[ii,jj],B[jj];N=Ns[ii],careful,careful_tol);N=Ns[ii])
        end
    end 

    C
end

# """
#     *(A::AbstractMatrix{<:AbstractPowerSeries}, B::AbstractMatrix{<:AbstractPowerSeries}; N::Union{Integer,AbstractMatrix}=-1) where {T,S}
# """
function *(A::AbstractMatrix{<:AbstractPowerSeries}, 
           B::AbstractMatrix{<:AbstractPowerSeries}; 
           N::Union{Integer,AbstractMatrix}=-1, 
           careful::Bool=false, careful_tol::Number=1e-12)
    p,q = size(A);
    q2, r = size(B)
    @assert q==q2

    C = Matrix{AbstractPowerSeries}(undef, p, q)
    C[:] = [ZeroPowerSeries() for ii = 1:p, kk = 1:q]

    Ns = ones(Integer, p, r)
    Ns[:] .= N

    for ii = 1:p
        for kk = 1:r
            for jj = 1:q
                C[ii,kk] = +(C[ii,kk], *(A[ii,jj],B[jj,kk];N=Ns[ii,kk],careful,careful_tol);N=Ns[ii,kk])
            end
        end
    end 

    C
end



# """
#     inv_j!(B::AbstractPowerSeries{T}, A::AbstractPowerSeries{T},
#                 jj::Integer) where {T}

# Get the `jj`th entry of the inverse of `A` (i.e. `B[jj] = (inv(A))[jj]`)
# """
function inv_j!(B::AbstractPowerSeries{T}, A::AbstractPowerSeries{T},
                jj::Integer) where {T}
    if jj == 1
        B[1] = one(T) ./ get_a(A[1])
    else
        B[jj] = B[jj] - A[2]*B[jj-1]
        for kk = 3:jj
            B[jj] = B[jj] - A[kk]*B[jj-kk+1]
        end

        B[jj] =  B[1]*B[jj]
    end
end

"""
    inv(A::AbstractPowerSeries; N::Integer = -1)

Find the multiplicative inverse of `A` to order `N`. Approximately satisfies `A*inv(A)=1`
"""
function inv(A::AbstractPowerSeries; N::Integer = -1)
    p0 = - get_p0(A);
    N = (N == -1) ? get_N(A) : N;
    B = similar(A; p0, N)
    for jj = 1:N
        inv_j!(B, A, jj);
    end

    B
end


# """
#     exp_j!(B::AbstractPowerSeries{T}, A::AbstractPowerSeries{T}, alpha::T,
#                 jj::Integer) where {T}

# Get the `jj`th entry of e^(alpha A)
# """
function exp_j!(B::AbstractPowerSeries{T}, A::AbstractPowerSeries{T}, alpha::T,
                jj::Integer) where {T}
    if jj == 1
        B[1] = exp.(alpha.*get_a(A[1]));
    else
        B[jj] = (alpha/(jj - 1))*A[2]*B[jj-1]
        for kk = 3:jj
            B[jj] = B[jj] + (alpha*(kk-1)/(jj - 1))*A[kk]*B[jj-kk+1]
        end
    end
end

"""
    exp(A::AbstractPowerSeries{T}[, alpha::T]; N::Integer = -1) where {T}

Exponentiate the power series `A` to get `e^(αA)` to order `N`.
"""
function exp(A::AbstractPowerSeries{T}, alpha::T; N::Integer = -1) where {T}
    p0 = get_p0(A);
    if p0>0
        return exp(distribute_p0(A), alpha)
    end
    @assert get_p0(A) == 0

    N = (N == -1) ? get_N(A) : N
    B = similar(A; N)
    for jj = 1:N
        exp_j!(B, A, alpha, jj);
    end

    B
end

function exp(A::AbstractPowerSeries{T}; N::Integer = -1) where {T}
    exp(A, one(T); N)
end

"""
    sinh(A::AbstractPowerSeries{T}[, alpha::T]) where {T}

Find `sinh(αA) = (e^(αA) - e^(-αA))/2`
"""
function sinh(A::AbstractPowerSeries{T}, alpha::T) where {T}
    0.5*(exp(A, alpha) - exp(A, -alpha))
end
function sinh(A::AbstractPowerSeries{T}) where {T}; sinh(A, one(T)); end

"""
    cosh(A::AbstractPowerSeries{T}[, alpha::T]) where {T}

Find `cosh(αA) = (e^(αA) + e^(-αA))/2`.
"""
function cosh(A::AbstractPowerSeries{T}, alpha::T) where {T}
    0.5*(exp(A, alpha) + exp(A, -alpha))
end
function cosh(A::AbstractPowerSeries{T}) where {T}; cosh(A, one(T)); end

# """
#     sincos_j!(S::AbstractPowerSeries{T}, C::AbstractPowerSeries{T},
#               A::AbstractPowerSeries{T}, ω::Number, jj::Integer) where {T}

# Simultaneously computes the `jj`th entry of sin(ωA) and cos(ωA).
# """
function sincos_j!(S::AbstractPowerSeries{T}, C::AbstractPowerSeries{T},
                   A::AbstractPowerSeries{T}, ω::Number, jj::Integer) where {T}
    if jj == 1
        S[1] = sin.(ω * get_a(A[1]))
        C[1] = cos.(ω * get_a(A[1]))
    else
        S[jj] = (ω/(jj - 1))*A[2]*C[jj-1]
        C[jj] = - (ω/(jj - 1))*A[2]*S[jj-1]
        for kk = 3:jj
            S[jj] = S[jj] + (ω*(kk-1)/(jj - 1))*A[kk]*C[jj-kk+1]
            C[jj] = C[jj] - (ω*(kk-1)/(jj - 1))*A[kk]*S[jj-kk+1]
        end
    end
end

"""
    sincos(A::AbstractPowerSeries[, ω::Number])

Find both the sine and cosine of `A`, i.e. `sin(ωA), cos(ωA) = sincos(A[, ω])`
(note: both `sin` and `cos` call `sincos`), so calling this can save time if you
require both)
"""
function sincos(A::AbstractPowerSeries, ω::Number)
    p0 = get_p0(A);
    if p0>0
        return exp(distribute_p0(A), alpha)
    end
    @assert get_p0(A) == 0

    N = get_N(A);
    S = similar(A);
    C = similar(A);

    for jj = 1:N
        sincos_j!(S, C, A, ω, jj)
    end

    S, C
end
function sincos(A::AbstractPowerSeries{T}) where {T}; sincos(A, one(T)); end

"""
    sin(A::AbstractPowerSeries[, ω::Number])

Finds the sine of a power series `sin(ωA)`. Note: if you also require `cos`, it
is twice as efficient to call `sincos`.
"""
function sin(A::AbstractPowerSeries, ω::Number)
    sincos(A, ω)[1]
end
function sin(A::AbstractPowerSeries{T}) where {T}; sin(A, one(T)); end

"""
    cos(A::AbstractPowerSeries[, ω::Number])

Finds the cosine of a power series `cos(ωA)`. Note: if you also require `sin`,
it is twice as efficient to call `sincos`.
"""
function cos(A::AbstractPowerSeries, ω::Number)
    sincos(A, ω)[2]
end
function cos(A::AbstractPowerSeries{T}) where {T}; cos(A, one(T)); end

# """
#     power_j!(B::AbstractPowerSeries{T}, A::AbstractPowerSeries{T}, alpha::Number,
#              jj::Integer) where {T}

# Get the `jj`th entry of A^α, i.e. `B[jj] = (A^α)[jj]`.
# """
function power_j!(B::AbstractPowerSeries{T}, A::AbstractPowerSeries{T},
                  alpha::Number, jj::Integer) where {T}
    if jj == 1
        B[1] = get_a(A[1]).^alpha
    else
        B[jj] = one(T)*alpha*(jj-1)*A[jj]*B[1];
        for kk = 1:jj-2
            B[jj] = B[jj] + one(T)*(alpha*(jj-kk-1) - (kk))*B[kk+1]*A[jj-kk];
        end

        if typeof(B)<:FluxPowerSeries
            B[jj] = B[jj] / ((jj-1) * A[1])
        else
            B[jj] = Diagonal((jj-1) .* vec(get_a(A[1]))) \ get_a(B[jj])
        end
    end
end

"""
    ^(A::AbstractPowerSeries, alpha::Number; N::Integer=-1)

Raise `A` to the power α to `N` terms in the power series, i.e. compute `A^α`.
"""
function ^(A::AbstractPowerSeries, alpha::Number; N::Integer=-1)
    p0 = get_p0(A);
    if !(typeof(alpha) <: Integer)
        @assert p0 == 0
    else
        p0 = p0*alpha
    end
    N = (N == -1) ? get_N(A) : N;
    B = similar(A; p0, N);

    for jj = 1:N
        power_j!(B, A, alpha, jj);
    end

    B
end

# """
#     LinearAlgebra.cross(a::Vector{S}, b::Vector{T}; Ns=[-1,-1,-1],
#             careful::Bool=false, careful_tol::Number=1e-12) where
#             {S <: Union{AbstractPowerSeries}, 
#             T <: Union{AbstractPowerSeries, Number}}

# Performs the cross product between two 3-vectors of power series.
# """
function LinearAlgebra.cross(a::Vector{S}, b::Vector{T}; Ns=[-1,-1,-1],
                             careful::Bool=false, careful_tol::Number=1e-12) where
             {S <: AbstractPowerSeries, 
              T <: Union{AbstractPowerSeries, Number}}
    @assert length(a) == 3
    @assert length(b) == 3
    c = Vector{AbstractPowerSeries}(undef, 3)
    N1,N2,N3 = Ns
    c[1] = *(a[2],b[3];N=N1,careful,careful_tol)-*(b[2],a[3];N=N1,careful,careful_tol)     
    c[2] = *(a[3],b[1];N=N2,careful,careful_tol)-*(b[3],a[1];N=N2,careful,careful_tol)
    c[3] = *(a[1],b[2];N=N3,careful,careful_tol)-*(b[1],a[2];N=N3,careful,careful_tol)
    
    c
end

function LinearAlgebra.cross(a::Vector{S}, b::Vector{T}; Ns=[-1,-1,-1],
    careful::Bool=false, careful_tol::Number=1e-12) where
            {S <: Number, 
             T <: AbstractPowerSeries}
    cross(b,a;Ns,careful,careful_tol)
end


# """
#     dot(a::Vector{<:AbstractPowerSeries}, b::Vector{<:AbstractPowerSeries})

# Performs the dot product between two 3-vectors of power series.
# """
function LinearAlgebra.dot(a::Vector{<:AbstractPowerSeries}, b::Vector{<:AbstractPowerSeries}; 
                           N::Integer=-1, careful::Bool=false, careful_tol::Number=1e-12)
    @assert length(a) == length(b)
    
    C = *(a[1], b[1]; N, careful, careful_tol)
    for ii = 2:length(a)
        C = C + *(a[ii], b[ii]; N, careful, careful_tol)
    end

    C 
end

# """
#     adjoint(A::AbstractPowerSeries{T}) where {T<:Complex}

# Finds the adjoint of a complex power series.
# """
function adjoint(A::AbstractPowerSeries{T}) where {T<:Complex}
    B = similar(A);
    for ii = 1:get_N(A)
        B[ii] = conj.(get_a(A[ii]));
    end

    B
end

# """
#     adjoint(A::AbstractPowerSeries{T}) where {T<:Complex}

# Finds the adjoint of a real power series.
# """
function adjoint(A::AbstractPowerSeries{T}) where {T<:Real}
    A
end


"""
    LinearAlgebra.norm(A::AbstractVector{<:AbstractPowerSeries})

Returns the norm of a power series `A`, i.e. `(A'*A)^0.5`.
For a norm of the elements of `A`, see [`Fnorm`](@ref).
"""
function LinearAlgebra.norm(A::AbstractVector{<:AbstractPowerSeries})
    (dot(A,A))^0.5
end

# """
#     *(a::S, b::AbstractVector) where {S <: AbstractPowerSeries}

# Multiplies each entry of a vector by a power series.
# """
function *(a::S, b::AbstractVector) where {S <: AbstractPowerSeries}
    [a*bi for bi in b]
end

"""
    surface_integrate(A::AbstractPowerSeries)

Performs the surface integral ∫A dθ dϕ.
"""
function surface_integrate(A::AbstractPowerSeries{T}) where {T}
    N = get_N(A);

    p0 = get_p0(A);
    B = FluxPowerSeries(zeros(N); p0) 

    for ii = 1:2:N
        surface_integrate!(B, A, ii)
    end

    B
end

# """
#     rho_deriv_j!(dA::AbstractPowerSeries{T}, A::AbstractPowerSeries{T},
#                jj::Integer) where {T}

# Finds the `jj`th entry of the derivative of `A` with respect to `ρ`, i.e.
# `dA[jj] = (dA/dρ)[jj]`.
# """
function rho_deriv_j!(dA::AbstractPowerSeries{T}, A::AbstractPowerSeries{T},
                    jj::Integer) where {T}
    p0 = get_p0(A);
    dA[jj] = (jj - one(T) + p0) * A[jj]
end

"""
    rho_deriv(A::AbstractPowerSeries{T}) where {T}

Finds the derivative of `A` with respect to `ρ`.
"""
function rho_deriv(A::AbstractPowerSeries{T}) where {T}
    # println("Untested!!!")
    p0 = get_p0(A)-1;
    dA = similar(A; p0)
    for jj = 1:get_N(A)
        rho_deriv_j!(dA, A, jj)
    end

    dA
end

# """
#     rho_integrate_j!(intA::AbstractPowerSeries{T}, A::AbstractPowerSeries{T},
#                         jj::Integer) where {T}
# """
function rho_integrate_j!(intA::AbstractPowerSeries{T}, A::AbstractPowerSeries{T},
                        jj::Integer) where {T}
    p0 = get_p0(A);
    # println("jj = $jj, one(T)/(jj + 1 + p0) = $(one(T)/(jj + 1 + p0)), A[jj] = $(A[jj])")
    intA[jj] = (one(T)/(jj + 1 + p0)) * A[jj]
end

"""
    rho_integrate(A::AbstractPowerSeries)

Perform the integral
    `B(r) = ∫₀ʳ A(ρ) ρ dρ`
See [`rho_antideriv`](@ref) for the integral without the integrating factor.
"""
function rho_integrate(A::AbstractPowerSeries)
    # println("Untested!!!")
    p0 = get_p0(A)+2;
    intA = similar(A; p0)
    for jj = 1:get_N(A)
        rho_integrate_j!(intA, A, jj)
    end

    intA
end

"""
    rho_antideriv(A::AbstractPowerSeries)

Perform the integral
    `B(r) = ∫₀ʳ A(ρ) dρ`
(note: unlike [`rho_integrate`](@ref), which includes an integrating factor for analytic integrals)
"""
function rho_antideriv(A::AbstractPowerSeries)
    # println("rho_antideriv still untested! Use at your own risk")
    p0 = get_p0(A)+1;
    intA = similar(A; p0)
    for jj = 1:get_N(A)
        rho_antideriv_j!(intA, A, jj)
    end

    intA
end


"""
    volume_integate(A::AbstractPowerSeries{T}; N::Integer=-1) where {T}

Integral over a volume. Implemented by composing [`rho_integrate`](@ref) with 
[`surface_integrate`](@ref).
"""
function volume_integrate(A::AbstractPowerSeries{T}; N::Integer=-1) where {T}
    return rho_integrate(surface_integrate(A;N))
end


"""
    Fnorm(A::AbstractPowerSeries{T}) where {T}

Frobenius norm of the elements of `A`. Useful for determining approximate equality of 
AbstractPowerSeries.
"""
function Fnorm(A::AbstractPowerSeries{T}) where {T}
    my_norm = zero(T);
    for ii = 1:get_N(A);
        my_norm = my_norm + LinearAlgebra.norm(vec(get_a(A[ii])))^2
    end

    sqrt(my_norm);
end

"""
    distribute_p0(A::AbstractPowerSeries, p0::Integer)

If `A` has an order `q0 ≥ p0`, creates a new AbstractPowerSeries of the same type with order offset
`p0`` with equal coefficients. For analyticity to be preserved, we additionally require that `q0` 
and `p0` are even offsets of each other (`p0-q0 = 0 mod 2`). If `p0` is not specified, the default
value is `p0=0`.
"""
function distribute_p0(A::AbstractPowerSeries, p0::Integer)
    p0A = get_p0(A);
    @assert p0 <= p0A;
    @assert mod(p0A - p0, 2) == 0

    offset = p0A - p0;
    N = get_N(A) + offset

    B = similar(A; p0, N);
    for ii = offset+1:N
        distribute_p0!(B, A, ii, offset)
    end

    B
end

function distribute_p0(A::AbstractPowerSeries); distribute_p0(A, 0); end


"""
    unsafe_distribute_p0(A::AbstractPowerSeries[, p0::Integer])

Creates a new PowerSeries with order offset p0 with equal coefficients to `A`. For analyticity to be 
preserved, we require that `q0` and `p0` are even offsets of each other (`p0-q0 = 0 mod 2`). If `p0` 
is not specified, the default value is `p0=0`.

If you are unsure about whether p0 can be increased without losing information, 
[`remove_zeros`](@ref) can automatically check whether the order can be changed.
"""
function unsafe_distribute_p0(A::AbstractPowerSeries, p0::Integer)
    p0A = get_p0(A);
    offset = p0A - p0;
    @assert mod(offset, 2) == 0

    if offset >= 0
        return distribute_p0(A, p0)
    end

    N = get_N(A) + offset

    B = similar(A; p0, N);
    for ii = 1:N
        unsafe_distribute_p0!(B, A, ii, offset)
    end

    B
end

"""
    remove_zeros(A::AbstractPowerSeries, tol=1e-10)

If one has access to a power series `A = ρ^2N B`, returns the power series `B`. For example, when 
the polar Laplacian is applied to a series, the output has zeros in the first two orders. This 
removes those zeros.
"""
function remove_zeros(A::AbstractPowerSeries; tol=1e-10)
    p0 = get_p0(A)
    Nr = get_N(A) 

    idx = 2
    if idx > Nr
        return A
    end

    tmp = unsafe_distribute_p0(A, idx+p0)
    # display(Fnorm(A-tmp))
    while idx < Nr && Fnorm(A-tmp)<tol
        
        A = tmp
        tmp = unsafe_distribute_p0(A, idx+p0)
        idx += 2
        # display(Fnorm(A-tmp))
    end

    return A
end

""" 
    L2norm(A::AbstractPowerSeries, δ::Number)

Find the L2 norm of a power series via near-axis operations.
NOTE: this is not the true L2 norm; just the near-axis approximation of it.
"""
function L2norm(A::AbstractPowerSeries, δ::Number)
    sqrt(evaluate(volume_integrate(A'*A),δ))
end
