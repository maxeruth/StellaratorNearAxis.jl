## SpectralPowerSeriesSlice functions
function get_a(A::SpectralPowerSeriesSlice)
    A.a
end

function zero_SpectralPowerSeriesSlice(M::Integer, jj::Integer)
    SpectralPowerSeriesSlice(zeros(M, jj))
end

function zero_SpectralPowerSeriesSlice(T::DataType, M::Integer, jj::Integer)
    SpectralPowerSeriesSlice(zeros(T, M, jj))
end

function *(a::T, A::SpectralPowerSeriesSlice{T}) where {T}
    SpectralPowerSeriesSlice(a * A.a);
end

function +(A::SpectralPowerSeriesSlice{T}, B::SpectralPowerSeriesSlice{T}) where {T}
    SpectralPowerSeriesSlice(A.a + B.a)
end

function -(A::SpectralPowerSeriesSlice{T}, B::SpectralPowerSeriesSlice{T}) where {T}
    SpectralPowerSeriesSlice(A.a - B.a)
end


## SpectralPowerSeries functions
function get_slice(A::SpectralPowerSeries{T}, n::Integer) where {T}
    A.a[n]
end

function getindex(A::SpectralPowerSeries, jj::Integer)
    (jj <= get_N(A)) ? A.a[jj] : ZeroPowerSeriesSlice()
end

function setindex!(A::SpectralPowerSeries, X, ind::Integer)
    set_slice!(A.a[ind], X)
end

""" zero_SpectralPowerSeries(M::Integer, N::Integer; p0::Integer = 0)

Return a SpectralPowerSeries object initialized with zero coefficients.

Inputs:
 - `M`: Number of Fourier coefficients in the poloidal s coordinate
 - `N`: Number of coefficients in radial ρ and toroidal θ coordinates
 - `p0`: The leading order ρ coefficient

Output:
 - A zero SpectralPowerSeries
"""
function zero_SpectralPowerSeries(M::Integer, N::Integer; p0::Integer = 0)
    a = [zero_SpectralPowerSeriesSlice(M, jj) for jj = 1:N];
    SpectralPowerSeries(a; p0)
end

function zero_SpectralPowerSeries(T::DataType, M::Integer, N::Integer; p0::Integer = 0)
    a = [zero_SpectralPowerSeriesSlice(T, M, jj) for jj = 1:N];
    SpectralPowerSeries(a; p0)
end


function get_N(A::SpectralPowerSeries{T}) where {T}
    length(A.a)
end

function get_M(A::SpectralPowerSeries{T}) where {T}
    get_M(A.a[1])
end

function Base.show(io::IO,  ::MIME"text/plain", A::SpectralPowerSeriesSlice{T}; showname=true, NM=11) where {T}
    M, jj = size(A.a);
    if showname; println(io, "SpectralPowerSeriesSlice{$(T)}"); end
    m=0
    if iseven(jj)
        m=1
    end

    Njj = 7;

    top = @sprintf "n=%-3s     " (@sprintf "%d" jj-1)
    for ii = 1:jj
        if (jj <= Njj) || ((ii < Njj-2) || (ii >= jj-1) && (jj > Njj))
            m = isodd(jj) ? 2*(ii÷2) : 2*((ii-1)÷2) + 1
            s = ""
            if iseven(ii)
                s = (m==1) ? "sin(θ)" : @sprintf "sin(%dθ)" m
            else
                s = m==0 ? "1" : ((m==1) ? "cos(θ)" : @sprintf "cos(%dθ)" m)
            end
            top = top * (@sprintf "%-11s" s)
        elseif ii == Njj-1
            top = top * "…  "
        end
    end
    printstyled(io, "$top\n"; bold=true)
    # println("$top2")

    for ii = 1:M
        if (M <= NM) || ((ii < NM-1) || (ii >= M-1) && (M > NM))
            line = ""
            m = ii÷2
            if iseven(ii)
                line = (m==1) ? "sin(s)" : @sprintf "sin(%ds)" m
            else
                line = m==0 ? "1" : ((m==1) ? "cos(s)" : @sprintf "cos(%ds)" m)
            end
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

function Base.show(io::IO, M::MIME"text/plain",A::SpectralPowerSeries{T}) where {T}
    N = get_N(A);
    println(io, "SpectralPowerSeries{$(T)}");
    println(io, "  Order:       N=$(N)")
    println(io, "  No. s-modes: M=$(get_M(A))")
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


function similar(A::SpectralPowerSeries{T}; M::Union{Integer, Nothing}=nothing,
        N::Union{Integer, Nothing}=nothing, p0::Union{Integer, Nothing}=nothing) where {T}
#
    M =  isnothing(M)  ? get_M(A)  : M;
    N =  isnothing(N)  ? get_N(A)  : N;
    p0 = isnothing(p0) ? get_p0(A) : p0;

    zero_SpectralPowerSeries(T, M, N; p0)
end

# # Get the composed basis
# # basis[2n:2n+1] = [sin(n*(twist*ϕ + B)), cos(n*(twist*ϕ + B))]
# function composition_basis(A::SpatialPowerSeries{T}, M::Integer;
#                            twist::Integer=1) where {T}
#     # println("Untested!!!")
#     N = get_N(A);
#     MA = get_M(A);
#     p0 = get_p0(A);
#     @assert p0 == 0;
#     # basis = [zero_SpatialPowerSeries(MA, N) for ii = 1:N]
#     basis = Vector{SpatialPowerSeries{T}}(undef, M);

#     basis[1] = zero_SpatialPowerSeries(T, MA, N; p0);
#     basis[1][1] = ones(T, MA);

#     Φ = zero_SpatialPowerSeries(T, MA, N; p0);
#     Φ[1] = (0:MA-1).*(twist*2π/MA);
#     twisted_A = A+Φ;

#     for m = 2:2:M
#         S, C = sincos(twisted_A, m/2.);
#         basis[m] = S;
#         if m < M
#             basis[m+1] = C;
#         end
#     end

#     basis
# end

# ## Function to compose with Fourier coefficients stored in A
# # f(ϕ) = A[0] + A[1] sin(ϕ) + A[2] cos(ϕ) + A[3] sin(2ϕ) + ...
# function ϕ_compose(A::Vector{T}, basis::Vector{SpatialPowerSeries{T}}) where {T}
#     #
#     MA = length(A);
#     @assert 1 <= MA <= length(basis);

#     C = A[1]*basis[1];
#     for ii = 2:MA
#         C = C + A[ii]*basis[ii]
#     end

#     C
# end



# Make an M × N fourier matrix F
function full_fourier_matrix(x::AbstractVector, N::Integer)
    M = length(x)
    F = zeros(M, N)

    # 
    Ncos = (N+1)÷2
    for n = 0:Ncos-1
        @avx for jj = 1:M
            @inbounds F[jj, 2n+1] = cos(n * x[jj]);
        end
    end

    #
    Nsin = N÷2
    #
    for n = 1:Nsin
        @avx for jj = 1:M
            F[jj, 2n] = sin(n * x[jj])
        end
    end

    F
end

# Make an M × N fourier matrix F
function full_fourier_matrix(M::Integer, N::Integer)
    x = (0:M-1) .* (2π/M)
    full_fourier_matrix(x, N)
end



function fourier_weight_matrix(MA::Integer, M::Integer)
    w = ones(M) .* MA/2;
    w[1] = MA;
    Diagonal(w);
end

function SpectralPowerSeries_j!(B::SpectralPowerSeries{T}, A::SpatialPowerSeries{T},
                                jj::Integer, Fϕ::AbstractArray{T}, W::AbstractArray{T}) where T
    Fθ = half_fourier_matrix(jj, jj)
    B[jj] = W \ (Fϕ'*get_a(A[jj])*inv(Fθ'))
end

"""
    SpectralPowerSeries(A::SpatialPowerSeries{T}; M::Integer=-1) where {T}

Transform a SpatialPowerSeries to a SpectralPowerSeries of the same order.
By default, the output uses the same number of Fourier modes as the number of spatial collocation nodes.
Alternatively, the user can set `M` to choose a different number of modes.
"""
function SpectralPowerSeries(A::SpatialPowerSeries{T}; M::Integer=-1) where {T}
    N = get_N(A);
    MA = get_M(A)
    M = (M == -1) ? MA : M;
    @assert mod(M, 2) == 1

    p0 = get_p0(A);

    B = zero_SpectralPowerSeries(T, M, N; p0)

    Fϕ = full_fourier_matrix(MA, M);
    W = fourier_weight_matrix(MA, M);

    for jj = 1:N
        SpectralPowerSeries_j!(B, A, jj, Fϕ, W)
    end
    B
end

function SpectralPowerSeries(::ZeroPowerSeries{T}; M::Union{Integer, Nothing}=nothing,
        N::Union{Integer, Nothing}=nothing, p0::Union{Integer, Nothing}=nothing) where {T}
    N =  isnothing(N)  ? 1 : N;
    M =  isnothing(M)  ? 1 : M;
    p0 = isnothing(p0) ? 0 : p0;

    zero_SpectralPowerSeries(T, M, N; p0)
end


function SpatialPowerSeries_j!(B::SpatialPowerSeries{T}, A::SpectralPowerSeries{T},
                                jj::Integer, Fϕ::AbstractMatrix) where T
    Fθ = half_fourier_matrix(jj, jj)
    B[jj] = Fϕ*get_a(A[jj])*Fθ'
end

"""
    SpatialPowerSeries(A::SpectralPowerSeries{T}; M::Integer = -1) where {T}

Transform a SpectralPowerSeries to a SpatialPowerSeries of the same order.
By default, the output uses the same number of collocation nodes as the Fourier modes.
Alternatively, the user can set `M` to choose a different number of nodes.
"""
function SpatialPowerSeries(A::SpectralPowerSeries{T}; M::Integer = -1) where {T}
    NA = get_N(A)
    N = NA

    MA = get_M(A);
    M = (M == -1) ? MA : M;

    p0 = get_p0(A);

    B = zero_SpatialPowerSeries(T, M, N; p0)

    Fϕ = full_fourier_matrix(M, MA);

    for jj = 1:min(N,NA)
        SpatialPowerSeries_j!(B, A, jj, Fϕ)
    end
    B
end

function s_deriv_operator(M)
    A = zeros(M,M);

    for m = 2:2:M
        k = m÷2
        A[m, m+1] = -k;
        A[m+1, m] =  k;
    end

    A
end

function s_deriv_j!(dA::SpectralPowerSeries{T}, A::SpectralPowerSeries{T}, jj::Integer) where {T}
    M = get_M(dA);
    dAslice = get_a(dA[jj]);
    Aslice = get_a(A[jj]);

    for m = 2:2:M
        k = m÷2
        dAslice[m,   :] = -k*Aslice[m+1, :];
        dAslice[m+1, :] =  k*Aslice[m,   :];
    end
end

"""
    s_deriv(A::SpectralPowerSeries{T}) where {T}

Takes the derivative of `A` with respect to `s`.

WARNING: currently this function is only defined for an odd number of Fourier modes.
"""
function s_deriv(A::SpectralPowerSeries{T}) where {T}
    @assert mod(get_M(A), 2) == 1 # Idk how to handle the even case
    N = get_N(A);
    dA = similar(A);
    for ii = 1:N
        s_deriv_j!(dA, A, ii)
    end
    dA
end

function theta_deriv_operator(n)
    A = zeros(n,n);
    
    odd = mod(n, 2);
    s = (odd==1 ? -1 : 1)

    for kk = 1:n÷2
        m = (1 + odd + 2*(kk-1))

        A[m, m+1] =  s*m;
        A[m+1, m] = -s*m;
    end

    A
end

function theta_deriv_j!(dA::SpectralPowerSeries{T}, A::SpectralPowerSeries{T}, jj::Integer) where {T}
    dAslice = dA[jj];
    Aslice = A[jj];

    odd = mod(jj , 2);
    if odd==1
        dAslice.a[:, 1] .= zero(T);
    end
    for k = 1:jj÷2
        m = (1 + odd + 2*(k-1))
        # Used to make sure that the minus sign of the derivative is in the right place
        s = (odd==1 ? -1 : 1)
        dAslice.a[:, m] = s*m*Aslice.a[:, m+1];
        dAslice.a[:, m+1] = -s*m*Aslice.a[:, m];
    end

end

"""
    theta_deriv(A::SpectralPowerSeries{T}) where {T}

Takes the derivative of `A` in the `θ` direction.
"""
function theta_deriv(A::SpectralPowerSeries{T}) where {T}
    N = get_N(A);
    dA = similar(A);
    for jj = 1:N
        theta_deriv_j!(dA, A, jj);
    end

    dA
end


"""
    grad(A::SpectralPowerSeries{T}) where {T}

Returns the "gradient" of `[dA/dρ, dA/dθ, dA/ds]`. 
Note that for a true gradient, it is necessary to use the metric.
"""
function grad(A::SpectralPowerSeries{T}) where {T}
    [rho_deriv(A), theta_deriv(A), s_deriv(A)]
end


"""
    div(A::SpectralPowerSeries{T}) where {T}

Returns the "divergence" `dA[1]/dρ + dA[2]/dθ + dA[3]/ds`. 
Note that for a true divergence, it is necessary to use the metric.
"""
function div(A::Vector{SpectralPowerSeries{T}}) where {T}
    rho_deriv(A[1]) + theta_deriv(A[2]) + s_deriv(A[3])
end

function distribute_p0!(B::SpectralPowerSeries, A::SpectralPowerSeries, 
                        ii::Integer, offset::Integer)
    B[ii].a[:, 1:ii-offset] = get_a(A[ii-offset])
end

function unsafe_distribute_p0!(B::SpectralPowerSeries, A::SpectralPowerSeries, 
                               ii::Integer, offset::Integer)
    B[ii] = get_a(A[ii-offset])[:, 1:ii]
end



function unsafe_distribute_p0(A::SpectralPowerSeries); unsafe_distribute_p0(A, 0); end

function remove_zeros(series::SpatialPowerSeries; tol=1e-10)
    return SpatialPowerSeries(remove_zeros(SpectralPowerSeries(series); tol))
end


# """
#     evaluate_j(A::SpectralPowerSeries, ρ::AbstractVector,
#                     θ::AbstractVector, ϕ::AbstractVector, jj::Integer)

# Evaluate the SpectralPowerSeriesSlice either on a grid.
# """
function evaluate_j(A::SpectralPowerSeries, ρ::AbstractVector,
                    θ::AbstractVector, Φ::AbstractArray, jj::Integer)
    Nρ = length(ρ); Nθ = length(θ); Nϕ = size(Φ,1);

    ajj = A[jj].a;
    Nr = size(ajj, 2);

    R = ρ.^(jj-1+get_p0(A));
    Θ = half_fourier_matrix(θ, Nr);

    reshape( R * reshape( Θ*ajj'*Φ' ,1,Nθ*Nϕ),Nρ,Nθ,Nϕ)
end

# This is probably typed poorly. We should probably be more careful here
"""
    evaluate(A::SpectralPowerSeries, rho::Union{Number, AbstractVector},
             theta::Union{Number, AbstractVector}, phi::Union{Number, AbstractVector})

Evaluate the SpectralPowerSeries either at a point or on a grid.
"""
function evaluate(A::SpectralPowerSeries, rho::Union{Number, AbstractVector},
                  theta::Union{Number, AbstractVector},
                  phi::Union{Number, AbstractVector})
    if typeof(rho) <: Number; return evaluate(A, [rho], theta, phi); end
    if typeof(theta) <: Number; return evaluate(A, rho, [theta], phi); end
    if typeof(phi) <: Number; return evaluate(A, rho, theta, [phi]); end

    Nrho = length(rho); Ntheta = length(theta); Nphi = length(phi);
    M = get_M(A);

    a = zeros(Nrho, Ntheta, Nphi)
    Phi = full_fourier_matrix(phi, M)
    for jj = 1:get_N(A)
        a[:,:,:] = a + evaluate_j(A, rho, theta, Phi, jj)
    end

    a
end

function section_j!(B::SpectralPowerSeries, A::SpectralPowerSeries, Φ::AbstractArray, jj::Integer)
    B[jj] = Φ*A[jj].a
end

# """
#     section(A::SpectralPowerSeries, s::Number)

# Find a polar expression for the SpectralPowerSeries `A` at a fixed value of `s`
# """
function section(A::SpectralPowerSeries, ϕ::Number)
    N = get_N(A);
    B = zero_SpectralPowerSeries(1,N);
    Φ = full_fourier_matrix([ϕ], get_M(A))

    for jj = 1:N
        section_j!(B, A, Φ, jj)
    end
    set_p0!(B, get_p0(A))

    B
end

"""
    composition_basis(F::AbstractVector{<:SpectralPowerSeries})

Get a basis of the form
-`[[1], [ρF cos(θF), ρF sin(θF)], [ρF^2, ρF^2 cos(2 θF), ρF^2 sin(2 θF)], ...]`
from the input of a SpectralPowerSeries
-`F = [ρF cos(θF), ρF sin(θF)]`
This can be used in conjunction with the function `compose` to find
the composition of another SpectralPowerSeries with `F`.
"""
function composition_basis(F::AbstractVector{<:SpectralPowerSeries}; M::Integer=-1)
    @assert length(F) == 2
    MF = get_M(F[1]);
    @assert get_M(F[2]) == MF;
    N = get_N(F[1])
    @assert get_N(F[2]) == N
    @assert iszero(norm(F[1][1].a))
    @assert iszero(norm(F[2][1].a))

    M = (M==-1) ? MF : M

    Ftree = [[zero_SpatialPowerSeries(M,N) for n = 0:m] for m = 0:N-1];
    Ftree[1][1][1].a[:] .= 1.0; # The leading order behavior is constant
    Ftree[2][:] .= SpatialPowerSeries.(F; M);

    c = Ftree[2][1];
    s = Ftree[2][2];

    for m=3:N
        if isodd(m)
            cm1 = Ftree[m-1][1]
            sm1 = Ftree[m-1][2]
            Ftree[m][1] = c*cm1 + s*sm1 # 1 = cos^2 + sin^2

            for n = 1:m÷2
                cm1 = Ftree[m-1][2n-1]
                sm1 = Ftree[m-1][2n];

                Ftree[m][2n]   = sm1*c + cm1*s # sin((n+1)θ) = sin(nθ)cos(θ) + cos(nθ)sin(θ)
                Ftree[m][2n+1] = cm1*c - sm1*s # cos((n+1)θ) = cos(nθ)cos(θ) - sin(nθ)sin(θ)
            end
        else
            cm1 = Ftree[m-1][1];
            Ftree[m][1] = c*cm1;
            Ftree[m][2] = s*cm1;

            for n = 1:(m÷2)-1
                sm1 = Ftree[m-1][2n]
                cm1 = Ftree[m-1][2n+1]

                Ftree[m][2n+1] = cm1*c - sm1*s # cos((n+1)θ) = cos(nθ)cos(θ) - sin(nθ)sin(θ)
                Ftree[m][2n+2] = sm1*c + cm1*s # sin((n+1)θ) = sin(nθ)cos(θ) + cos(nθ)sin(θ)
            end
        end
    end


    # [SpectralPowerSeries.(Ftree_i) for Ftree_i in Ftree]
    Ftree
end

"""
    compose(A::SpectralPowerSeries, Ftree::AbstractVector)

Compose a SpectralPowerSeries `A` with the basis `Ftree`. See
[`composition_basis`](@ref) for a method to obtain `Ftree` from a
change-of-coordinates SpectralPowerSeries `F`.
"""
function compose(A::SpectralPowerSeries, Ftree::AbstractVector)
    N = get_N(A)
    Ms = get_M(A);
    Mc = get_M(Ftree[1][1])
    B = zero_SpatialPowerSeries(Mc,N)
    F = full_fourier_matrix(Mc, Ms); 
    
    Aij = zero_SpatialPowerSeries(Mc,1)
    for ii = 1:N
        aii = F*A[ii].a
        for jj = 1:ii
            Aij[1] = aii[:,jj]
            B = B + *(Aij, Ftree[ii][jj]; N=N)
        end
    end

    SpectralPowerSeries(B, M=Ms);
end

"""
    invert_coordinates(F::AbstractVector)

Invert the SpectralPowerSeries coordinate transform `F`. That is, 
the composition via [`compose`](@ref) of `F` with the output `G` 
satisfies `F∘G(x) ≈ x`
"""
function invert_coordinates(F::AbstractVector; M::Integer=-1)
    Ms = get_M(F[1]);
    Mc = (M==-1) ? Ms : M
    Nρ = get_N(F[1]);

    Ftree = composition_basis(F; M=Mc)
    G_s = [zero_SpectralPowerSeries(Ms,Nρ) for ii = 1:2]
    G_c = [zero_SpatialPowerSeries(Mc,Nρ) for ii = 1:2]

    # Get the leading order behavior
    Fθ2 = half_fourier_matrix(2, 2)

    Fs = zeros(Mc,2,2)
    for ii = 1:2
        Fs[:,ii,:] = Ftree[2][ii][2].a
    end
    Fs = reshape(reshape(Fs,2Mc,2) * inv(Fθ2'),Mc,2,2)

    G0s = similar(Fs);
    for ii = 1:Mc
        G0s[ii,:,:] = inv(Fs[ii,:,:])
    end
    G0s = reshape(reshape(G0s,2Mc,2)*Fθ2',Mc,2,2)

    for ii = 1:2
        G_c[ii][2] = G0s[:,ii,:]
        G_s[ii] = SpectralPowerSeries(G_c[ii], M=Ms)
    end

    # Get the nth order stuff
    # TODO: Should be able to work out a way of doing this that turns it into a problem
    # of perturbing around the identity, removing the $A$ step
    for n = 3:Nρ
    # for n = 3:3
        Fθ = half_fourier_matrix(n,n);
        res = [compose(G_si, Ftree) for G_si in G_s]
        res[1][2].a[1,1] = res[1][2].a[1,1] - 1
        res[2][2].a[1,2] = res[2][2].a[1,2] - 1
        resc = [SpatialPowerSeries(resi, M=Mc) for resi in res]

        R = zeros(Mc,2,n)
        for ii = 1:2
            R[:,ii,:] = resc[ii][n].a #* inv(Fθ')
        end

        As = zeros(Mc, n, n)
        for ii = 1:n
            As[:, ii, :] = Ftree[n][ii][n].a
        end

        Gns = zeros(Mc, 2, n)
        for ii = 1:Mc
            Gns[ii,:,:] = -R[ii,:,:] * inv(As[ii,:,:]) * Fθ'
        end

        for ii = 1:2
            G_c[ii][n] = Gns[:,ii,:]
            G_s[ii] = SpectralPowerSeries(G_c[ii], M=Ms)
        end
        # [println("component $ii, slice $jj error = $(norm(res[ii][jj].a))") for ii=1:2, jj = 1:Nρ]
    end

    G_s
end

"""
    fit_SpectralPowerSeries(a::AbstractArray, r::AbstractVector, M::Integer, N::Integer)

Fit a SpectralPowerSeries to the samples of some function.

input:
- `a`: A Nr × Nθ × Nϕ array of function samples. It is assumed that
   `a[i,j,k] = f(r[i], θ[j], ϕ[k])`, where `r` is given as input, and `θ` and
   `ϕ` are on a uniform grid from `0` to `2π`
- `r`: The radial coordinates sampled at
- `N`: The order in ρ and θ the series is interpolated to
- `M`: The order in ϕ the series is interpolated to
"""
function fit_SpectralPowerSeries(a::AbstractArray, r::AbstractVector,
                                 N::Integer, M::Integer)

#
    Nr, Nθ, Nϕ = size(a)
    @assert (Nr ≥ (any(r .== 0.) ? N+1 : N) ) && (Nϕ≥M)

    θ = (0:Nθ-1).*(2π/Nθ)

    # We will perform the fit separately in the toroidal direction and
    # in the cross section (hopefully this actually solves the least-squares
    # problem)
    ϕ_modes = full_fourier_matrix(Nϕ, M)
    # Nmax = maximum([Nr, Nθ])
    Nmax = N
    rθ_modes = zeros(Nr*Nθ, Nmax*(Nmax+1)÷2)
    k = 1
    for ii = 1:Nmax
        rθ_modes[:, k:k+ii-1] = kron(half_fourier_matrix(θ, ii), r.^(ii-1))
        k = k+ii
    end

    Wϕ = fourier_weight_matrix(Nϕ, M)
    c1 = reshape(a, Nr*Nθ, Nϕ)*ϕ_modes*inv(Wϕ)
    coeffs = rθ_modes\c1

#     coeffs = rθ_modes \ (reshape(a, Nr*Nθ, Nϕ)*ϕ_modes*inv(Wϕ))

    ### Begin Testing
    # ϕ_resid = LinearAlgebra.norm(reshape(a, Nr*Nθ, Nϕ) - c1*ϕ_modes')/LinearAlgebra.norm(reshape(a, Nr*Nθ, Nϕ))
    # println("ϕ resid: ", ϕ_resid)
    # rθ_resid = LinearAlgebra.norm(c1 - rθ_modes*coeffs)
    # println("rθ resid: ", rθ_resid)
    ### End Testing

    A = zero_SpectralPowerSeries(M, N; p0 = 0)
    k = 1
    for ii = 1:N
        A[ii] = coeffs[k:k+ii-1, :]'
        k = k+ii
    end

    A
end



"""
    fit_SpectralPowerSeries(a::AbstractArray, r::AbstractVector, M::Integer, N::Integer)

Fit a SpectralPowerSeries to the samples of some function.

input:
- `f`: A function from R^3→R, with signature `f(r,θ,ϕ)`
- `r`: The radial coordinates sampled at
- `Nθ`: Number of θ samples
- `Nϕ`: Number of ϕ samples
- `N`: The order in ρ and θ the series is interpolated to
- `M`: The order in ϕ the series is interpolated to
"""
function fit_SpectralPowerSeries(f::Function, r::AbstractVector, Nθ::Integer,
                                 Nϕ::Integer, N::Integer, M::Integer)

#
    θ = (0:Nθ-1) .* (2π/Nθ)
    ϕ = (0:Nϕ-1) .* (2π/Nϕ)
    a = [f(ri, θi, ϕi) for ri in r, θi in θ, ϕi in ϕ]
    fit_SpectralPowerSeries(a, r, N, M)
end


function surface_integrate!(B::FluxPowerSeries{T}, A::SpectralPowerSeries{T}, ii::Integer) where{T}
    @assert mod(ii, 2) == 1
    B.a[ii] = (2π)^2 * A[ii].a[1,1];
end
