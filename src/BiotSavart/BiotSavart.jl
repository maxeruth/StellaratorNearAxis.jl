"""
    Coil

A coil consists of a filament and a current. It can be constructed by
`Coil(r0, J)`, where `r0` is a 3-vector of SpectralPowerSeries and 
`J` is a number expressing the current per unit length.
"""
struct Coil
    r0::Vector # The coil location
    J::Number  # The current in the coil
end

"""
    evaluate(x::AbstractArray, c::Coil, M::Number)

Get the magnetic field from the coil `c` at points `x` in a 3×N array.
Uses Biot-Savart discretized via the trapezoidal rule with `M` quadrature
points.
"""
function evaluate(x::AbstractArray{T}, c::Coil, M::Integer) where {T}
    to_Spatial = (f) -> SpatialPowerSeries(f; M=M);

    K,N = size(x);
    @assert K == 3

    r0_s = c.r0;
    dr0_s = s_deriv.(r0_s);

    r0_c = to_Spatial.(r0_s);
    dr0_c = to_Spatial.(dr0_s);

    X = zeros(3, M);
    dX = zeros(3, M);

    for ii = 1:3
        X[ii,:] = r0_c[ii][1].a
        dX[ii, :] = dr0_c[ii][1].a;
    end

    Xn = zeros(T, 3, M);
    B = zeros(T, 3, N)
    Bn = zeros(T, 3, M)
    fac = c.J * 10^-7;

    for (n, xn) in enumerate(eachcol(x))
        for ii = 1:3
            Xn[ii,:] = X[ii,:] .- xn[ii]
        end

        for jj = 1:M
            norm_jj = sqrt(Xn[:,jj]'*Xn[:,jj])
            Bn[:,jj] = LinearAlgebra.cross(Xn[:,jj], dX[:,jj])/norm_jj^3
        end
        B[:,n] = sum(Bn, dims=2) .* (fac * (2π/M));
    end
    
    B
end

function evaluate(x::AbstractVector, c::Coil, M::Integer)
    vec(evaluate(reshape(x,3,1), c, M))
end

function evaluate(x::AbstractArray, cs::Vector{Coil}, M::Integer)
    B = evaluate(x,cs[1],M);
    for c in cs[2:end]
        B = B + evaluate(x,c,M);
    end

    B
end

# function axis_residual(r0::AbstractVector, cs::Vector{Coil}, Mc::Integer, Mcoil::Integer, x0::Number)
#     r0_c = [SpatialPowerSeries(ri, M=Mc) for ri in r0]
#     r0p_c = [SpatialPowerSeries(s_deriv(ri), M=Mc) for ri in r0]

#     X = vcat([ri[1].a[:]' for ri in r0_c]...)
#     dX = vcat([ri[1].a[:]' for ri in r0p_c]...)
#     B = evaluate(X, cs, Mcoil)
    
#     # dX_dot_B = sum(dX .* B; dims=1)
#     # dX_norm2 = sum(dX .* dX, dims=1)
#     B_norm  = sqrt.(sum(B .* B, dims=1))
#     B_over_Bnorm = B * Diagonal(1 ./ B_norm[:])

#     # angle = acos.(dX_dot_B ./ sqrt.(dX_norm2 .* B_norm2))
#     # diff_dX_norm2 = dX_norm2[:] - vcat(dX_norm2[2:end], dX_norm2[1])

#     res = sum((dX - B_over_Bnorm).^2) + (X[1,1] - x0)^2
#     println("res = $res")#, angle_res = $(sum(angle.^2))")
#     res
# end

# function wrapped_axis_residual(cs::Vector{Coil}, Ms::Integer, Mc::Integer, Mcoil::Integer, x0::Number)
#     function f(x::AbstractVector{T}) where {T}
#         r0 = [zero_SpectralPowerSeries(T, Ms, 1) for ii = 1:3]
#         y = reshape(x, Ms, 3);
#         for ii = 1:3
#             r0[ii][1].a[:] = y[:,ii]
#         end
#         return axis_residual(r0, cs, Mc, Mcoil, x0)
#     end

#     return f
# end

# """
#     find_magnetic_axis(r0::Vector{SpectralPowerSeries}, cs::Vector{Coil}, Mc::Integer, Mcoil::Integer)

# Find the magnetic axis of a coil set `cs`. Requires an initial guess `r0`. `Mc` gives the number of
# toroidal points ti 
# """
# function find_magnetic_axis(r0::AbstractVector, cs::Vector{Coil}, Mc::Integer, Mcoil::Integer)
#     Ms = get_M(r0[1]);
#     x0 = evaluate(r0[1], 0., 0., 0.)[1]
#     println("$(axis_residual(r0, cs, Mc, Mcoil,x0))")
#     f = wrapped_axis_residual(cs, Ms, Mc, Mcoil, x0);
#     x0 = vcat([ri[1].a[:] for ri in r0]...);

#     return optimize(f, x0, LBFGS(), Optim.Options(iterations = 10); autodiff = :forward)
# end

"""
    magnetic_trajectory(c::Vector{Coil}, x0, T; tol=1e-8)

Evolve the magnetic field generated by the coil set `cs` with `M` quadrature points
starting at the point `x0` a distance `S` to tolerance `tol`. Outputs a OrdinaryDiffEq
solution object. It also returns the times at which the trajectory passes through the 
Poincare surface intersects with the initial toroidal coordinate.
"""
function magnetic_trajectory(cs::Vector{Coil}, x0::AbstractVector, S::Number, M::Number; tol::Number=1e-8)
    function f!(dx::AbstractVector, x::AbstractVector, p, s)
        M = p[1]
        B = evaluate(x, cs, M);
        dx[:] = B ./ sqrt(B'*B)
    end

    # int is the integrator
    function condition(x, t, int)
        angle((x[1] + im*x[2])*exp(-im*int.p[2]))
    end

    function affect!(int)
        x = int.u
        ϕ = angle(x[1] + im*x[2])
        if (ϕ-int.p[2])^2 < 0.01
            push!(int.p[3], [sqrt(x[1]^2 + x[2]^2), x[3], angle(x[1] + im*x[2])])
        end
    end

    ϕ0 = angle(x0[1] + im * x0[2])
    params = [M, ϕ0, [[sqrt(x0[1]^2 + x0[2]^2), x0[3], ϕ0]]]
    sspan = (0, S);

    cb = ContinuousCallback(condition, affect!)

    prob = ODEProblem(f!, x0, sspan, params);
    
    sol = solve(prob, Tsit5(), reltol = tol, abstol = tol, callback=cb)
    sol , hcat(sol.prob.p[3]...)
end

function magnetic_map_trajectory(cs::Vector{Coil}, x0::AbstractVector, M::Number; tol::Number=1e-8, Smax::Number=100.)
    function f!(dx::AbstractVector, x::AbstractVector, p, s)
        M = p[1]
        B = evaluate(x, cs, M);
        dx[:] = B ./ sqrt(B'*B)
    end

    # int is the integrator
    function condition(x, t, int)
        angle((x[1] + im*x[2])*exp(-im*int.p[2]))
    end

    function affect!(int)
        x = int.u
        ϕ = angle(x[1] + im*x[2])
        if (ϕ-int.p[2])^2 < 0.01
            push!(int.p[3], [sqrt(x[1]^2 + x[2]^2), x[3]])
            terminate!(int)
        end
    end

    ϕ0 = angle(x0[1] + im * x0[2])
    params = [M, ϕ0, [[sqrt(x0[1]^2 + x0[2]^2), x0[3]]]]
    sspan = (0, Smax);

    cb = ContinuousCallback(condition, affect!)

    prob = ODEProblem(f!, x0, sspan, params);
    
    solve(prob, Tsit5(), reltol = tol, abstol = tol, callback=cb)
end

function magnetic_map(cs::Vector{Coil}, x0::AbstractVector, M::Number; tol::Number=1e-8, Smax::Number=100.)
    sol = magnetic_map_trajectory(cs, x0, M; tol, Smax)
    sol.u[end]
end

"""
    find_magnetic_axis(x0::AbstractVector, cs::Vector{Coil}, Mcoil::Integer)

Find a point on the magnetic axis of a coil set `cs` with `Mcoil` quadrature points. 
Requires an initial guess of a point on the axis `x0`. Note: This is 
"""
function find_magnetic_axis(x0::AbstractVector, cs::Vector{Coil}, Mcoil::Integer; tol::Number=1e-8)
    ϕ0 = angle( x0[1] + im * x0[2] );
    x0 = [sqrt(x0[1]^2 + x0[2]^2), x0[3]]
    function f!(F, x)
        xF =  magnetic_map(cs, [x[1]*cos(ϕ0), x[1]*sin(ϕ0), x[2]], Mcoil; tol)
        F[:] = [sqrt(xF[1]^2 + xF[2]^2), xF[3]] - x
    end
    
    return nlsolve(f!, x0, autodiff = :forward, ftol=tol)
end

function axis_from_point(x0::AbstractVector, cs::Vector{Coil}, Ms::Integer, Mc::Integer, 
                         Mcoil::Integer; tol=1e-8, Smax::Number=100.)
    sol = magnetic_map_trajectory(cs, x0, Mcoil; tol, Smax)

    S = sol.t[end]
    ss = (0:Mc-1) .* (S/Mc)
    r0_c = [zero_SpatialPowerSeries(Mc, 1) for ii = 1:3]
    for ii = 1:Mc
        s = ss[ii];
        ri = sol(s)
        for jj = 1:3
            r0_c[jj][1].a[ii] = ri[jj]
        end
    end

    r0 = [SpectralPowerSeries(ri, M=Ms) for ri in r0_c]

    sol, r0
end

function get_field_on_axis(r0::AbstractVector, cs::Vector{Coil}, Mc::Number, N::Number, Mcoil::Number)
    # Get the coordinate system
    Ms = get_M(r0[1])

    to_Spectral = (x) -> SpectralPowerSeries(x, M=Ms);
    to_Spatial = (x) -> SpatialPowerSeries(x, M=Mc);

    r0_s = r0;
    r0_c, ℓp_s, ℓp_c, ℓp_ave, κ_s, κ_c, τ_s, τ_c, Q_s, Q_c = get_κτ(r0_s, Mc)
    

    x_s = zero_SpectralPowerSeries(Ms, 2);
    x_s[2].a[1,1] = 1.;
    x_c = to_Spatial(x_s);
    y_s = similar(x_s);
    y_s[2].a[1,2] = 1.
    y_c = to_Spatial(y_s);
    ρ = PowerSeriesRho()

    hinv = inv(*(ℓp_c, 1. - *(κ_c, x_c; N=2); N=N))
    display(hinv)

    ρ_tst = 0.1; θ_tst = 0.1; s_tst = 0.1;
    println("hinv = $(evaluate(to_Spectral(hinv), ρ_tst, θ_tst, s_tst)[1])")
    κ_tst = evaluate(κ_s, ρ_tst, θ_tst, s_tst)[1]
    ℓp_tst = evaluate(ℓp_s, ρ_tst, θ_tst, s_tst)[1]
    x_tst = evaluate(x_s, ρ_tst, θ_tst, s_tst)[1]
    println("κ = $κ_tst, ℓp = $ℓp_tst")
    println("hinv_tst = $(1/(ℓp_tst*(1-x_tst*κ_tst)))")

    r = [+(r0_c[ii], *(Q_c[ii, 1],x_c;N) + *(Q_c[ii,2],y_c;N) ; N) for ii = 1:3]
    


    # Initialize the field
    B = [zero_SpatialPowerSeries(Mc, N) for ii = 1:3]

    for c in cs
        # Get the necessary coil information
        c0_s = c.r0;
        dc0_s = s_deriv.(c0_s);

        c0_c = [SpatialPowerSeries(ci; M=Mcoil) for ci in c0_s];
        dc0_c = [SpatialPowerSeries(ci; M=Mcoil) for ci in dc0_s];
        
        fac = c.J * 10^-7;
    
        # Evaluate the field
        for ii = 1:Mcoil
            yi = [ci[1].a[ii] for ci in c0_c]
            dyi = [dci[1].a[ii] for dci in dc0_c]

            xi = yi - r

            num = cross(xi, dyi)
            den = (xi'*xi)^(3/2)

            B = B + (fac * (2π/Mcoil)) * num ./ den
        end
    end

    QB = [zero_SpatialPowerSeries(Mc, N) for ii = 1:3]
    for ii=1:3, jj=1:3
        QB[ii] = QB[ii] + *(Q_c[jj,ii],B[jj];N)
    end

    zps = ZeroPowerSeries() # zero_SpatialPowerSeries(Mc, N)
    Q_polar = [x_c/ρ       y_c/ρ     zps;
               -y_c/(ρ^2)  x_c/(ρ^2) (*(-ℓp_c*τ_c,hinv;N));
               zps         zps        hinv]
    QB_polar = [zero_SpatialPowerSeries(Mc, N) for ii = 1:3]

    # Test stuff
    Qinv_polar = [x_c/ρ,  (-1. * y_c),  *(ℓp_c * τ_c,-1. *  y_c; N=2),
                  y_c/ρ,  x_c,          *(ℓp_c*τ_c, x_c; N=2),
                  zps,    zps,          *(ℓp_c, 1. - *(κ_c, x_c; N=2); N=2)];
    Qinv_polar = reshape(Qinv_polar, 3, 3)'

    set_p0!(QB_polar[1], -1)
    set_p0!(QB_polar[2], -2)
    # QB_polar[3] = QB[3];
    
    for ii = 1:3, jj = 1:3
        # println("ii=$ii, jj=$jj")
        # display(QB_polar[ii].p0)
        # display(Q_polar[ii,jj].p0)
        # display(QB[jj].p0)
        QB_polar[ii] = QB_polar[ii] + *(Q_polar[ii,jj], QB[jj]; N=N)
    end

    to_Spectral.(QB_polar), to_Spectral.(r), to_Spectral.(Q_c), to_Spectral.(Q_polar) , to_Spectral.(Qinv_polar)
    # B
end

function potential_from_field(Bsup_s::AbstractVector, g_c::AbstractArray, Mc::Integer)
    @assert get_p0(Bsup_s[1]) == -1
    Ms = get_M(Bsup_s[1])
    Nρ = get_N(Bsup_s[1])

    Bsup_c = [SpatialPowerSeries(Bi, M=Mc) for Bi in Bsup_s]

    # Get the magnetic field differential form
    B_c = g_c*Bsup_c
    B_s = [SpectralPowerSeries(B_ci; M = Ms) for B_ci in B_c]

    # Get the constant part
    B0_s = zero_SpectralPowerSeries(Ms, 1);
    B0_s[1].a[:] = B_s[3][1].a

    ϕ_s = deepcopy(B_s[1]) 
    set_p0!(ϕ_s, 0)
    ϕ_s[1].a[:] .= 0.;
    ϕ_s[2].a[:] .= 0.;
    for n = 3:Nρ
        ϕ_s[n] = ϕ_s[n].a .* (1/(n-1))
    end

    ϕ_c = SpatialPowerSeries(ϕ_s, M=Mc)

    B_c, B_s, ϕ_c, ϕ_s, B0_s
end


function field_to_nae(r0_s::AbstractVector{SpectralPowerSeries{T}}, 
                      Bsup_s::AbstractVector{SpectralPowerSeries{T}}, Mc::Integer, 
                      K_reg::Integer, N_reg::Integer) where {T}
    r0_s = to_arclength(r0_s, Mc); # Set 

    Ms = get_M(r0_s[1]);
    Nρ = get_N(Bsup_s[1])
    r0_c, ℓp_s, ℓp_c, ℓp_ave, κ_s, κ_c, τ_s, τ_c, Q_s, Q_c = get_κτ(r0_s, Mc)
    g_s, g_c, ginv_s, ginv_c, rootg_s, rootg_c, rootginv_s, rootginv_c = (
       get_FS_metric(Nρ,ℓp_c,κ_c,τ_c,Ms))

    B_c, B_s, ϕ_c, ϕ_s, B0_s = potential_from_field(Bsup_s, g_c, Mc)

    dϕ_s = [zero_SpectralPowerSeries(T, Ms, Nρ) for ii = 1:3]
    BK_s = [zero_SpectralPowerSeries(T, Ms, Nρ) for ii = 1:3];
    BK_c = [SpatialPowerSeries(B_si; M = Mc) for B_si in B_s]

    ψ_s = zero_SpectralPowerSeries(T, Ms, Nρ)
    ψ_c = zero_SpatialPowerSeries(T, Mc, Nρ)

    ξ_s = [zero_SpectralPowerSeries(T, Ms, Nρ) for ii = 1:2]
    ξ_c = [zero_SpatialPowerSeries(T, Ms, Nρ)  for ii = 1:2]

    ι = zero_FluxPowerSeries(T, Nρ)

    # display(typeof.([]))
    nae = DirectNearAxisEquilibrium(Ms, Mc, Nρ, r0_s, r0_c, ℓp_s, ℓp_c, ℓp_ave, κ_s, κ_c, τ_s,
             τ_c, Q_s, Q_c,
             g_s, g_c, ginv_s, ginv_c, rootg_s, rootg_c, rootginv_s, rootginv_c,
             B0_s, B_s, B_c, ϕ_s, ϕ_c, dϕ_s, BK_s, BK_c, ψ_s, ψ_c, ξ_s, ξ_c, ι,
             K_reg, N_reg)

    
    # update_BK!(nae)
    nae.BK_c[:] = nae.ginv_c*nae.B_c
    nae.BK_s[:] = [SpectralPowerSeries(BK_i, M=Ms) for BK_i in nae.BK_c]
    # nae.BK_s[:] .= nae.B_s
    # nae.BK_c[:] .= nae.B_c
    
    nae
end
