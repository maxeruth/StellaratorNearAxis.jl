"""
    Coil

A coil consists of a filament and a current. It can be constructed by `Coil(r0, J, Mcoil)`, where 
`r0` is a 3-vector of SpectralPowerSeries, `J` is a number expressing the current per unit length, 
and `Mcoil` is the number of quadrature nodes for Biot-Savart.
"""
struct Coil
    r0_s::Vector        # The coil location (SpectralPowerSeries)
    X_c::AbstractArray  # The coil location (SpatialPowerSeries, for `evaluate`)
    dX_c::AbstractArray # The coil location derivative (SpatialPowerSeries, for `evaluate`)
    J::Number  # The current in the coil

    function Coil(r0_s, J::Number, Mcoil::Integer)
        to_Spatial = (f) -> SpatialPowerSeries(f; M=Mcoil);
        r0_c = to_Spatial.(r0_s)
        dr0_c = to_Spatial.(s_deriv.(r0_s))
        
        X = Vector{SVector{3,Float64}}(undef, Mcoil);
        dX = Vector{SVector{3,Float64}}(undef, Mcoil);

        for jj = 1:Mcoil
            X[jj]  = @SVector [ r0_c[1][1].a[jj], r0_c[2][1].a[jj], r0_c[3][1].a[jj]]
            dX[jj] = @SVector [dr0_c[1][1].a[jj],dr0_c[2][1].a[jj],dr0_c[3][1].a[jj]]
        end

        new(r0_s,X,dX,J)
    end
end

function get_current(current, scale, objs)
    cls = current["@class"]
    if cls == "ScaledCurrent"
        scale = scale * current["scale"]
        current = objs[current["current_to_scale"]["value"]]
        return get_current(current, scale, objs)
    elseif cls == "CurrentSum"
        current_a = objs[current["current_a"]["value"]]
        current_b = objs[current["current_b"]["value"]]
        return get_current(current_a, scale, objs) + get_current(current_b, scale, objs)
    end

    return scale * current["current"]
end

"""

Input:
- `coil_file`: A JSON file holding the output of a Simsopt coil optimization
- `Mcoil`: The number of coil quadrature nodes for Biot-Savart

Output:
- The coil set contained in `coil_file`, which can be used to initialize a near axis expansion.
"""
function load_coils(coil_file::AbstractString, Mcoil::Number)
    coils_py = JSON.parsefile(coil_file);
    coils = Coil[]
    objs = coils_py["simsopt_objs"]

    for (key, value) in objs
        if key[1:4] == "Coil"
            current = objs[value["current"]["value"]]
            J = get_current(current, 1.0, objs)
            
            val = value["curve"]["value"]
            curve = objs[val]
            flip = false;
            phi = 0.0;
            
            if val[1:12] == "RotatedCurve"
                phi = curve["phi"]
                flip = curve["flip"]
                val = curve["curve"]["value"]
            end
            
            curve = objs[val]
            dofs = objs[curve["dofs"]["value"]]
            coeffs = dofs["x"]["data"]
            Ms_coil = length(coeffs)÷3;
            coeffs = reshape(coeffs, Ms_coil, 3)
            
            rotmat = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

            if flip
                rotmat = rotmat * Diagonal([1, -1, -1])
            end
            coeffs = coeffs * rotmat
            
            r0_coil = [zero_SpectralPowerSeries(Ms_coil, 1) for ii = 1:3]
            for ii = 1:3
                r0_coil[ii][1].a[:] = coeffs[:,ii]
            end

            push!(coils, Coil(r0_coil, J, Mcoil))
        end
    end

    coils
end




"""
    evaluate(x::AbstractArray, c::Coil)

Get the magnetic field from the coil `c` at points `x` in a 3×N array. Uses Biot-Savart discretized 
via the trapezoidal rule.
"""
function evaluate(x::AbstractArray{T}, c::Coil) where {T}
    K,N = size(x);
    @assert K == 3

    X  = c.X_c;
    dX = c.dX_c
    M = length(X)

    Xn = zeros(T, 3, M);
    B = zeros(T, 3, N)
    Bn = zeros(T, 3, M)
    fac = c.J * 10^-7;

    for (n, xn) in enumerate(eachcol(x))
        Xn = @SVector [xn[1],xn[2],xn[3]]
        Bn = @SVector [0., 0., 0.]
        for jj = 1:M
            v = X[jj] - Xn
            norm_jj = norm(v)
            Bn = Bn + StaticArrays.cross(v, dX[jj])/norm_jj^3
        end
        B[:,n] = Bn .* (fac*(2π/M))
    end
    
    B # .* (fac .* (2π/M))
end

function evaluate(x::AbstractVector, c::Coil)
    vec(evaluate(reshape(x,3,1), c))
end

function evaluate(x::AbstractArray, cs::Vector{Coil})
    B = evaluate(x,cs[1]);
    for c in cs[2:end]
        B = B + evaluate(x,c);
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
    magnetic_trajectory(cs::Vector{Coil}, x0::AbstractVector, S::Number; tol::Number=1e-8)

Evolve the magnetic field generated by the coil set `cs`
starting at the point `x0` a distance `S` to tolerance `tol`. Outputs a OrdinaryDiffEq
solution object. It also returns the times at which the trajectory passes through the 
Poincare surface intersects with the initial toroidal coordinate.
"""
function magnetic_trajectory(cs::Vector{Coil}, x0::AbstractVector, S::Number; tol::Number=1e-8)
    function f!(dx::AbstractVector, x::AbstractVector, p, s)
        B = evaluate(x, cs);
        dx[:] = B ./ sqrt(B'*B)
    end

    # int is the integrator
    function condition(x, t, int)
        angle((x[1] + im*x[2])*exp(-im*int.p[1]))
    end

    function affect!(int)
        x = int.u
        phi = angle(x[1] + im*x[2])
        if (phi-int.p[1])^2 < 0.01
            push!(int.p[2], [sqrt(x[1]^2 + x[2]^2), x[3], angle(x[1] + im*x[2])])
        end
    end

    phi0 = angle(x0[1] + im * x0[2])
    params = (phi0, [[sqrt(x0[1]^2 + x0[2]^2), x0[3], phi0]])
    sspan = (0, S);

    cb = ContinuousCallback(condition, affect!)

    prob = ODEProblem(f!, x0, sspan, params);
    
    sol = solve(prob, Tsit5(), reltol = tol, abstol = tol, callback=cb)
    sol , hcat(sol.prob.p[2]...)
end

function magnetic_map_trajectory(cs::Vector{Coil}, x0::AbstractVector; tol::Number=1e-8, Smax::Number=100.)
    function f!(dx::AbstractVector, x::AbstractVector, p, s)
        B = evaluate(x, cs);
        dx[:] = B ./ sqrt(B'*B)
    end

    # int is the integrator
    function condition(x, t, int)
        angle((x[1] + im*x[2])*exp(-im*int.p[1]))
    end

    function affect!(int)
        x = int.u
        phi = angle(x[1] + im*x[2])
        if (phi-int.p[1])^2 < 0.01
            push!(int.p[2], [sqrt(x[1]^2 + x[2]^2), x[3]])
            terminate!(int)
        end
    end

    phi0 = angle(x0[1] + im * x0[2])
    params = (phi0, [[sqrt(x0[1]^2 + x0[2]^2), x0[3]]])
    sspan = (0, Smax);

    cb = ContinuousCallback(condition, affect!)

    prob = ODEProblem(f!, x0, sspan, params);
    
    solve(prob, Tsit5(), reltol = tol, abstol = tol, callback=cb)
end

function magnetic_map(cs::Vector{Coil}, x0::AbstractVector; tol::Number=1e-8, Smax::Number=100.)
    sol = magnetic_map_trajectory(cs, x0; tol, Smax)
    sol.u[end]
end

"""
    find_magnetic_axis(x0::AbstractVector, cs::Vector{Coil})

Find a point on the magnetic axis of a coil set `cs` with `Mcoil` quadrature points. Requires an 
initial guess of a point on the axis `x0`. 
"""
function find_magnetic_axis(x0::AbstractVector, cs::Vector{Coil}; tol::Number=1e-12)
    phi0 = angle( x0[1] + im * x0[2] );
    x0 = [sqrt(x0[1]^2 + x0[2]^2), x0[3]]
    function f!(F, x)
        xF =  magnetic_map(cs, [x[1]*cos(phi0), x[1]*sin(phi0), x[2]]; tol)
        F[:] = [sqrt(xF[1]^2 + xF[2]^2), xF[3]] - x
    end
    
    sol = nlsolve(f!, x0, autodiff = :forward, ftol=tol)

    xFf = sol.zero
    x0 = [xFf[1]*cos(phi0), xFf[1]*sin(phi0), xFf[2]]

    return x0
end

function axis_from_point(x0::AbstractVector, cs::Vector{Coil}, Ms::Integer, Mc::Integer; tol=1e-8, Smax::Number=100.)
    sol = magnetic_map_trajectory(cs, x0; tol, Smax)

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

    r0
end


"""
    get_field_on_axis(r0::AbstractVector, cs::Vector{Coil}, Mc::Number, N::Number, Mcoil::Number)

Expand the magnetic field from the coil set `cs` about the axis `r0` (found, e.g., via 
[`find_magnetic_axis`](@ref)). The field is computed on `Mc` collocation points and converted to
a SpectralPowerSeries of order `Ms` dictated by the value of `Ms` used to define `r0`. The field
is computed to order `N` and with `Mcoil` quadrature points on the coils.

Output:
- `B`: The magnetic field in Frenet-Serret coordinates
- `r`: The Frenet-Serret coordinates
"""
function get_field_on_axis(r0::AbstractVector, cs::Vector{Coil}, Mc::Number, N::Number, Mcoil::Number)
    # Get the coordinate system
    Ms = get_M(r0[1])

    to_Spectral = (x) -> SpectralPowerSeries(x, M=Ms);
    to_Spatial = (x) -> SpatialPowerSeries(x, M=Mc);

    r0_s = r0;
    r0_c, ellp_s, ellp_c, ellp_ave, kappa_s, kappa_c, tau_s, tau_c, Q_s, Q_c = get_kappatau(r0_s, Mc)
    

    x_s = zero_SpectralPowerSeries(Ms, 2);
    x_s[2].a[1,1] = 1.;
    x_c = to_Spatial(x_s);
    y_s = similar(x_s);
    y_s[2].a[1,2] = 1.
    y_c = to_Spatial(y_s);
    ρ = PowerSeriesRho()

    hinv = inv(*(ellp_c, 1. - *(kappa_c, x_c; N=2); N=N))

    # ρ_tst = 0.1; θ_tst = 0.1; s_tst = 0.1;
    # println("hinv = $(evaluate(to_Spectral(hinv), ρ_tst, θ_tst, s_tst)[1])")
    # kappa_tst = evaluate(kappa_s, ρ_tst, θ_tst, s_tst)[1]
    # ellp_tst = evaluate(ellp_s, ρ_tst, θ_tst, s_tst)[1]
    # x_tst = evaluate(x_s, ρ_tst, θ_tst, s_tst)[1]
    # println("kappa = $kappa_tst, ellp = $ellp_tst")
    # println("hinv_tst = $(1/(ellp_tst*(1-x_tst*kappa_tst)))")

    r = [+(r0_c[ii], *(Q_c[ii, 1],x_c;N) + *(Q_c[ii,2],y_c;N) ; N) for ii = 1:3]
    


    # Initialize the field
    B = [zero_SpatialPowerSeries(Mc, N) for ii = 1:3]

    for c in cs
        # Get the necessary coil information
        c0_s = c.r0_s;
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
            den = dot(xi,xi)^(3/2)

            B = B + (fac * (2π/Mcoil)) * num ./ den
        end
    end

    QB = [zero_SpatialPowerSeries(Mc, N) for ii = 1:3]
    for ii=1:3, jj=1:3
        QB[ii] = QB[ii] + *(Q_c[jj,ii],B[jj];N)
    end

    zps = ZeroPowerSeries() # zero_SpatialPowerSeries(Mc, N)
    Q_polar = [x_c/ρ       y_c/ρ     zps;
               -y_c/(ρ^2)  x_c/(ρ^2) (*(-ellp_c*tau_c,hinv;N));
               zps         zps        hinv]
    QB_polar = [zero_SpatialPowerSeries(Mc, N) for ii = 1:3]

    # Test stuff
    Qinv_polar = [x_c/ρ,  (-1. * y_c),  *(ellp_c * tau_c,-1. *  y_c; N=2),
                  y_c/ρ,  x_c,          *(ellp_c*tau_c, x_c; N=2),
                  zps,    zps,          *(ellp_c, 1. - *(kappa_c, x_c; N=2); N=2)];
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

    to_Spectral.(QB_polar), to_Spectral.(r) #, to_Spectral.(Q_c), to_Spectral.(Q_polar) , to_Spectral.(Qinv_polar)
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

    phi_s = deepcopy(B_s[1]) 
    set_p0!(phi_s, 0)
    phi_s[1].a[:] .= 0.;
    phi_s[2].a[:] .= 0.;
    for n = 3:Nρ
        phi_s[n] = phi_s[n].a .* (1/(n-1))
    end

    phi_c = SpatialPowerSeries(phi_s, M=Mc)

    B_c, B_s, phi_c, phi_s, B0_s
end

"""
    field_to_nae(r0_s::AbstractVector{SpectralPowerSeries{T}}, 
        Bsup_s::AbstractVector{SpectralPowerSeries{T}}, Mc::Integer, K_reg::Integer, 
        N_reg::Integer) where {T}

Create a DirectNearAxisEquilibrium object with the magnetic field from a coil expansion. 
"""
function field_to_nae(r0_s::AbstractVector{SpectralPowerSeries{T}}, 
                      Bsup_s::AbstractVector{SpectralPowerSeries{T}}, Mc::Integer, 
                      K_reg::Number, N_reg::Integer) where {T}
    r0_s = to_arclength(r0_s, Mc); # Set 

    Ms = get_M(r0_s[1]);
    Nρ = get_N(Bsup_s[1])
    r0_c, ellp_s, ellp_c, ellp_ave, kappa_s, kappa_c, tau_s, tau_c, Q_s, Q_c = get_kappatau(r0_s, Mc)
    g_s, g_c, ginv_s, ginv_c, rootg_s, rootg_c, rootginv_s, rootginv_c = (
       get_FS_metric(Nρ,ellp_c,kappa_c,tau_c,Ms))

    B_c, B_s, phi_c, phi_s, B0_s = potential_from_field(Bsup_s, g_c, Mc)

    dphi_s = [zero_SpectralPowerSeries(T, Ms, Nρ) for ii = 1:3]
    BK_s = [zero_SpectralPowerSeries(T, Ms, Nρ) for ii = 1:3];
    BK_c = [SpatialPowerSeries(B_si; M = Mc) for B_si in B_s]

    psi_s = zero_SpectralPowerSeries(T, Ms, Nρ)
    psi_c = zero_SpatialPowerSeries(T, Mc, Nρ)

    xi_s = [zero_SpectralPowerSeries(T, Ms, Nρ) for ii = 1:2]
    xi_c = [zero_SpatialPowerSeries(T, Ms, Nρ)  for ii = 1:2]

    iota = zero_FluxPowerSeries(T, Nρ)

    # display(typeof.([]))
    nae = DirectNearAxisEquilibrium(Ms, Mc, Nρ, r0_s, r0_c, ellp_s, ellp_c, ellp_ave, kappa_s, kappa_c, tau_s,
             tau_c, Q_s, Q_c,
             g_s, g_c, ginv_s, ginv_c, rootg_s, rootg_c, rootginv_s, rootginv_c,
             B0_s, B_s, B_c, phi_s, phi_c, dphi_s, BK_s, BK_c, psi_s, psi_c, xi_s, xi_c, iota,
             K_reg, N_reg)

    
    # update_BK!(nae)
    nae.BK_c[:] = nae.ginv_c*nae.B_c
    nae.BK_s[:] = [SpectralPowerSeries(BK_i, M=Ms) for BK_i in nae.BK_c]
    # nae.BK_s[:] .= nae.B_s
    # nae.BK_c[:] .= nae.B_c
    
    nae
end

