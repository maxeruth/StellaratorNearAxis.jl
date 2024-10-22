## Code to perform the direct Frenet-Serret Expansion

## The solution struct
# _s : Spectral quantity
# _c : On collocation nodes
struct DirectNearAxisEquilibrium
    # Orders
    Ms::Integer;
    Mc::Integer;
    Nρ::Integer;

    # Geometric Quantities
    r0_s::Vector   # Axis
    r0_c::Vector
    ℓp_s::SpectralPowerSeries           # Rate of travel along axis
    ℓp_c::SpatialPowerSeries
    ℓp_ave::Number                      # Working in arclength coordinates, constant rate of travel
    κ_s::SpectralPowerSeries            # Axis curvature
    κ_c::SpatialPowerSeries
    τ_s::SpectralPowerSeries            # Axis torsion
    τ_c::SpatialPowerSeries
    Q_s::Matrix       # Frenet-Serret Unit Vectors
    Q_c::Matrix
    g_s::Matrix   # Metric
    g_c::Matrix
    ginv_s::Matrix # Inverse metric
    ginv_c::Matrix
    rootg_s::SpectralPowerSeries    # Metric determinant
    rootg_c::SpatialPowerSeries
    rootginv_s::SpectralPowerSeries # Inverse metric determinant
    rootginv_c::SpatialPowerSeries

    # # Input Quantities
    B0_s::SpectralPowerSeries # Magnetic field strength on axis
    # B0_c::SpatialPowerSeries
    # # p::FluxPowerSeries      # Pressure
    #
    # # Equilibrium quantities
    B_s::Vector # Magnetic field
    B_c::Vector
    ϕ_s::SpectralPowerSeries         # Magnetic potential (for divergence free part)
    ϕ_c::SpatialPowerSeries

    dϕ_s::Vector
    BK_s::Vector # Regularized magnetic field
    BK_c::Vector 
    #
    # Derived quantities
    ψ_s::SpectralPowerSeries
    ψ_c::SpatialPowerSeries

    ξ_s::Vector
    ξ_c::Vector
    ι::FluxPowerSeries                # Rotation number

    # Regularization Quantities
    K_reg::Integer # Characteristic number of Fourier modes to resolve
    N_reg::Integer # Power law of the regularization
end

function to_arclength(r0_s::Vector{SpectralPowerSeries{T}}, Mc::Integer) where {T}
    # if Ms=-1
        
    # end
    Ms = get_M(r0_s[1])

    r0p_s = ϕ_deriv.(r0_s)
    r0p_c = [SpatialPowerSeries(ri; M = Mc) for ri in r0p_s]
    ℓp_c = norm(r0p_c)
    ℓp_s = SpectralPowerSeries(ℓp_c; M=Ms)
    
    # Solve ∂ₜ(ℓ) = ℓ'(t) to find the length of the curve
    function f_ℓp!(dt, ℓ, p, t)
        dt[:] = evaluate(ℓp_s, 0.0, 0.0, t)
    end
    ℓ0 = [0.0];
    tspan = (0.0, 2π);
    prob = ODEProblem(f_ℓp!, ℓ0, tspan)

    sol = solve(prob, Tsit5(), reltol = 1e-12, abstol = 1e-12);
    L = sol(2π);
    # display(sol)
    # display(L)
    # display(ℓp_c)


    # Solve ∂ₛt = (L/2π) * (∂ₜℓ)⁻¹
    function f_ℓp_inv!(dt, t, p, s)
        dt[:] = (L/(2π)) ./ evaluate(ℓp_s, 0.0, 0.0, t);
    end
    t0 = [0.0];
    sspan = (0.0, 2π);
    prob = ODEProblem(f_ℓp_inv!, t0, sspan);
    
    sol = solve(prob, Tsit5(), reltol = 1e-12, abstol = 1e-12);

    # Sample r0(t) at equal s-points
    # display(sol(2π).-2π)
    pts = [sol(2π*(ii/Mc))[1] * (2π/sol(2π)[1]) for ii = 0:Mc-1]
    # display(pts)
    r0_arclength_c = [zero_SpatialPowerSeries(Mc, 1) for ii=1:3]
    # display(r0_arclength_c[1])
    for ii = 1:3
        r0_arclength_c[ii][1].a[:] = evaluate(r0_s[ii], 0.0, 0.0, pts);
    end

    [SpectralPowerSeries(ri; M=Ms) for ri in r0_arclength_c]
end

function get_κτ(r0_s::Vector{SpectralPowerSeries{T}}, Mc::Integer) where {T}
    Ms = get_M(r0_s[1])
    r0_c = [SpatialPowerSeries(ri; M = Mc) for ri in r0_s]
    r0p_s = ϕ_deriv.(r0_s)
    r0p_c = [SpatialPowerSeries(ri; M = Mc) for ri in r0p_s]
    r0pp_s = ϕ_deriv.(r0p_s)
    r0pp_c = [SpatialPowerSeries(ri; M = Mc) for ri in r0pp_s]
    r0ppp_s = ϕ_deriv.(r0pp_s)
    r0ppp_c = [SpatialPowerSeries(ri; M = Mc) for ri in r0ppp_s]


    ℓp_c = norm(r0p_c);
    ℓp_s = SpectralPowerSeries(ℓp_c; M=Ms)
    ℓp_ave = sum(ℓp_c[1].a)/length(ℓp_c[1].a);

    rp_times_rpp = cross(r0p_c, r0pp_c);
    norm_rp_times_rpp = norm(rp_times_rpp);

    # κ = ||r' × r''|| / ||r'||^3
    κ_c = norm_rp_times_rpp * ℓp_c^(-3.)
    κ_s = SpectralPowerSeries(κ_c; M=Ms)

    # τ = ((r' × r'') ⋅ r''') / ||r' × r''||^2
    τ_c = dot(rp_times_rpp, r0ppp_c) * norm_rp_times_rpp^(-2.)
    τ_s = SpectralPowerSeries(τ_c; M=Ms)

    # Orthogonal Frame
    Q_c = Matrix{AbstractPowerSeries}(undef, 3, 3)
    Q_c[:, 3] = [r0p_c[ii] * inv(ℓp_c) for ii = 1:3]
    Q_c[:, 1] = cross(r0p_c, cross(r0pp_c, r0p_c))
    Q_c[:, 1] = [ni * inv(norm_rp_times_rpp) * inv(ℓp_c) for ni in Q_c[:,1]]
    Q_c[:, 2] = cross(Q_c[:,3], Q_c[:,1])

    Q_s = [SpectralPowerSeries(Q_c[ii,jj]; M=Ms) for ii = 1:3, jj = 1:3]

    r0_c, ℓp_s, ℓp_c, ℓp_ave, κ_s, κ_c, τ_s, τ_c, Q_s, Q_c
end

function initial_metric_matrix()
    g = Matrix{AbstractPowerSeries}(undef, 3, 3);
    g[1, 2] = ZeroPowerSeries();
    g[1, 3] = ZeroPowerSeries();
    g[2, 1] = ZeroPowerSeries();
    g[3, 1] = ZeroPowerSeries();
    g[1, 1] = IdentityPowerSeries();
    g
end

function get_FS_metric(Nρ::Integer, ℓp_c::SpatialPowerSeries{T},
               κ_c::SpatialPowerSeries{T}, τ_c::SpatialPowerSeries{T},
               Ms::Integer) where {T}
    ρ = PowerSeries_ρ()
    ρ2 = ρ^2;
    g_c = initial_metric_matrix()
    g_s = initial_metric_matrix()
    ginv_c = initial_metric_matrix()
    ginv_s = initial_metric_matrix()

    # g_c
    Mc = get_M(ℓp_c)
    hs = zero_SpatialPowerSeries(Mc, 2);
    hs[1] = ones(T, Mc);
    hs[2] = -κ_c[1].a * (cos.([0, 2π/3]))'

    M = get_M(κ_c);
    # g = [ρ²    , ℓ'τρ²
    #      ℓ'τρ² , (ℓ')²(hₛ² + ρ²τ²)    ]
    g_c[2,2] = ρ2*zero_SpatialPowerSeries(T, M, Nρ)
    g_c[2,2][1] = ones(T, M);
    g_c[2,3] = ρ2*(*(ℓp_c, τ_c; N=Nρ))
    g_c[3,2] = ρ2*(*(ℓp_c, τ_c; N=Nρ))
    g_c[3,3] = *(*(hs, hs; N=3) + ρ2*(*(τ_c,τ_c; N=3)), ℓp_c*ℓp_c; N=Nρ)

    # rootg
    rootg_c = ρ*((ρ^-2)*(g_c[2,2]*g_c[3,3] - g_c[2,3]*g_c[3,2]))^0.5;
    rootg_s = SpectralPowerSeries(rootg_c, M=Ms)

    # rootginv_c
    rootginv_c = inv(rootg_c)
    rootginv_s = SpectralPowerSeries(rootginv_c, M=Ms)
    rootginv_c2 = rootginv_c^2


    # ginv_c
    ginv_c[2,2] = g_c[3,3]*rootginv_c2;
    ginv_c[3,3] = g_c[2,2]*rootginv_c2;
    ginv_c[2,3] = - g_c[2,3]*rootginv_c2
    ginv_c[3,2] = deepcopy(ginv_c[2,3]);

    for ii = 2:3, jj = 2:3
        # g_s
        g_s[ii,jj] = SpectralPowerSeries(g_c[ii,jj]; M=Ms)

        # ginv_s
        ginv_s[ii,jj] = SpectralPowerSeries(ginv_c[ii,jj]; M=Ms)
    end

    g_s, g_c, ginv_s, ginv_c, rootg_s, rootg_c, rootginv_s, rootginv_c
end

function InitialVacuumNearAxisEquilibrium(r0_s::Vector{SpectralPowerSeries{T}},
               Nρ::Integer, Mc::Integer, B0_s::SpectralPowerSeries{T}, 
               K_reg::Integer, N_reg::Integer) where {T}
    # @assert mod(Nρ, 2) == 1
    r0_s = to_arclength(r0_s, Mc); # Set 

    Ms = get_M(r0_s[1]);
    r0_c, ℓp_s, ℓp_c, ℓp_ave, κ_s, κ_c, τ_s, τ_c, Q_s, Q_c = get_κτ(r0_s, Mc)
    # display(ℓp_s)
    # [display(r0_ci.a) for r0_ci in r0_c]
    g_s, g_c, ginv_s, ginv_c, rootg_s, rootg_c, rootginv_s, rootginv_c = (
       get_FS_metric(Nρ,ℓp_c,κ_c,τ_c,Ms))


    B_s = [zero_SpectralPowerSeries(T, Ms, Nρ) for ii = 1:3]
    B_s[3][1] = B0_s[1];
    B_c = [SpatialPowerSeries(B_si; M = Mc) for B_si in B_s]

    ϕ_s = zero_SpectralPowerSeries(T, Ms, Nρ)
    ϕ_c = zero_SpatialPowerSeries(T, Mc, Nρ)

    dϕ_s = [zero_SpectralPowerSeries(T, Ms, Nρ) for ii = 1:3]
    BK_s = [zero_SpectralPowerSeries(T, Ms, Nρ) for ii = 1:3]
    B_s[3][1] = B0_s[1];
    BK_c = [SpatialPowerSeries(B_si; M = Mc) for B_si in B_s]

    ψ_s = zero_SpectralPowerSeries(T, Ms, Nρ)
    ψ_c = zero_SpatialPowerSeries(T, Mc, Nρ)

    ξ_s = [zero_SpectralPowerSeries(T, Ms, Nρ) for ii = 1:2]
    ξ_c = [zero_SpatialPowerSeries(T, Ms, Nρ)  for ii = 1:2]

    ι = zero_FluxPowerSeries(T, Nρ)

    # display(typeof.([]))
    DirectNearAxisEquilibrium(Ms, Mc, Nρ, r0_s, r0_c, ℓp_s, ℓp_c, ℓp_ave, κ_s, κ_c, τ_s,
             τ_c, Q_s, Q_c,
             g_s, g_c, ginv_s, ginv_c, rootg_s, rootg_c, rootginv_s, rootginv_c,
             B0_s, B_s, B_c, ϕ_s, ϕ_c, dϕ_s, BK_s, BK_c, ψ_s, ψ_c, ξ_s, ξ_c, ι,
             K_reg, N_reg)
end

function vacuum_update_B!(nae::DirectNearAxisEquilibrium)
    # Update B0
    nae.B_s[:] = grad(nae.ϕ_s);
    nae.B_s[3][1] = nae.B_s[3][1] + nae.B0_s[1];
    nae.B_c[:] = [SpatialPowerSeries(Bi, M=nae.Mc) for Bi in nae.B_s];
end

function vacuum_residual(nae::DirectNearAxisEquilibrium)
    to_Spectral = (x) -> SpectralPowerSeries(x, M=nae.Ms);
    to_Spatial = (x) -> SpatialPowerSeries(x, M=nae.Mc);
    # co = (x)->change_order(x, nae.Nρ+1)

    B_contra_s = to_Spectral.(nae.rootg_c * (nae.ginv_c*nae.B_c));
    divB_c = to_Spatial(div(B_contra_s))
    res = PowerSeries_ρ()^(-1) * to_Spectral(*(divB_c, inv(nae.ℓp_c); N=nae.Nρ+1))
    # display(res[1].a)
    # display(res[2].a)
    # display(res[3].a)
    unsafe_distribute_p0(res)
end

function regularized_residual(nae::DirectNearAxisEquilibrium)
    println("regularized_residual: still don't trust this routine...")
    to_Spectral = (x) -> SpectralPowerSeries(x, M=nae.Ms);
    to_Spatial = (x) -> SpatialPowerSeries(x, M=nae.Mc);

    Brootg_s = to_Spectral.(nae.rootg_c * nae.BK_c);
    divB_c = to_Spatial(div(Brootg_s))
    res = PowerSeries_ρ()^(-1) * to_Spectral(*(divB_c, inv(nae.ℓp_c); N=nae.Nρ+1))

    for ii = 1:3
        println("ii=$ii")
        display(Brootg_s[ii])
    end

    # display(res[1].a)
    # display(res[2].a)
    # display(res[3].a)
    # display(res)
    unsafe_distribute_p0(res)
end

function vacuum_diagonal_operator(n::Integer)
    v = zeros(n-1);
    if mod(n, 2) == 0
        v[1] = n^2
        for jj = 2:2:n-1
            v[jj:jj+1] .= n^2 - jj^2
        end
    else
        for jj = 1:2:n-1
            v[jj:jj+1] .= n^2 - jj^2
        end
    end

    Diagonal(v)
end

# M = number of s-modes
# N = number of θ-modes
function regularization_diagonal_operator(M::Integer, N::Integer, K_reg::Integer, N_reg::Integer)
    v = ones(M,N-1); 
    
    for m = 2:2:M
        k = m÷2
        v[m:m+1,:] = v[m:m+1,:] .+ (k/K_reg)^(2*N_reg);
    end

    if mod(N, 2) == 0
        for jj = 2:2:N-1
            v[:, jj:jj+1] .= v[:, jj:jj+1] .+ (jj/K_reg)^(2*N_reg);
        end
    else
        for jj = 1:2:N-1
            v[:, jj:jj+1] .= v[:, jj:jj+1] .+ (jj/K_reg)^(2*N_reg);
        end
    end

    v
end

function update_BK!(nae)
    to_Spectral = (x) -> SpectralPowerSeries(x, M=nae.Ms);
    to_Spatial = (x) -> SpatialPowerSeries(x, M=nae.Mc);
    ρ = PowerSeries_ρ()

    dϕds_s = deepcopy(nae.ϕ_s)
    for ii = 1:nae.N_reg
        dϕds_s = -(nae.K_reg*nae.ℓp_s[1].a[1])^(-2) * ϕ_deriv(ϕ_deriv(dϕds_s))
    end

    
    ## Potential alternate method for taking the derivative
    # dϕds_s2= deepcopy(nae.ϕ_s)
    # for ii = 1:nae.Nρ
    #     dϕds_s[ii].a[1,:] .= 0.
    #     for jj = 2:2:nae.Ms
    #         dϕds_s[ii].a[jj:jj+1,:] = dϕds_s[ii].a[jj:jj+1,:] .* ((jj÷2)/nae.K_reg)^(2. * nae.N_reg)
    #     end
    # end

    # display(Fnorm(dϕds_s - dϕds_s2))

    dϕ_s = Diagonal([ρ, inv(ρ), ZeroPowerSeries()]) * grad(dϕds_s)

    nae.dϕ_s[1:2] = grad(dϕds_s)[1:2]
    dϕ_c = [to_Spatial(dϕ_s[1]), to_Spatial(dϕ_s[2]), ZeroPowerSeries()]

    # nae.BK_c[:] = dϕ_c + nae.ginv_c*nae.B_c
    nae.BK_c[:] = *(nae.ℓp_c, nae.rootginv_c, N=nae.Nρ) * dϕ_c + nae.ginv_c*nae.B_c


    # nae.BK_c[:] = nae.ginv_c*nae.B_c
    nae.BK_s[:] = to_Spectral.(nae.BK_c)
end

function get_IC_from_eq(nae::DirectNearAxisEquilibrium)
    Nρ = nae.Nρ
    Ms = nae.Ms

    ϕ_IC = zeros(Ms, 2Nρ-1)
    ϕ_s = nae.ϕ_s

    ϕ_IC[:, 1]   = ϕ_s[1].a[:];
    ϕ_IC[:, 2:3] = ϕ_s[2].a[:];
    
    for ii = 3:Nρ
        if mod(ii,2) == 1
            ϕ_IC[:, 2ii-2:2ii-1] = ϕ_s[ii].a[:, ii:-1:ii-1];
        else
            ϕ_IC[:, 2ii-2:2ii-1] = ϕ_s[ii].a[:, ii-1:ii];
        end
    end

    ϕ_IC
end

"""
    vacuum_solve(nae::DirectNearAxisEquilibrium, ϕ_IC::AbstractArray)

Solve for a vacuum field. The boundary conditions at each order are given boundary 
 - `ϕ_IC = [ϕ_0c0, ϕ_1c1, ϕ_1s1, ϕ_2c2, ϕ_2s2, ..., ϕ_NcN, ϕ_NsN]`
"""
function vacuum_solve(nae::DirectNearAxisEquilibrium, ϕ_IC::AbstractArray)
    # Ms = nae.Ms;
    # Mc = nae.Mc;
    # Nρ = nae.Nρ;

    for ii = 1:nae.Nρ
        nae.ϕ_s[ii].a[:] .= 0.;
        nae.ϕ_c[ii].a[:] .= 0.;
        for jj = 1:3
            nae.B_s[jj][ii].a[:] .= 0.;
            nae.B_c[jj][ii].a[:] .= 0.
        end
    end


    nae.ϕ_s[1].a[:] = ϕ_IC[:, 1];
    nae.ϕ_s[2].a[:] = ϕ_IC[:, 2:3];
    
    # Update B0
    vacuum_update_B!(nae)
    # update_BK!(nae)

    for ii = 3:nae.Nρ
        # Set ϕ_s to zero for the residual

        # Get residual
        res = vacuum_residual(nae)[ii-2];
        
        # Solve diagonal part of system
        D = vacuum_diagonal_operator(ii-1)
        D_reg = regularization_diagonal_operator(nae.Ms, ii-1, nae.K_reg, nae.N_reg)
        (nae.ϕ_s[ii]).a[:,1:ii-2] = -(res.a * inv(D)) ./ D_reg;

        if mod(ii,2) == 1
            nae.ϕ_s[ii].a[:, ii:-1:ii-1] = ϕ_IC[:, 2ii-2:2ii-1];
        else
            nae.ϕ_s[ii].a[:, ii-1:ii]    = ϕ_IC[:, 2ii-2:2ii-1];
        end

        # Update B0
        vacuum_update_B!(nae)
        # update_BK!(nae)
    end
    update_BK!(nae)

end



function ψ_residual(nae::DirectNearAxisEquilibrium)
    to_Spectral = (x) -> SpectralPowerSeries(x, M=nae.Ms);
    to_Spatial = (x) -> SpatialPowerSeries(x, M=nae.Mc);

    dψ_s = nabla(nae.ψ_s)
    dψ_c = to_Spatial(dψ_s)
    res_c = dot(nae.B_c, dψ_c)

    return [to_Spectral(res_ci) for res_ci in res_c]
end

function ψ_matrices(nae::DirectNearAxisEquilibrium)

    Ms = nae.Ms;
    Mc = nae.Mc;
    x = (0:Ms-1).*(2π/Mc);
    F = full_fourier_matrix(x, Ms)
    # ϕ20, ϕ21, ϕ22 and derivative
end


function magnetic_trajectory(nae::DirectNearAxisEquilibrium, ρ0, θ0, s0, sf; tol=1e-8)
    BK_s = nae.BK_s

    function f!(dx, x, p, s)
        ρ = sqrt(x[1]^2 + x[2]^2);
        θ = angle(x[1] + im*x[2]);
        B = [evaluate(BKi, ρ, θ, s)[1] for BKi in BK_s];
        ρdot = B[1]/B[3]
        θdot = B[2]/B[3]

        dx[:] = [ρdot*cos(θ)-θdot*ρ*sin(θ), ρdot*sin(θ)+θdot*ρ*cos(θ)]
        # println("s=$s, ρdot=$ρdot, θdot = $θdot") ;flush(stdout)
        # println("s=$s, x=$(x), B=$(B), dx=$(dx)");flush(stdout)
    end

    sspan = (s0, sf);
    x0 = [ρ0*cos(θ0), ρ0*sin(θ0)]
    prob = ODEProblem(f!, x0, sspan);
    
    solve(prob, Tsit5(), reltol = tol, abstol = tol)
end

function evolve_BK_map(nae::DirectNearAxisEquilibrium, s0::Number, sf::Number; tol=1e-12)
    BK = nae.BK_s;
    N = get_N(BK[1]);

    NF = ((N-2)*(N-1))÷2 - 1;   

    # Initialize position as identity
    x0 = zeros(NF, 2);
    F0_s = [zero_SpectralPowerSeries(1,2) for ii = 1:2];
    F0_s[1][2].a[1] = 1.;
    F0_s[2][2].a[2] = 1.;

    x0[1:2,1] = F0_s[1][2].a[:]; 
    x0[1:2,2] = F0_s[2][2].a[:];

    x0 = x0[:];

    
    xF = zero_SpectralPowerSeries(1,2)
    xF[2].a[1] = 1.;
    xF = SpatialPowerSeries(xF)
    y = zero_SpectralPowerSeries(1,2)
    y[2].a[2] = 1.;
    y = SpatialPowerSeries(y)
    ρ = PowerSeries_ρ();


    function f!(dx, x, p, t)
        # Cast x as a SpatialPowerSeries on a section
        xt = reshape(x, NF, 2)
        Ft = [zero_SpectralPowerSeries(1,N-2) for ii = 1:2]
        for ii = 1:2
            kk = 0;
            for jj = 2:N-2
                Ft[ii][jj].a[:] = xt[kk+1:kk+jj, ii]
                kk = kk+jj
            end
        end
        
        # Evaluate the magnetic field dynamics at the slice
        Bpolar = [SpatialPowerSeries(section(BKi, t)) for BKi in BK]
        Bρ = Bpolar[1]; Bθ = Bpolar[2]; Bs = Bpolar[3];
        Bsinv = inv(Bs);
        Bcart = [(*(Bρ/ρ, xF; N=N) - *(Bθ, y; N=N))*Bsinv, 
                 (*(Bρ/ρ, y; N=N) + *(Bθ, xF; N=N))*Bsinv]
        Bcart_s = SpectralPowerSeries.(Bcart)
        
        # Compose the magnetic with the current position to get the derivative
        Ftree = section_composition_basis(Ft)
        f_c = [section_compose(unsafe_distribute_p0(SpectralPowerSeries(Bcarti)), Ftree) for Bcarti in Bcart]

        # Fill dx with the derivative
        ft = zeros(NF,2);
        for ii = 1:2
            kk = 0;
            for jj = 2:N-2
                ft[kk+1:kk+jj, ii] = f_c[ii][jj].a[:]
                kk = kk+jj
            end
        end

        dx[:] = ft[:]
    end

    sspan = (s0, sf);
    prob = ODEProblem(f!, x0, sspan);
    
    sol = solve(prob, Tsit5(), reltol = tol, abstol = tol)

    xf = reshape(sol(sf), NF, 2);
    F = [zero_SpectralPowerSeries(1,N-2) for ii = 1:2]
    for ii = 1:2
        kk = 0;
        for jj = 2:N-2
            F[ii][jj].a[:] = xf[kk+1:kk+jj, ii]
            kk = kk+jj
        end
    end


    F
end


function α_map_residual(G::AbstractVector, ι::FluxPowerSeries, ψ::SpatialPowerSeries, Ftree::AbstractVector)
    res = [SpatialPowerSeries(section_compose(Gi, Ftree)) for Gi in G];
    N = get_N(ψ);
    ρ = PowerSeries_ρ();
    
    ιψ = zero_SpectralPowerSeries(1,N+1)
    ιψ[1] = [2π*ι[1]]
    ιψ = SpatialPowerSeries(ιψ)

    G_c = SpatialPowerSeries.(G; N)

    ψn = ψ
    for n = 3:2:N-2
        ιψ = ιψ + 2π*ι[n] * ψn
        ψn = ψn*ψ
    end
    s, c = sincos(ιψ);

    RG = [c*G_c[1] - s*G_c[2], s*G_c[1] + c*G_c[2]]
    # println("res = $(hcat([resi[4].a for resi in res]...))")
    # println("RG  = $(hcat([Gi[2].a for Gi in RG]...))")
    # println("c[1] = $(c[1].a[1]), s[1] = $(s[1].a[1])")
    res = res - RG
    # res[1] = res[1] - (c*G_c[1] - s*G_c[2]);
    # res[2] = res[2] - (s*G_c[1] + c*G_c[2]);
    
    SpectralPowerSeries.(res)
end

"""
    find_α_map(F::AbstractVector)

Find a map to coordinates where `F` is conjugate to a rotation.
"""
function find_α_map(F::AbstractVector)
    N = get_N(F[1]);
    Ftree = section_composition_basis(F)
    G = [zero_SpectralPowerSeries(1,N) for ii = 1:2]
    ι = zero_FluxPowerSeries(N-1)
    ψ = zero_SpatialPowerSeries(1,N+1);
    
    F0p = zeros(2,2);
    F0p[1,:] = Ftree[2][1][2].a;
    F0p[2,:] = Ftree[2][2][2].a;

    λ, P = eigen(F0p);
    ι0 = -angle(λ[1])/(2π)
    ι[1] = ι0
    
    QU = [1 -im
          1  im] ./sqrt(2)
    U = real.(P*QU);
    R = [cos(2π*ι0) -sin(2π*ι0)
         sin(2π*ι0)  cos(2π*ι0)];
    Rp = [-sin(2π*ι0) -cos(2π*ι0)
           cos(2π*ι0) -sin(2π*ι0)];
    Uinv = inv(U);
    # println("U*R*Uinv - F0p = $(U*R*Uinv-F0p)")
    # println("R*Uinv = $(R*Uinv)")
    G[1][2] = Uinv[1,:];
    G[2][2] = Uinv[2,:];
    
    G_c = SpatialPowerSeries.(G);
    ψ = *(G_c[1],G_c[1];N=N+1) + *(G_c[2],G_c[2];N=N+1)
    
    RpG0 = SpatialPowerSeries.([Rp[1,1]*G[1]+Rp[1,2]*G[2], Rp[2,1]*G[1] + Rp[2,2]*G[2]])
    # println("Rp  - $(Rp)")
    # println("F   = $(vcat([Fi[2].a for Fi in F]...))")
    # println("G   = $(vcat([Gi[2].a for Gi in G]...))")
    # println("G∘F = $(vcat([section_compose(Gi,Ftree)[2].a for Gi in G]...))")
    
    for n = 3:N
        # println("n=$n")
        sz = isodd(n) ? 2n : 2n+1
        A = zeros(2n,sz)
        b = zeros(2n);

        A[1:n,1:n] = hcat([Fi[n].a[:] for Fi in Ftree[n]]...)  - R[1,1]*I
        A[n+1:2n, n+1:2n] = A[1:n, 1:n]
        A[1:n, n+1:2n] = -R[1,2]*Matrix(1.0*I, n, n)
        A[n+1:2n, 1:n] = -R[2,1]*Matrix(1.0*I, n, n)

        res = α_map_residual(G, ι, ψ, Ftree)
        b[1:n] = -res[1][n].a;
        b[n+1:2n] = -res[2][n].a;
        
        if isodd(n)
            coef = A\b
            G[1][n] = coef[1:n];
            G[2][n] = coef[n+1:2n];
        else
            ψn = ψ
            for jj = 6:2:n
                # println("jj=$jj")
                ψn = ψn*ψ
            end
            RpG0ψn = [SpectralPowerSeries(Ai*ψn) for Ai in RpG0];
            AK[1:n,2n+1]    = - 2π*RpG0ψn[1][n].a
            A[n+1:2n,2n+1] = - 2π*RpG0ψn[2][n].a

            coef = A\b
            G[1][n] = coef[1:n];
            G[2][n] = coef[n+1:2n];
            ι[n-1]  = coef[2n+1];
            # display(coef)
            λ,P = eigen(A[1:2n,1:2n]);
            # display(eigen(A[1:2n,1:2n]))
            # display(inv(P)*b)
            # display(inv(P)*A[:,2n+1])
            # display((inv(P)*b) ./ (inv(P)*A[:,2n+1]))
        end

        G_c = SpatialPowerSeries.(G);
        ψ = *(G_c[1],G_c[1];N=N+1) + *(G_c[2],G_c[2];N=N+1)
    end

    G, ι, SpectralPowerSeries(ψ)
end

function flux_leading_order(hx_c, hy_c, nae)
    Ms = nae.Ms;
    Mc = nae.Mc;
    to_Spectral = (x) -> SpectralPowerSeries(x, M=Ms);
    to_Spatial = (x) -> SpatialPowerSeries(x, M=nae.Mc);

    F = full_fourier_matrix(Mc, Ms)
    W = fourier_weight_matrix(Mc, Ms)
    D = ϕ_deriv_operator(Ms)

    Fθ = half_fourier_matrix(2,2)
    hx_c_a2 = hx_c[2].a * inv(Fθ')
    hy_c_a2 = hy_c[2].a * inv(Fθ')

    A = zeros(4Ms, 4Ms);
    Hxx = W\F'*Diagonal(hx_c_a2[:,1])*F
    Hxy = W\F'*Diagonal(hx_c_a2[:,2])*F
    Hyx = W\F'*Diagonal(hy_c_a2[:,1])*F
    Hyy = W\F'*Diagonal(hy_c_a2[:,2])*F

    J  = [0. -1.; 1. 0.]

    i1 = 1:2Ms;
    i2 = 2Ms+1:4Ms;
    A[i1,i1] = kron(J', Hxx+D)
    A[i1,i2] = kron(J', Hyx)
    A[i2,i1] = kron(J', Hxy)
    A[i2,i2] = kron(J', Hyy+D)
    # A = -A

    λ, P = eigen(A);
    p = sortperm(abs.(λ))
    λ = λ[p[1:4]];
    P = P[:,p[1:4]];

    p = sortperm(real.(λ))
    λ = λ[p]
    P = P[:,p]
    v1 = P[:,1]
    v2 = P[:,2]

    G0 = reshape(real(v1[1:Ms:end]),2,2);
    if det(G0) > 0
        λ = λ[1]
        v1 = reshape(v1,Ms,2,2)
        v2 = reshape(v2,Ms,2,2)
    else
        λ = λ[3]
        v1 = reshape(P[:,3],Ms,2,2)
        v2 = reshape(P[:,4],Ms,2,2)
    end
    A = [v1[1,1,1] v2[1,1,1]; v1[1,1,2] v2[1,1,2]]
    # display(A)
    # display(norm(v1-v2))
    # display(P[:,1])
    # display(P[:,2])
    a = A\[1., 0.]
    v = real.(a[1]*v1 + a[2]*v2)

    G0_s = [zero_SpectralPowerSeries(Ms, 1) for ii = 1:2, jj = 1:2]
    for ii = 1:2, jj = 1:2
        G0_s[ii,jj][1].a[:] = v[:,ii,jj]
    end
    G0_c = to_Spatial.(G0_s);
    detG0_c = G0_c[1,1]*G0_c[2,2] - G0_c[1,2]*G0_c[2,1];
    γ = π * sum((nae.B_c[3] / detG0_c)[1].a) / Mc;
    v = v*sqrt(γ);

    # G0_s = [zero_SpectralPowerSeries(Ms, 1) for ii = 1:2, jj = 1:2]
    # for ii = 1:2, jj = 1:2
    #     G0_s[ii,jj][1].a[:] = v[:,ii,jj]
    # end
    # G0_c = to_Spatial.(G0_s);
    # detG0_c = G0_c[1,1]*G0_c[2,2] - G0_c[1,2]*G0_c[2,1];
    # display( Float64(π) * nae.B_c[3]/detG0_c)



    # TODO: what if it's not elliptic?

    nae.ξ_s[1][2] = v[:,1,:];
    nae.ξ_s[2][2] = v[:,2,:];
    nae.ξ_c[:] .= to_Spatial.(nae.ξ_s)
    ψ_c= dot(nae.ξ_c,nae.ξ_c);
    ψ_s = to_Spectral(ψ_c)
    for ii = 1:get_N(ψ_c)
        nae.ψ_c[ii] = ψ_c[ii]
        nae.ψ_s[ii] = ψ_s[ii]
    end
    nae.ι[1] = real(λ)
end

function flux_resid(hx_c, hy_c, nae)
    Mc = nae.Mc;
    Ms = nae.Ms;
    to_Spectral = (x) -> SpectralPowerSeries(x, M=Ms);
    to_Spatial = (x) -> SpatialPowerSeries(x, M=Mc);
    Nρ = nae.Nρ
    
    ξ_s = nae.ξ_s
    ξ_c = nae.ξ_c
    ψ_c = nae.ψ_c
    ι = nae.ι

    Gpolar = to_Spatial.([ρ_deriv(ξ_s[1]) θ_deriv(ξ_s[1]); 
                          ρ_deriv(ξ_s[2]) θ_deriv(ξ_s[2])])


    x_s = zero_SpectralPowerSeries(Ms, 2);
    x_s[2].a[1,1] = 1.;
    x_c = to_Spatial(x_s);
    y_s = similar(x_s);
    y_s[2].a[1,2] = 1.
    y_c = to_Spatial(y_s);
    ρ = PowerSeries_ρ()

    polartocart = [x_c/ρ     y_c/ρ; 
                  -y_c/(ρ^2) x_c/(ρ^2)];

    G = to_Spatial.(unsafe_distribute_p0.(to_Spectral.(
                    [*(Gpolar[ii,1],polartocart[1,jj],N=Nρ+1) for ii = 1:2, jj=1:2] 
                  + [*(Gpolar[ii,2],polartocart[2,jj],N=Nρ+1) for ii = 1:2, jj=1:2])));


    Gh = [*(G[ii,1],hx_c,N=Nρ) + 
          *(G[ii,2],hy_c,N=Nρ) for ii = 1:2]

    Jξ = [-ξ_c[2], ξ_c[1]]
    ιJξ = [to_Spectral(flux_compose(ι, ψ_c) * Jξi) for Jξi in Jξ]

    res = ϕ_deriv.(ξ_s) + to_Spectral.(Gh) - ιJξ

    ξinv_s = invert_coordinates(nae.ξ_s; M=Mc)
    ξinvtree = composition_basis(ξinv_s; M=Mc)
    
    # hξ = ϕ_deriv.(ξ_s) + to_Spectral.(Gh);
    # [compose(hξi, ξinvtree) for hξi in hξ]
    
    [compose(resi, ξinvtree) for resi in res]
end

function get_flux_coordinates(nae::DirectNearAxisEquilibrium) 
    ## Setup   
    Ms = nae.Ms;
    Mc = nae.Mc;
    Nρ = nae.Nρ;
    
    to_Spectral = (x) -> SpectralPowerSeries(x, M=nae.Ms);
    to_Spatial = (x) -> SpatialPowerSeries(x, M=nae.Mc);

    BK = nae.BK_c
    # BK = nae.ginv_c*nae.B_c
    hρ = BK[1] / BK[3];
    hθ = BK[2] / BK[3];
    
    x_s = zero_SpectralPowerSeries(Ms, Nρ);
    x_s[2].a[1,1] = 1.;
    x_c = to_Spatial(x_s);
    y_s = similar(x_s);
    y_s[2].a[1,2] = 1.
    y_c = to_Spatial(y_s);
    ρ = PowerSeries_ρ();

    identity_s = [x_s, y_s]

    hx_c = remove_zeros(*(x_c,hρ,N=Nρ+1)*inv(ρ) - *(y_c,hθ,N=Nρ+1))
    hy_c = remove_zeros(*(y_c,hρ,N=Nρ+1)*inv(ρ) + *(x_c,hθ,N=Nρ+1))

    hx_s = to_Spectral(hx_c);
    hy_s = to_Spectral(hy_c);
    hx_c = to_Spatial(hx_s)
    hy_c = to_Spatial(hy_s)

    for ii = 1:2
        for jj = 1:Nρ
            nae.ξ_s[ii][jj].a[:] .= 0.;
            nae.ξ_c[ii][jj].a[:] .= 0.;
            nae.ι[jj] = 0.
        end
    end


    ## Get the leading order behavior
    flux_leading_order(hx_c, hy_c, nae)

    res = flux_resid(hx_c, hy_c, nae)
    # println("n=2")
    # for jj = 1:Nρ
    #     for ii = 1:2
    #         println("norm(res[$ii][$jj]) = $(norm(res[ii][jj].a))")
    #     end
    # end
    
    ## Get higher order behavior
    ι0 = nae.ι[1]
    Dϕ = ϕ_deriv_operator(Ms)
    J = [0. -1.; 1. 0.]
    minus_ι0J = -ι0*J
    I2 = [1. 0.; 0. 1.]
    for n = 3:Nρ
        Dθ = θ_deriv_operator(n)
        minus_resn = zeros(2, Ms, n); # index resn[ξdim, ϕmode, θmode]
        minus_resn[1, :, :] = -res[1][n].a
        minus_resn[2, :, :] = -res[2][n].a
        ξupdate = similar(minus_resn)
        if isodd(n)
            # Do the O(1) stuff
            # jj = 1
            ξupdate[:,1,1] = minus_ι0J\minus_resn[:,1,1];

            for jj = 2:2:Ms
                ind = jj:jj+1
                Dϕjj = Dϕ[ind,ind]
                op = kron(Dϕjj, I2) + kron(I2, minus_ι0J)
                ξupdate[:,ind,1] = reshape(op\vec(minus_resn[:,ind,1]),2,2)
            end
        end

        odd = mod(n, 2);
        
        for ii = 1:n÷2
            m = (1 + odd + 2*(ii-1))
            indi = m:m+1
            Dθkk = ι0 .* Dθ[indi, indi]
            
            # Constant in ϕ contribution
            if (ii == 1) && iseven(n)
                P = [1. 0. 0. 1. ; 0. 1. -1. 0.; 0. 1. 1. 0.; 1. 0. 0. -1.]' ./ sqrt(2)
                Λ = zeros(4,4); Λ[3:4, 3:4] = 2ι0 .* J;
                resii = P'*vec(minus_resn[:,1,indi])
                nae.ι[n-1] = -resii[2]/sqrt(2)
                ξupdate[:,1,indi] = reshape(P[:,3:4]*((2ι0 .* J)\resii[3:4]),2,2)
                
                # println("$(n-1), ι[$(n-1)] = $(nae.ι[n-1])")
                # println("minus_resn[:,1,indi] = $(norm(minus_resn[:,1,indi]))")
                # println("resii = $(norm(resii))")
                # op = kron(Dθkk, I2) + kron(I2, minus_ι0J)
                # op = P*Λ*P'
                # println("After Residual = $(norm(op*vec(ξupdate[:,1,indi]) - vec(minus_resn[:,1,indi])-nae.ι[n-1]*vec(J) ))")
            else
                op = kron(Dθkk, I2) + kron(I2, minus_ι0J)
                ξupdate[:,1,indi] = reshape(op\vec(minus_resn[:,1,indi]),2,2)
            end
            
            for jj = 2:2:Ms
                indj = jj:jj+1
                Dϕjj = Dϕ[indj, indj]
                op = kron(Dθkk, I2, I2) + kron(I2, Dϕjj, I2) + kron(I2, I2, minus_ι0J)
                ξupdate[:,indj,indi] = reshape(op\vec(minus_resn[:,indj,indi]),2,2,2)
            end
        end
        
        ξmap_n = deepcopy(identity_s)
        ξmap_n[1][n] = ξupdate[1,:,:]
        ξmap_n[2][n] = ξupdate[2,:,:]
        # println("typeof(ξmap_n[1]) = $(typeof(ξmap_n[1]))")
        
        ξ_tree = composition_basis(nae.ξ_s)
        nae.ξ_s[:] = [compose(ξmap_ni, ξ_tree) for ξmap_ni in ξmap_n]
        nae.ξ_c[:] = to_Spatial.(nae.ξ_s)
        ψ_c= dot(nae.ξ_c,nae.ξ_c);
        ψ_s = to_Spectral(ψ_c)
        for ii = 1:get_N(ψ_c)
            nae.ψ_c[ii] = ψ_c[ii]
            nae.ψ_s[ii] = ψ_s[ii]
        end

        res = flux_resid(hx_c, hy_c, nae)
        # println("n=$n")
        # for jj = 1:Nρ
        #     println("norm(res[:][$jj]) = $(sum([norm(res[ii][jj].a) for ii = 1:2]))")
        # end
    end
end