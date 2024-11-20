"""
    flux_surface_plot!(ax, nae::DirectNearAxisEquilibrium, u::SpectralPowerSeries, rho0::Number, 
                           theta0::Number, s0::Number; Ntheta::Integer=50, Ns::Integer=50, 
                           colormap=:diverging_isoluminant_cjo_70_c25_n256, kwargs...)

3D plot a flux surface of a near-axis equilibrium file. 
Input:
- `ax`: A 3D CairoMakie `CairoMakie.Axis3` struct
- `u`: The SpectralPowerSeries to be plotted in direct coordinates (e.g. nae.phi_s)
- `(rho0,theta0,s0)`: A spatial point the flux surface (approximately) intersects
- `(Ntheta,Ns)`: The resolution of the plot
- `(colormap,kwargs)`: Additional arguments for CairoMakie.mesh!
"""
function flux_surface_plot!(ax, nae::DirectNearAxisEquilibrium, u::SpectralPowerSeries, rho0::Number, 
                           theta0::Number, s0::Number; Ntheta::Integer=50, Ns::Integer=50, 
                           colormap=:diverging_isoluminant_cjo_70_c25_n256, kwargs...)
    flux_surface_plot!(ax, nae.r0_s,nae.Q_s,nae.xi_s,u,rho0,theta0,s0;Ntheta,Ns,colormap,kwargs...)
end

function flux_surface_plot!(ax, r0_s::AbstractVector, Q_s::AbstractMatrix, 
                           xi_s::AbstractVector, u::SpectralPowerSeries, rho0::Number, 
                           theta0::Number, s0::Number; Ntheta::Integer=50, Ns::Integer=50, 
                           colormap=:diverging_isoluminant_cjo_70_c25_n256, kwargs...)
    # Get position as a function of flux coordinates
    xi_inv = invert_coordinates(xi_s)
    xi_basis = composition_basis(xi_inv)
    Nrho = get_N(xi_inv[1])
    
    r_s = r0_to_r(r0_s, Q_s)
    r_xi = [compose(change_order(ri, Nrho),xi_basis) for ri in r_s]
    u_xi = compose(u, xi_basis);

    # Evaluate position
    R0     = norm([evaluate(xi_i,rho0,theta0,s0)[1] for xi_i in xi_s])
    Thetas = (0:Ntheta-1).*(2π/Ntheta)
    Ss     = (0:Ns    -1).*(2π/Ns    )
    r_eval = zeros(3,Ntheta,Ns)
    for ii = 1:3
        r_eval[ii,:,:] = evaluate(r_xi[ii],R0,Thetas,Ss)[1,:,:]
    end
    u_eval = evaluate(u_xi,R0,Thetas,Ss)[1,:,:]

    # Build mesh
    faces = zeros(Integer, 3, 2, Ntheta, Ns)
    for ii = 1:Ntheta, jj = 1:Ns
        iip1 = mod1(ii+1,Ntheta);
        jjp1 = mod1(jj+1,Ns);
        faces[:,1,ii,jj] = [ii+Ntheta*(jj-1), iip1+Ntheta*(jjp1-1), ii+Ntheta*(jjp1-1)]
        faces[:,2,ii,jj] = [iip1+Ntheta*(jjp1-1), ii+Ntheta*(jj-1), iip1+Ntheta*(jj-1)]
        # faces[:,1,ii,jj] = [1,1,1]
    end

    faces = reshape(faces,3,2*Ntheta*Ns)';
    r_eval = reshape(r_eval,3,Ntheta*Ns)';
    
    mesh!(ax, r_eval, faces; colormap, color=u_eval[:], kwargs...)
end

"""
    flux_contours!(ax, xi_s::AbstractVector, rho0::Number, theta0::Number, s0::Number, 
                        Nsurf::Integer; Ntheta = 101, color=:black, kwargs...)

Contour plot of flux surfaces on a cross section
Input:
- `ax`: A 3D CairoMakie `CairoMakie.Axis` struct
- `xi_s`: The SpectralPowerSeries vector defining the flux surfaces (e.g. `nae.xi_s`)
- `(rho0,theta0)`: A point on a Poincare section which the outermost surface intersects
- `s0`: The Poincare section taken
- `Nsurf`: The number of surfaces to plot
- `Ntheta`: Number of points used to plot each surface
- `(color,kwargs)`: Additional arguments for CairoMakie.lines!
"""
function flux_contours!(ax, xi_s::AbstractVector, rho0::Number, theta0::Number, s0::Number, 
                        Nsurf::Integer; Ntheta = 101, color=:black, kwargs...)    
    # Get evaluation coordinates
    Rhos = zeros(Nsurf)
    Thetas = (0:Ntheta-1) .* (2π/(Ntheta-1))
    
    for ii = 1:Nsurf
        X,Y = [evaluate(xi_si, ii*rho0/Nsurf, theta0, s0)[1] for xi_si in xi_s];
        Rhos[ii] = sqrt(X^2 + Y^2);
    end
    
    # Evaluate surface on cross section
    xs = zeros(2, Nsurf, Ntheta);
    xi_inv_s = invert_coordinates(xi_s);
    xs[1,:,:] = evaluate(xi_inv_s[1], Rhos, Thetas, s0)[:,:,1];
    xs[2,:,:] = evaluate(xi_inv_s[2], Rhos, Thetas, s0)[:,:,1];

    # Plot the surfaces
    for ii = 1:Nsurf
        lines!(ax, xs[1,ii,:], xs[2,ii,:]; color, kwargs...);
    end
end


"""
    plot_curve!(ax, r0::Vector, M::Number; color=:black)

Plot a 3D curve from a vector of SpectralPowerSeries `r0`. Uses `M` mesh points. Can be used for 
both near-axis equilibria (`r0 = nae.r0_s`) or coils (`r0 = coil.r0_s`).
"""
function plot_curve!(ax, r0::Vector, M::Integer; kwargs...)
    r0_c = [SpatialPowerSeries(ri; M=M) for ri in r0]
    f_periodic = (x) -> vcat(x, x[1])

    x = f_periodic(r0_c[1][1].a)
    y = f_periodic(r0_c[2][1].a)
    z = f_periodic(r0_c[3][1].a)

    lines!(x,y,z; kwargs...)
end