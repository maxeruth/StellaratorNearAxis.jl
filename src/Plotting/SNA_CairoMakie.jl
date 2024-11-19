# """
#     surface_heatmap(r::Vector, Q::AbstractArray, u::SpectralPowerSeries, 
#                     ρ::Number, Nθ::Integer, Ns::Integer)
# """
# function surface_heatmap(ax, r::Vector, Q::AbstractArray,
#                          u::SpectralPowerSeries, ρ::Number, Nθ::Integer, Ns::Integer;
#                          colormap=:diverging_isoluminant_cjo_70_c25_n256)
#     θs = (0:Nθ-1).*(2π/Nθ);
#     ss = (0:Ns-1).*(2π/Ns);

#     f_r = (x)->begin
#         ρ,θ,s=x
#         r0 = [evaluate(ri,ρ,θ,s)[1] for ri in r]
#         Q0 = [evaluate(qi,ρ,θ,s)[1] for qi in Q]
#         r0 + Q0*[ρ*cos(θ),ρ*sin(θ),0]
#     end

#     rs = zeros(3,Nθ,Ns);
#     for (ii,θ) in enumerate(θs), (jj,s) = enumerate(ss)
#         rs[:,ii,jj] = f_r([ρ,θ,s])
#     end
#     u = evaluate(u,ρ,θs,ss)[1,:,:];

#     faces = zeros(Integer, 3, 2, Nθ, Ns)
#     for ii = 1:Nθ, jj = 1:Ns
#         iip1 = mod1(ii+1,Nθ);
#         jjp1 = mod1(jj+1,Ns);
#         faces[:,1,ii,jj] = [ii+Nθ*(jj-1), iip1+Nθ*(jjp1-1), ii+Nθ*(jjp1-1)]
#         faces[:,2,ii,jj] = [iip1+Nθ*(jjp1-1), ii+Nθ*(jj-1), iip1+Nθ*(jj-1)]
#         # faces[:,1,ii,jj] = [1,1,1]
#     end
#     faces = reshape(faces,3,2*Nθ*Ns)';
#     rs = reshape(rs,3,Nθ*Ns)';


#     cr = cgrad(colormap)
#     Ncr = length(cr);
#     u = u .- minimum(u);
#     u = u./maximum(u);
#     u = round.(Integer, u .* (Ncr-1)) .+ 1

#     mesh!(ax, rs, faces, color=cr[u[:]])
# end


# """
#     plot_curve!(ax, r0::Vector, M::Number; color=:black)

# Plot a 3D curve from a vector of SpectralPowerSeries `r0`. Uses `M`
# mesh points.
# """
# function plot_curve!(ax, r0::Vector, M::Integer; kwargs...)
#     r0_c = [SpatialPowerSeries(ri; M=M) for ri in r0]
#     f_periodic = (x) -> vcat(x, x[1])

#     x = f_periodic(r0_c[1][1].a)
#     y = f_periodic(r0_c[2][1].a)
#     z = f_periodic(r0_c[3][1].a)

#     lines!(x,y,z; kwargs...)
# end