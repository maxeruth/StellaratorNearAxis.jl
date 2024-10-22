function *(::IdentityPowerSeries, B::AbstractPowerSeries; kwargs...); B; end
function *(A::AbstractPowerSeries, ::IdentityPowerSeries; kwargs...); A; end
function *(::IdentityPowerSeries, B::ZeroPowerSeries; kwargs...); B; end
function *(A::ZeroPowerSeries, ::IdentityPowerSeries; kwargs...); A; end
# function get_p0(A::IdentityPowerSeries); 0; end
# function get_N(A::IdentityPowerSeries); Inf; end
function similar(A::IdentityPowerSeries{T}; optargs...) where {T}
    IdentityPowerSeries(;T)
end
function get_p0(A::IdentityPowerSeries{T}) where {T}; zero(T); end
function change_order(A::IdentityPowerSeries{T}, N::Integer) where {T}
    IdentityPowerSeries(;T)
end
function evaluate(A::IdentityPowerSeries{T}, ρ, θ, s) where {T}; ones(T,length(ρ),length(θ),length(s)); end