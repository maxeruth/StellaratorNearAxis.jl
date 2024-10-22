## ZeroPowerSeries Functions
function +(A::ZeroPowerSeries, B::AbstractPowerSeries; kwargs...); B; end
function +(A::AbstractPowerSeries, B::ZeroPowerSeries; kwargs...); A; end
function +(A::ZeroPowerSeries, B::ZeroPowerSeries; kwargs...); B; end
function -(A::ZeroPowerSeries, B::AbstractPowerSeries; kwargs...); B; end
function *(A::ZeroPowerSeries, B::AbstractPowerSeries; kwargs...); ZeroPowerSeries(); end
function *(A::AbstractPowerSeries, B::ZeroPowerSeries; kwargs...); ZeroPowerSeries(); end
function *(A::ZeroPowerSeries, B::ZeroPowerSeries; kwargs...); ZeroPowerSeries(); end
# function get_p0(A::ZeroPowerSeries); 0; end
# function get_N(A::ZeroPowerSeries); Inf; end
function similar(A::ZeroPowerSeries{T}; optargs...) where {T}
    ZeroPowerSeries(;T)
end
function zero(A::AbstractPowerSeries{T}) where {T}; ZeroPowerSeries(;T); end
function Fnorm(A::ZeroPowerSeries{T}) where {T}; zero(T); end
function get_p0(A::ZeroPowerSeries); 0; end
function change_order(A::ZeroPowerSeries{T}, N::Integer) where {T}
    ZeroPowerSeries(;T)
end
function evaluate(A::ZeroPowerSeries{T}, ρ, θ, s) where {T}; zeros(T,length(ρ),length(θ),length(s)); end
