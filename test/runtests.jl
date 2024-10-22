using Test
using StellaratorNearAxis
using LinearAlgebra

# Check d/dρ ( int_ρ (A) ) = A
@testset "Test Power Series" begin
    Na = 10;
    Nb = 11
    Nc = min(Na, Nb);

    ainv = (-1 .^ (0:Nb-1)) ./ factorial.(0:Nb-1);
    Ainv = FluxPowerSeries(ainv);

    a = 1 ./ factorial.(0:Na-1);
    A = FluxPowerSeries(a); # e^x

    b = (2 .^ (0:Nb-1)) ./ factorial.(0:Nb-1);
    B = FluxPowerSeries(b); # e^2x

    c = (3 .^ (0:Nc-1)) ./ factorial.(0:Nc-1);
    c = FluxPowerSeries(b); # e^3x

    # FluxPowerSeries / Operations

    # # Times
    # @test Fnorm(c - a*b) < 1e-10
    # @test Fnorm(b - a*a) < 1e-10
    #
    # # inv
    # @test Fnorm(inv(A) - Ainv) < 1e-10
    #
    # # pow
    # @test Fnorm(b - a^2) < 1e-10

    #
    # @test

    # IdentityPowerSeries

    # ZeroPowerSeries

    # SpatialPowerSeries

    # SpectralPowerSeries

end
