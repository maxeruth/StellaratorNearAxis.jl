using Documenter
using StellaratorNearAxis
using Pkg; Pkg.add("CairoMakie")
using CairoMakie
using LinearAlgebra

makedocs(
    sitename = "StellaratorNearAxis",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold = 10240000
    ),
    modules = [StellaratorNearAxis],
    pages = [
        "StellaratorNearAxis.jl" => "index.md",
        "Documentation" => "lib/PowerSeries.md",
        "Example" => "examples/example.md"
    ],
    # warnonly=true
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/maxeruth/StellaratorNearAxis.jl.git"
)

