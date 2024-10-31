using Documenter
using StellaratorNearAxis
using Pkg; Pkg.add("CairoMakie")
using CairoMakie

makedocs(
    sitename = "StellaratorNearAxis",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold = 10240000
    ),
    modules = [StellaratorNearAxis],
    pages = [
        "StellaratorNearAxis.jl" => "index.md"
        # "Examples" => [
        #     "Birkhoff Averages" => "examples/birkhoff_averaging/birkhoff_averaging.md",
        #     "Birkhoff Extrapolation" => "examples/extrapolation/extrapolation.md",
        #     "Approximately Invariant Kernel Functions" => "examples/kernel/kernel.md"
        # ],
        "Documentation" => [
            "Power Series" => "lib/PowerSeries.md"
        ]
        # "Internal Documentation" => "lib/Internal.md"
    ],
    warnonly=true
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/maxeruth/StellaratorNearAxis.jl.git"
)

