using Literate
using Base
execute = true
# ENV["JULIA_DEBUG"]="Literate"
Literate.notebook("./example.jl", pwd()*"/../examples/"; execute)
Literate.markdown("./example.jl", pwd()*"/../docs/src/examples/"; execute)
