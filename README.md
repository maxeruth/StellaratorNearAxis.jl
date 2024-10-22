# StellaratorNearAxis.jl
This is a package that computes the stellarator near-axis expansion to high order. It will be in flux over the next few weeks (written 10/22/2024). On the to-do list are:
1. Create a more comprehensive test suite and the corresponding CI and code coverage bots.
2. Incorporate `Documenter.jl` and make a more general documentation.
3. Include a couple example notebooks.
4. Creation of a couple quality-of-life macros for automatically checking series operations
5. Transition to a greek-free typography (for those who don't know how to type $\rho$ in Julia)

On the longer-term to-do list are:
1. Implementation of more metrics (such as quasisymmetry)
2. Incorporation of other solvers besides direct vacuum near-axis
3. Addition of AD tools for optimization
4. Hooks to Simsopt and/or StellaratorOptimization.jl

If you would like to use this package in your research and would particularly like one of these features, contact me at `maximilian.ruth@austin.utexas.edu`. I'm happy to collaborate!