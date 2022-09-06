# Double Mach reflection of a strong shock

This folder provides setups and results for the classical double Mach
reflection problem of Woodward and Colella for the compressible Euler equations
with Trixi.jl. To reproduce the results, you need to execute

```bash
julia -e 'include("double_mach_reflection.jl"); run_experiments()'
```

You can also add some flags to Julia to increase the runtime performance, e.g.,

```bash
julia --check-bounds=no --threads=4 -e 'include("double_mach_reflection.jl"); run_experiments()'
```

Within Julia, you can verify that CFL-based step size control needs very small
CFL factors to run the first few time step without crashing by executing

```julia
julia> include("double_mach_reflection.jl")

julia> run_double_mach_reflection(SSPRK43, adaptive=false, cfl=0.01, maxiters=100); # crashes

julia> run_double_mach_reflection(SSPRK43, adaptive=false, cfl=0.001, maxiters=100); # runs
```
