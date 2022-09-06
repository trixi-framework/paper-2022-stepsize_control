# Astrophysical jet

This folder provides setups and results for an astrophysical Mach 2000 jet
problem for the compressible Euler equations with Trixi.jl. To reproduce the
results, you need to execute

```bash
julia -e 'include("astro_jet.jl"); run_experiments()'
```

You can also add some flags to Julia to increase the runtime performance, e.g.,

```bash
julia --check-bounds=no --threads=4 -e 'include("astro_jet.jl"); run_experiments()'
```

Within Julia, you can verify that CFL-based step size control needs very small
CFL factors to run the first few time step without crashing by executing

```julia
julia> include("astro_jet.jl")

julia> run_astro_jet(SSPRK43, adaptive=false, cfl=1.0e-8, maxiters=100); # crashes

julia> run_astro_jet(SSPRK43, adaptive=false, cfl=1.0e-9, maxiters=100); # runs
```
