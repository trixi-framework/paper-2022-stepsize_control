# Kelvin-Helmholtz instability

This folder provides setups and results of a Kelvin-Helmholtz instability
problem for the compressible Euler equations with Trixi.jl. To reproduce the
results, you need to execute

```bash
julia -e 'include("kelvin_helmholtz.jl"); run_experiments()'
```

You can also add some flags to Julia to increase the runtime performance, e.g.,

```bash
julia --check-bounds=no --threads=4 -e 'include("kelvin_helmholtz.jl"); run_experiments()'
```
