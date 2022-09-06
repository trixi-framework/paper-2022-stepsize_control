# Orszag-Tang vortex

This folder provides setups and results of the classical Orszag-Tang vortex
problem for the ideal GLM MHD equations with Trixi.jl. To reproduce the results,
you need to execute

```bash
julia -e 'include("orszag_tang.jl"); run_experiments()'
```

You can also add some flags to Julia to increase the runtime performance, e.g.,

```bash
julia --check-bounds=no --threads=4 -e 'include("orszag_tang.jl"); run_experiments()'
```
