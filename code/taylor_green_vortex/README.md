# Inviscid Taylor-Green vortex

This folder provides setups and results of the classical inviscid Taylor-Green
vortex for the compressible Euler equations of an ideal gas with Trixi.jl.
To reproduce the results, you need to execute

```bash
julia -e 'include("taylor_green_vortex.jl"); run_experiments()'
```

You can also add some flags to Julia to increase the runtime performance, e.g.,

```bash
julia --check-bounds=no --threads=4 -e 'include("taylor_green_vortex.jl"); run_experiments()'
```

This will also generate some output files in the directory `out_StructuredMesh`
that you can visualize, e.g., with Paraview. The Paraview state
`taylor_green_vortex_StructuredMesh_pressure.pvsm` contains all settings to
reproduce the plot shown in the paper.
