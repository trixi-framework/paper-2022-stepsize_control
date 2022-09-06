# Shallow water Exner equations with moving bottom topography

This folder provides setups and results for the shallow water equations with
a moving bottom topography with Trixi.jl. To reproduce the time step monitoring
results, you need to execute

```shell
julia -e 'include("shallow_water.jl"); run_experiments()'
```

You can also add some flags to Julia to increase the runtime performance, e.g.,

```shell
julia --check-bounds=no --threads=4 -e 'include("shallow_water.jl"); run_experiments()'
```

Note, with either command the time step monitoring runs take several hours. It will produce
three `.txt` files containing the timing and effective CFL numbers values. To reproduce the
figures present in the paper, move these `.txt` files into the folder `timing_figs` and
execute the python script `timing_figs/timing_plots.py`.

To create the solution data for the figures and angle approximation plots first execute

```shell
julia -e 'include("shallow_water.jl"); run_shallow_water_exner()'
```

Again, the performance can be increased by adding some Julia flags, e.g.,

```shell
julia --check-bounds=no --threads=4 -e 'include("shallow_water.jl"); run_shallow_water_exner()'
```

After the run is completed, the solution output files are saved in the folder `out`. These HDF5
files must be processed by the tool `Trixi2Vtk` to create `.vtu` files that can be plotted with ParaView.
For this execute in a Julia REPL in this directory

```julia
using Pkg; Pkg.activate("."); Pkg.instantiate()
using Trixi2Vtk
trixi2vtk("out/solution_000*", output_directory="out")
```

This saves the `.vtu` files into the same `out` folder that contains the raw data.
Next, launch ParaView.
To reproduce the angle comparison figure, first load the ParaView state by clicking
through `File -> Load State`
and open `paraview_figs/spread_angle_paraview_config.pvsm`. Next, from the
prompt "Load State Data File Options"
select "Choose File Names", navigate to the `out` directory and select
the `solution_000.pvd` file.

The process to recreate the star shaped bottom figure is nearly identical except one
loads the ParaView configuration state file `paraview_figs/warped_viz_paraview_config.pvsm`.

This code was developed with Julia v1.7.
