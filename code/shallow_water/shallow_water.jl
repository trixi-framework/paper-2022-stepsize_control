# Shallow water equations with moving bottom topography

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using DelimitedFiles
using LinearAlgebra: norm

using OrdinaryDiffEq
using Trixi


include("../setup_pyplot.jl")
include("../monitor_stepsize.jl")


@doc raw"""
    ShallowWaterExnerEquations2D(gravity, H0, xi, Ag)

Shallow water Exner equations in two space dimensions given by
```math
\begin{aligned}
  \frac{\partial h}{\partial t} + \frac{\partial}{\partial x}(h v_1)
    + \frac{\partial}{\partial y}(h v_2) &= 0 \\
    \frac{\partial}{\partial t}(h v_1) + \frac{\partial}{\partial x}\left(h v_1^2 + \frac{g}{2}h^2\right)
    + \frac{\partial}{\partial y}(h v_1 v_2) + g h \frac{\partial b}{\partial x} &= 0 \\
    \frac{\partial}{\partial t}(h v_2) + \frac{\partial}{\partial x}(h v_1 v_2)
    + \frac{\partial}{\partial y}\left(h v_2^2 + \frac{g}{2}h^2\right) + g h \frac{\partial b}{\partial y} &= 0 \\
    \frac{\partial b}{\partial t} + \xi\frac{\partial}{\partial x}(q_{b1}) + \xi\frac{\partial}{\partial y}(q_{b2}) &= 0.
\end{aligned}
```
The Exner equation for the time evolution of the bottom topography contains the flux of the bottom sediment
``\mathbf{q}_b`` with the constant ``\xi = 1/(1-\simga), \sigma\in(0,1)`` being the porosity of the bed material.
There exist many closure models to determine the form of the sediment flux ``\mathbf{q}_b``.
Herein, a simple model due to Grass (1981) is considered:
```math
\mathbf{q}_b = A_g \mathbf{v} (v_1^2+v_2^2)^{\frac{m-1}{2}}.
```
So, the bottom flux is a function of the flow velocities and a constant coefficient ``A_g\in(0,1)`` and ``1\leq m \leq 4``.
This implementation makes a further assumption that ``m = 3`` such that the Grass model becomes
```math
\mathbf{q}_b = A_g \mathbf{v} (v_1^2+v_2^2).
```

The unknown quantities are the water height ``h`` and the velocities ``\mathbf{v} = (v_1, v_2)^T``.
The gravitational constant is denoted by `g` and the (possibly) variable bottom topography function ``b(x,y,t)``.
Conservative variable water height ``h`` is measured from the bottom topography ``b``, therefore one
also defines the total water height as ``H = h + b``.

The additional quantity ``H_0`` is also available to store a reference value for the total water height that
is useful to set initial conditions or test the "lake-at-rest" well-balancedness.

The bottom topography function ``b(x,y,t)`` is set inside the initial condition routine
for a particular problem setup.

Two good starting references for the shallow water Exner model:
- Justin Hudson and Peter K. Sweby (2001)
  Formulations for Numerically Approximating Hyperbolic Systems Governing Sediment Transport
  [DOI: 10.1023/A:1025304008907](https://doi.org/10.1023/A:1025304008907)
- Fayssal Benkhaldoun, Slah Sahmim, and Mohammed Seaid
  A two-dimensional finite volume morphodynamic model on unstructured triangular grids
  [DOI: 10.1002/fld.2129](https://doi.org/10.1002/fld.2129)
"""
struct ShallowWaterExnerEquations2D{RealT<:Real} <: Trixi.AbstractEquations{2, 4} # NDIMS, NVARS
  gravity::RealT # gravitational constant
  H0::RealT      # constant "lake-at-rest" total water height
  xi::RealT      # 1 / (1 - porosity) scaling factor on the sediment discharge
  Ag::RealT      # constant in [0,1]; strength of interaction between fluid and sediment
end

function ShallowWaterExnerEquations2D(; gravity_constant, H0=0.0, porosity=0.0, Ag)
  ShallowWaterExnerEquations2D(gravity_constant, H0, inv(1.0 - porosity), Ag)
end

Trixi.have_nonconservative_terms(::ShallowWaterExnerEquations2D) = Val(true)

Trixi.varnames(::typeof(cons2cons), ::ShallowWaterExnerEquations2D) = ("h", "h_v1", "h_v2", "b")
# Note, we use the total water height, H = h + b, as the first primitive variable for easier
# visualization and setting initial conditions
Trixi.varnames(::typeof(cons2prim), ::ShallowWaterExnerEquations2D) = ("H", "v1", "v2", "b")


# Pointwise fluxes and numerical two-point fluxes
@inline function Trixi.flux(u, normal_direction::AbstractVector,
                            equations::ShallowWaterExnerEquations2D)
  h = waterheight(u, equations)
  v1, v2 = velocity(u, equations)

  v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
  h_v_normal = h * v_normal
  p = 0.5 * equations.gravity * h^2

  f1 = h_v_normal
  f2 = h_v_normal * v1 + p * normal_direction[1]
  f3 = h_v_normal * v2 + p * normal_direction[2]
  f4 = equations.xi * equations.Ag * v_normal * (v1^2 + v2^2)
  return SVector(f1, f2, f3, f4)
end

@inline function Trixi.flux_wintermeyer_etal(u_ll, u_rr,
                                             normal_direction::AbstractVector,
                                             equations::ShallowWaterExnerEquations2D)
  # Unpack left and right state
  h_ll, h_v1_ll, h_v2_ll, _ = u_ll
  h_rr, h_v1_rr, h_v2_rr, _ = u_rr

  # Get the velocities on either side
  v1_ll, v2_ll = velocity(u_ll, equations)
  v1_rr, v2_rr = velocity(u_rr, equations)

  # Average each factor of products in flux
  h_v1_avg = 0.5 * (h_v1_ll + h_v1_rr )
  h_v2_avg = 0.5 * (h_v2_ll + h_v2_rr )
  v1_avg   = 0.5 * (v1_ll   + v1_rr   )
  v2_avg   = 0.5 * (v2_ll   + v2_rr   )
  p_avg    = 0.5 * equations.gravity * h_ll * h_rr

  # Calculate fluxes depending on normal_direction
  f1 = h_v1_avg * normal_direction[1] + h_v2_avg * normal_direction[2]
  f2 = f1 * v1_avg + p_avg * normal_direction[1]
  f3 = f1 * v2_avg + p_avg * normal_direction[2]
  # the splitting the f4 flux is ad hoc currently
  v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
  v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
  # v_dot_n_avg = 0.5 * (v_dot_n_ll + v_dot_n_rr)
  # f4 = equations.xi * equations.Ag * v_dot_n_avg * (v1_avg^2 + v2_avg^2)
  f4_ll = equations.xi * equations.Ag * v_dot_n_ll * (v1_ll^2 + v2_ll^2)
  f4_rr = equations.xi * equations.Ag * v_dot_n_rr * (v1_rr^2 + v2_rr^2)
  f4 = 0.5 * (f4_ll + f4_rr)

  return SVector(f1, f2, f3, f4)
end

@inline function Trixi.flux_nonconservative_wintermeyer_etal(
    u_ll, u_rr,
    normal_direction_ll::AbstractVector,
    normal_direction_average::AbstractVector,
    equations::ShallowWaterExnerEquations2D)
  # Pull the necessary left and right state information
  h_ll = waterheight(u_ll, equations)
  b_rr = u_rr[4]
  # Note this routine only uses the `normal_direction_average` and the average of the
  # bottom topography to get a quadratic split form DG gradient on curved elements
  return SVector(zero(eltype(u_ll)),
                 normal_direction_average[1] * equations.gravity * h_ll * b_rr,
                 normal_direction_average[2] * equations.gravity * h_ll * b_rr,
                 zero(eltype(u_ll)))
end

@inline function Trixi.flux_nonconservative_fjordholm_etal(
    u_ll, u_rr,
    normal_direction_ll::AbstractVector,
    normal_direction_average::AbstractVector,
    equations::ShallowWaterExnerEquations2D)
  # Pull the necessary left and right state information
  h_ll, _, _, b_ll = u_ll
  h_rr, _, _, b_rr = u_rr

  # Comes in two parts:
  #   (i)  Diagonal (consistent) term from the volume flux that uses `normal_direction_average`
  #        but we use `b_ll` to avoid cross-averaging across a discontinuous bottom topography
  f2 = normal_direction_average[1] * equations.gravity * h_ll * b_ll
  f3 = normal_direction_average[2] * equations.gravity * h_ll * b_ll

  #   (ii) True surface part that uses `normal_direction_ll`, `h_average` and `b_jump`
  #        to handle discontinuous bathymetry
  h_average = 0.5 * (h_ll + h_rr)
  b_jump = b_rr - b_ll

  f2 += normal_direction_ll[1] * equations.gravity * h_average * b_jump
  f3 += normal_direction_ll[2] * equations.gravity * h_average * b_jump

  # First and last equations do not have a nonconservative flux
  f1 = f4 = zero(eltype(u_ll))

  return SVector(f1, f2, f3, f4)
end


@inline function Trixi.max_abs_speeds(u, equations::ShallowWaterExnerEquations2D)
  h = waterheight(u, equations)
  v1, v2 = velocity(u, equations)

  c = equations.gravity * sqrt(h)
  return abs(v1) + c, abs(v2) + c
end


# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                           equations::ShallowWaterExnerEquations2D)
  # Extract and compute the velocities in the normal direction
  v1_ll, v2_ll = velocity(u_ll, equations)
  v1_rr, v2_rr = velocity(u_rr, equations)
  v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
  v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

  # Compute the wave celerity on the left and right
  h_ll = waterheight(u_ll, equations)
  h_rr = waterheight(u_rr, equations)
  c_ll = sqrt(equations.gravity * h_ll)
  c_rr = sqrt(equations.gravity * h_rr)

  # The normal velocities are already scaled by the norm
  return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end


# Helper function to extract variables from the conservative variables
@inline function Trixi.velocity(u, equations::ShallowWaterExnerEquations2D)
  h, h_v1, h_v2, _ = u

  v1 = h_v1 / h
  v2 = h_v2 / h
  return SVector(v1, v2)
end

@inline function waterheight(u, equations::ShallowWaterExnerEquations2D)
  return u[1]
end


# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::ShallowWaterExnerEquations2D)
  h, _, _, b = u

  H = h + b
  v1, v2 = velocity(u, equations)
  return SVector(H, v1, v2, b)
end


# Entropy function for the shallow water equations is the total energy
@inline Trixi.entropy(cons, equations::ShallowWaterExnerEquations2D) = energy_total(cons, equations)

@inline function Trixi.energy_total(cons, equations::ShallowWaterExnerEquations2D)
  h, h_v1, h_v2, b = cons

  e = (h_v1^2 + h_v2^2) / (2 * h) + 0.5 * equations.gravity * h^2 + equations.gravity * h * b
  return e
end

@inline function Trixi.cons2entropy(u, equations::ShallowWaterExnerEquations2D)
  h, _, _, b = u

  v1, v2 = velocity(u, equations)
  v_square = v1^2 + v2^2

  w1 = equations.gravity * (h + b) - 0.5 * v_square
  w2 = v1
  w3 = v2
  w4 = equations.gravity * h
  return SVector(w1, w2, w3, w4)
end



# Example taken from Section 4.3 of Benkhaldoun et al. (2010)
# [DOI: 10.1002/fld.2129](https://doi.org/10.1002/fld.2129)
function initial_condition_cone_island(x, t, equations::ShallowWaterExnerEquations2D)
  x1, x2 = x
  b = 0.0
  if 300 <= x1 <= 500 && 400 <= x2 <= 600
    b += (sin(pi * (x1 - 300) / 200))^2 * (sin(pi * (x2 - 400) / 200))^2
  end

  # Set the background values
  h = equations.H0 - b
  h_v1 = 10.0
  h_v2 = 0.0

  return SVector(h, h_v1, h_v2, b)
end


function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector, x, t,
                                      surface_flux_function, equations::ShallowWaterExnerEquations2D)
  # normalize the outward pointing direction
  normal = normal_direction / norm(normal_direction)

  # compute the normal velocity
  u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

  # create the "external" boundary solution state
  u_boundary = SVector(u_inner[1],
                       u_inner[2] - 2.0 * u_normal * normal[1],
                       u_inner[3] - 2.0 * u_normal * normal[2],
                       u_inner[4])

  # calculate the boundary flux
  flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)

  return flux
end


function boundary_condition_subcritical_inflow(u_inner, normal_direction::AbstractVector, x, t,
                                               surface_flux_function, equations::ShallowWaterExnerEquations2D)
  # Impulse and bottom from outside, height from inside
  u_outer = SVector(u_inner[1], 10.0, 0.0, 0.0)
  # calculate the boundary flux
  flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)

  return flux
end


function boundary_condition_subcritical_outflow(u_inner, normal_direction::AbstractVector, x, t,
                                                surface_flux_function, equations::ShallowWaterExnerEquations2D)
  # Impulse from inside, height and bottom from outside
  u_outer = SVector(equations.H0, u_inner[2], u_inner[3], 0.0)
  # calculate the boundary flux
  flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)

  return flux
end


# Particular setup for the conical sand dune run useful for comparing different time stepping algorithms
function run_shallow_water_exner_timings(alg; cfl=nothing, kwargs...)

  equations = ShallowWaterExnerEquations2D(gravity_constant=9.8, H0=10.0,
                                           porosity=0.4, Ag=0.001)

  volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
  solver = DGSEM(polydeg=5,
                 surface_flux=(flux_lax_friedrichs, flux_nonconservative_fjordholm_etal),
                 volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

  mesh_file = joinpath(@__DIR__, "box_mesh.mesh")
  mesh = UnstructuredMesh2D(mesh_file)

  initial_condition = initial_condition_cone_island

  boundary_conditions = Dict( :Bottom => boundary_condition_slip_wall,
                              :Top    => boundary_condition_slip_wall,
                              :Right  => boundary_condition_subcritical_outflow,
                              :Left   => boundary_condition_subcritical_inflow )

  semi = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver; boundary_conditions)

  # final time is for 100 hrs = 100 hours * 60 minutes/hour * 60 seconds/minute
  tspan = (0.0, 360_000.0)
  ode = semidiscretize(semi, tspan)

  summary_callback = SummaryCallback()

  analysis_interval = 2000
  analysis_callback = AnalysisCallback(semi, interval=analysis_interval)
  alive_callback = AliveCallback(analysis_interval=analysis_interval)

  basic_callbacks = (summary_callback, analysis_callback, alive_callback)

  if cfl !== nothing
    stepsize_callback = StepsizeCallback(cfl=cfl)
    basic_callbacks = (basic_callbacks..., stepsize_callback)
  end
  times, effective_cfl_numbers, monitor_stepsize_callback = monitor_stepsize()
  callbacks = CallbackSet(basic_callbacks..., monitor_stepsize_callback)

  sol = solve(ode, alg, save_everystep=false, callback=callbacks; kwargs...)
  summary_callback()

  return sol, times, effective_cfl_numbers
end


# Main functions for this setup
function run_experiments()

  let alg = BS3()
    # Controller suggested by Ranocha, Dalcin, Parsani, Ketcheson (2021)
    controller = PIDController(0.6, -0.2)
    name = string(nameof(typeof(alg)))

    sol_err, times_err, effective_cfl_numbers = run_shallow_water_exner_timings(
      alg; adaptive=true, controller,
      abstol=1.0e-7, reltol=1.0e-7,
      maxiters=10^8)

    @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
    open(joinpath(@__DIR__, "shallow_water_exner_$(name)_err.txt"), "w") do io
      println(io, "# nf = ",      sol_err.destats.nf)
      println(io, "# naccept = ", sol_err.destats.naccept)
      println(io, "# nreject = ", sol_err.destats.nreject)
      println(io, "# t cfl")
      writedlm(io, hcat(times_err, effective_cfl_numbers))
    end
  end

  let alg = RDPK3SpFSAL35()
    # Controller suggested by Ranocha, Dalcin, Parsani, Ketcheson (2021) is
    # already set as default in OrdinaryDiffEq.jl
    name = string(nameof(typeof(alg)))

    sol_err, times_err, effective_cfl_numbers = run_shallow_water_exner_timings(
      alg; adaptive=true,
      abstol=1.0e-7, reltol=1.0e-7,
      maxiters=10^8)

    @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
    open(joinpath(@__DIR__, "shallow_water_exner_$(name)_err.txt"), "w") do io
      println(io, "# nf = ",      sol_err.destats.nf)
      println(io, "# naccept = ", sol_err.destats.naccept)
      println(io, "# nreject = ", sol_err.destats.nreject)
      println(io, "# t cfl")
      writedlm(io, hcat(times_err, effective_cfl_numbers))
    end
  end

  let alg = RDPK3SpFSAL49()
    # Controller suggested by Ranocha, Dalcin, Parsani, Ketcheson (2021) is
    # already set as default in ordinaryDiffEq.jl
    name = string(nameof(typeof(alg)))

    sol_err, times_err, effective_cfl_numbers = run_shallow_water_exner_timings(
      alg; adaptive=true,
      abstol=1.0e-7, reltol=1.0e-7,
      maxiters=10^8)

    @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
    open(joinpath(@__DIR__, "shallow_water_exner_$(name)_err.txt"), "w") do io
      println(io, "# nf = ",      sol_err.destats.nf)
      println(io, "# naccept = ", sol_err.destats.naccept)
      println(io, "# nreject = ", sol_err.destats.nreject)
      println(io, "# t cfl")
      writedlm(io, hcat(times_err, effective_cfl_numbers))
    end
  end
end


# Particular setup for the conical sand dune run useful for creating figures
function run_shallow_water_exner()

  equations = ShallowWaterExnerEquations2D(gravity_constant=9.8, H0=10.0,
                                           porosity=0.4, Ag=0.001)

  volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
  solver = DGSEM(polydeg=5,
                 surface_flux=(flux_lax_friedrichs, flux_nonconservative_fjordholm_etal),
                 volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

  mesh_file = joinpath(@__DIR__, "box_mesh.mesh")
  mesh = UnstructuredMesh2D(mesh_file)

  initial_condition = initial_condition_cone_island

  boundary_conditions = Dict( :Bottom => boundary_condition_slip_wall,
                              :Top    => boundary_condition_slip_wall,
                              :Right  => boundary_condition_subcritical_outflow,
                              :Left   => boundary_condition_subcritical_inflow )

  semi = SemidiscretizationHyperbolic(
    mesh, equations, initial_condition, solver; boundary_conditions)

  # final time is for 100 hrs = 100 hours * 60 minutes/hour * 60 seconds/minute
  tspan = (0.0, 360_000.0)
  ode = semidiscretize(semi, tspan)

  summary_callback = SummaryCallback()

  analysis_interval = 2000
  analysis_callback = AnalysisCallback(semi, interval=analysis_interval)
  alive_callback = AliveCallback(analysis_interval=analysis_interval)

  save_solution = SaveSolutionCallback(interval=5000,
                                      save_initial_solution=true,
                                      save_final_solution=true)

  callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution)

  sol = solve(ode, RDPK3SpFSAL49(), abstol=1.0e-7, reltol=1.0e-7,
              save_everystep=false, callback=callbacks, adaptive=true,
              maxiters=10^7)
  summary_callback()
end
