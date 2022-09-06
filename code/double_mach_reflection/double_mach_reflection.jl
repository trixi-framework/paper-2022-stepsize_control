# Cold startup of the classical double Mach reflection simulation of the
# compressible Euler equations discretized in space with entropy-dissipative
# shock capturing DG methods, positivity-preserving limiters, and adaptive mesh
# refinement in Trixi.jl.

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using DelimitedFiles

using OrdinaryDiffEq
using Trixi


include("../setup_pyplot.jl")
include("../monitor_stepsize.jl")


"""
    initial_condition_double_mach_reflection(x, t, equations::CompressibleEulerEquations2D)

Compressible Euler setup for a double Mach reflection problem.
Involves strong shock interactions as well as steady / unsteady flow structures.
Also exercises special boundary conditions along the bottom of the domain that
is a mixture of Dirichlet and slip wall.
See Section IV c on the paper below for details.

- Paul Woodward and Phillip Colella (1984)
  The Numerical Simulation of Two-Dimensional Fluid Flows with Strong Shocks.
  [DOI: 10.1016/0021-9991(84)90142-6](https://doi.org/10.1016/0021-9991(84)90142-6)
"""
@inline function initial_condition_double_mach_reflection(x, t, equations::CompressibleEulerEquations2D)

  if x[1] < 1 / 6 + (x[2] + 20 * t) / sqrt(3)
    phi = pi / 6
    sin_phi, cos_phi = sincos(phi)

    rho =  8
    v1  =  8.25 * cos_phi
    v2  = -8.25 * sin_phi
    p   =  116.5
  else
    rho = 1.4
    v1  = 0
    v2  = 0
    p   = 1
  end

  prim = SVector(rho, v1, v2, p)
  return prim2cons(prim, equations)
end

# Supersonic outflow boundary condition. Solution is taken entirely from the
# internal state.
@inline function boundary_condition_outflow(
    u_inner, normal_direction::AbstractVector, x, t,
    surface_flux_function, equations::CompressibleEulerEquations2D)

  # NOTE: Only for the supersonic outflow is this strategy valid
  # Calculate the boundary flux entirely from the internal solution state
  return flux(u_inner, normal_direction, equations)
end

# Special mixed boundary condition type for the :Bottom of the domain.
# It is Dirichlet when x < 1/6 and a slip wall when x >= 1/6
@inline function boundary_condition_mixed_dirichlet_wall(
    u_inner, normal_direction::AbstractVector,
    x, t, surface_flux_function, equations::CompressibleEulerEquations2D)

  if x[1] < 1 / 6
    # From the BoundaryConditionDirichlet
    # get the external value of the solution
    u_boundary = initial_condition_double_mach_reflection(x, t, equations)
    # Calculate boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
  else # x[1] >= 1 / 6
    # Use the free slip wall BC otherwise
    flux = boundary_condition_slip_wall(u_inner, normal_direction, x, t,
              surface_flux_function, equations)
  end

  return flux
end


function run_double_mach_reflection(_alg; cfl=nothing, kwargs...)
  equations = CompressibleEulerEquations2D(1.4)

  boundary_condition_inflow = BoundaryConditionDirichlet(
    initial_condition_double_mach_reflection)
  boundary_conditions = Dict(
    :Bottom => boundary_condition_mixed_dirichlet_wall,
    :Top    => boundary_condition_inflow,
    :Right  => boundary_condition_outflow,
    :Left   => boundary_condition_inflow)

  surface_flux = flux_lax_friedrichs
  volume_flux  = flux_ranocha_turbo
  basis = LobattoLegendreBasis(4)
  indicator_sc = IndicatorHennemannGassner(equations, basis,
                                           alpha_max=0.5,
                                           alpha_min=0.001,
                                           alpha_smooth=true,
                                           variable=density_pressure)
  volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                   volume_flux_dg=volume_flux,
                                                   volume_flux_fv=surface_flux)
  solver = DGSEM(basis, surface_flux, volume_integral)

  mesh_file = joinpath(@__DIR__, "abaqus_double_mach.inp")
  mesh = P4estMesh{2}(mesh_file)

  semi = SemidiscretizationHyperbolic(mesh, equations,
    initial_condition_double_mach_reflection, solver; boundary_conditions)

  tspan = (0.0, 0.2)
  ode = semidiscretize(semi, tspan)

  summary_callback = SummaryCallback()
  analysis_interval = 500
  analysis_callback = AnalysisCallback(semi, interval=analysis_interval)
  alive_callback = AliveCallback(; analysis_interval)

  times, effective_cfl_numbers, monitor_stepsize_callback = monitor_stepsize()

  amr_indicator = IndicatorLöhner(semi, variable=density)
  amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                        base_level=0,
                                        med_level =3, med_threshold=0.05,
                                        max_level =6, max_threshold=0.1)
  amr_callback = AMRCallback(semi, amr_controller,
                            interval=2,
                            adapt_initial_condition=true,
                            adapt_initial_condition_only_refine=true)

  basic_callbacks = (summary_callback, analysis_callback, alive_callback,
                     monitor_stepsize_callback, amr_callback)

  if cfl !== nothing
    stepsize_callback = StepsizeCallback(cfl=cfl)
    basic_callbacks = (basic_callbacks..., stepsize_callback)
    kwargs = (; kwargs..., dt=1.0 #= will be overwritten by the CFL callback =#)
  end
  callbacks = CallbackSet(basic_callbacks...)

  stage_limiter! = PositivityPreservingLimiterZhangShu(
    thresholds=(1.0e-12, 1.0e-12),
    variables=(pressure, density))
  alg = _alg(stage_limiter!)

  sol = solve(ode, alg, save_everystep=false, callback=callbacks; kwargs...)
  summary_callback()

  return sol, times, effective_cfl_numbers
end



# Main function for this setup
function run_experiments()
  let alg = SSPRK43
    name = nameof(alg)
    # Controller suggested by Ranocha, Dalcin, Parsani, Ketcheson (2021)
    controller = PIDController(0.55, -0.27, 0.05)

    sol_err, times_err, effective_cfl_numbers = run_double_mach_reflection(
      alg; adaptive=true, abstol=1.0e-4, reltol=1.0e-4, controller)

    @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
    open(joinpath(@__DIR__, "double_mach_reflection_$(name)_err.txt"), "w") do io
      println(io, "# nf = ",      sol_err.destats.nf)
      println(io, "# naccept = ", sol_err.destats.naccept)
      println(io, "# nreject = ", sol_err.destats.nreject)
      println(io, "# t cfl")
      writedlm(io, hcat(times_err, effective_cfl_numbers))
    end

    plt.gcf().set_size_inches(6.3, 4.8, forward=true)
    plt.clf()
    ax1 = plt.gca()
    plt.loglog(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Time step size")
    plt.legend(loc="lower left", bbox_to_anchor=(0.1, 0.0))
    plt.xlabel("Time")
    plt.ylabel("Time step size")
    ax2 = ax1.twinx()
    plt.plot([], [])
    plt.loglog(times_err, effective_cfl_numbers, label="Effective CFL number")
    plt.legend(loc="lower right", bbox_to_anchor=(1.0, 0.13))
    plt.ylabel("Effective CFL number")
    plt.savefig(joinpath(@__DIR__, "double_mach_reflection_$(name)_time.png"), bbox_inches="tight")
    plt.savefig(joinpath(@__DIR__, "double_mach_reflection_$(name)_time.pdf"), bbox_inches="tight")

    plt.clf()
    ax1 = plt.gca()
    plt.loglog(times_err[2:end]-times_err[1:end-1], label="Time step size")
    plt.legend(loc="lower left", bbox_to_anchor=(0.1, 0.0))
    plt.xlabel("Time step")
    plt.ylabel("Time step size")
    ax2 = ax1.twinx()
    plt.plot([], [])
    plt.loglog(effective_cfl_numbers, label="Effective CFL number")
    plt.legend(loc="lower right", bbox_to_anchor=(1.0, 0.13))
    plt.ylabel("Effective CFL number")
    plt.savefig(joinpath(@__DIR__, "double_mach_reflection_$(name)_step.png"), bbox_inches="tight")
    plt.savefig(joinpath(@__DIR__, "double_mach_reflection_$(name)_step.pdf"), bbox_inches="tight")

    # hacks using internal knowledge of the plotting data structures
    pd = PlotData2D(sol_err)

    plt.clf()
    x_face, y_face = map(x->vec(vcat(x, fill(NaN, 1, size(x, 2)))), (pd.x_face, pd.y_face))
    plt.plot(x_face, y_face, color="black", linewidth=2)
    plt.xlabel(L"x")
    plt.ylabel(L"y")
    plt.xlim(extrema(pd.x))
    plt.ylim(extrema(pd.y))
    plt.gcf().set_size_inches(18.0, 6.3, forward=true)
    plt.gca().set_aspect("equal")
    plt.savefig(joinpath(@__DIR__, "double_mach_reflection_$(name)_mesh.pdf"), bbox_inches="tight")

    plt.clf()
    data = first.(pd.data)
    pcm = plt.tripcolor(vec(pd.x), vec(pd.y), vec(data), cmap="inferno",
                         linewidth=0, rasterized=true)
    pcm.set_edgecolor("face")
    plt.colorbar()
    plt.xlabel(L"x")
    plt.ylabel(L"y")
    plt.xlim(extrema(pd.x))
    plt.ylim(extrema(pd.y))
    plt.gcf().set_size_inches(18.0, 6.3, forward=true)
    plt.gca().set_aspect("equal")
    plt.savefig(joinpath(@__DIR__, "double_mach_reflection_$(name)_density.pdf"), bbox_inches="tight")
  end
end
