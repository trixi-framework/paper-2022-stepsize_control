# Cold startup of an astrophysical jet simulation of the compressible Euler
# equations discretized in space with entropy-dissipative shock capturing DG
# methods, positivity-preserving limiters, and adaptive mesh refinement in
# Trixi.jl.

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using DelimitedFiles

using OrdinaryDiffEq
using Trixi


include("../setup_pyplot.jl")
include("../monitor_stepsize.jl")


"""
    initial_condition_astro_jet(x, t, equations::IdealGlmMhdEquations2D)

An initial condition of a Mach = 2000 jet adopted from
- Yong Liu, Jianfang Lu, and Chi-Wang Shu
  An oscillation free discontinuous Galerkin method for hyperbolic systems
  https://tinyurl.com/c76fjtx4
"""
function initial_condition_astro_jet(x, t, equations::CompressibleEulerEquations2D)
  rho = 0.5
  v1 = 0
  v2 = 0
  p =  0.4127
  # add inflow for t>0 at x=-0.5
  # domain size is [-0.5,+0.5]^2
  if (t > 0) && (x[1] ≈ -0.5) && (abs(x[2]) < 0.05)
    rho = 5
    v1 = 800 # about Mach number Ma = 2000
    v2 = 0
    p = 0.4127
  end
  return prim2cons(SVector(rho, v1, v2, p), equations)
end


function run_astro_jet(_alg; cfl=nothing, kwargs...)
  equations = CompressibleEulerEquations2D(5/3)

  boundary_conditions = (
    x_neg = BoundaryConditionDirichlet(initial_condition_astro_jet),
    x_pos = BoundaryConditionDirichlet(initial_condition_astro_jet),
    y_neg = boundary_condition_periodic,
    y_pos = boundary_condition_periodic,
  )

  surface_flux = flux_lax_friedrichs
  volume_flux  = flux_ranocha_turbo
  basis = LobattoLegendreBasis(3)
  indicator_sc = IndicatorHennemannGassner(equations, basis,
                                           alpha_max=0.3,
                                           alpha_min=0.001,
                                           alpha_smooth=true,
                                           variable=density_pressure)
  volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                   volume_flux_dg=volume_flux,
                                                   volume_flux_fv=surface_flux)
  solver = DGSEM(basis, surface_flux, volume_integral)

  coordinates_min = (-0.5, -0.5)
  coordinates_max = ( 0.5,  0.5)
  mesh = TreeMesh(coordinates_min, coordinates_max,
                  initial_refinement_level=6,
                  periodicity=(false,true),
                  n_cells_max=10^5)

  semi = SemidiscretizationHyperbolic(mesh, equations,
    initial_condition_astro_jet, solver; boundary_conditions)

  tspan = (0.0, 0.001)
  ode = semidiscretize(semi, tspan)

  summary_callback = SummaryCallback()
  analysis_interval = 5000
  analysis_callback = AnalysisCallback(semi, interval=analysis_interval)
  alive_callback = AliveCallback(; analysis_interval)

  times, effective_cfl_numbers, monitor_stepsize_callback = monitor_stepsize()

  amr_indicator = IndicatorHennemannGassner(semi,
                                            alpha_max=1.0,
                                            alpha_min=0.0001,
                                            alpha_smooth=false,
                                            variable=density)
  amr_controller = ControllerThreeLevelCombined(semi, amr_indicator, indicator_sc,
                                    base_level=2,
    #= med_level = current level =# med_level =0, med_threshold=0.0003,
                                    max_level =8, max_threshold=0.003,
                                    max_threshold_secondary=indicator_sc.alpha_max)
  amr_callback = AMRCallback(semi, amr_controller,
                             interval=1,
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
    thresholds=(5.0e-6, 5.0e-6),
    variables=(density, pressure))
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

    sol_err, times_err, effective_cfl_numbers = run_astro_jet(
      alg; adaptive=true, dt=1.0e-12, abstol=1.0e-4, reltol=1.0e-4, controller)

    @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
    open(joinpath(@__DIR__, "astro_jet_$(name)_err.txt"), "w") do io
      println(io, "# nf = ",      sol_err.destats.nf)
      println(io, "# naccept = ", sol_err.destats.naccept)
      println(io, "# nreject = ", sol_err.destats.nreject)
      println(io, "# t cfl")
      writedlm(io, hcat(times_err, effective_cfl_numbers))
    end

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
    plt.savefig(joinpath(@__DIR__, "astro_jet_$(name)_time.png"), bbox_inches="tight")
    plt.savefig(joinpath(@__DIR__, "astro_jet_$(name)_time.pdf"), bbox_inches="tight")

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
    plt.savefig(joinpath(@__DIR__, "astro_jet_$(name)_step.png"), bbox_inches="tight")
    plt.savefig(joinpath(@__DIR__, "astro_jet_$(name)_step.pdf"), bbox_inches="tight")

    # hacks using internal knowledge of the plotting data structures
    pd = PlotData2D(sol_err)

    plt.clf()
    plt.plot(pd.mesh_vertices_x, pd.mesh_vertices_y, color="black", linewidth=2)
    plt.xlabel(L"x")
    plt.ylabel(L"y")
    plt.xlim(extrema(pd.x))
    plt.ylim(extrema(pd.y))
    plt.savefig(joinpath(@__DIR__, "astro_jet_$(name)_mesh.pdf"), bbox_inches="tight")

    colors = PyCall.pyimport("matplotlib").colors
    data = pd.data[1] # density
    plt.clf()
    pcm = plt.pcolormesh(pd.x, pd.y, data, shading="auto", cmap="inferno",
                linewidth=0, rasterized=true,
                norm=colors.LogNorm(vmin=minimum(data), vmax=maximum(data)))
    pcm.set_edgecolor("face")
    plt.colorbar()
    plt.xlabel(L"x")
    plt.ylabel(L"y")
    plt.xlim(extrema(pd.x))
    plt.ylim(extrema(pd.y))
    plt.savefig(joinpath(@__DIR__, "astro_jet_$(name)_density.pdf"), bbox_inches="tight")
  end
end
