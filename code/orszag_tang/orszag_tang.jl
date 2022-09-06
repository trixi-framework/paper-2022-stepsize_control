# Compare CFL-based and error-based step size control for entropy-dissipative
# semidiscretizations with shock capturing of the ideal generalized Lagrange
# multiplier (GLM) magnetohydrodynamics (MHD) equations for the Orszag-Tang
# vortex.

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using DelimitedFiles

using OrdinaryDiffEq
using Trixi


include("../setup_pyplot.jl")
include("../monitor_stepsize.jl")

function monitor_glm_speed()
  times = Vector{Float64}()
  glm_speeds = Vector{Float64}()

  cb = DiscreteCallback(
    (u, t, integrator) -> true, # condition
    (integrator) -> begin  # affect!
      t = integrator.t
      u_ode = integrator.u
      semi = integrator.p
      _, equations, _, _ = Trixi.mesh_equations_solver_cache(semi)
      c_h = equations.c_h
      push!(times, t)
      push!(glm_speeds, c_h)
      u_modified!(integrator, false)
    end,
    save_positions=(false,false),
    initialize=(cb, u, t, integrator) -> cb.affect!(integrator)
  )

  return times, glm_speeds, cb
end


"""
    initial_condition_orszag_tang(x, t, equations::IdealGlmMhdEquations2D)

The classical Orszag-Tang vortex test case. Here, the setup is taken from
- Dominik Derigs, Gregor J. Gassner, Stefanie Walch & Andrew R. Winters (2018)
  Entropy Stable Finite Volume Approximations for Ideal Magnetohydrodynamics
  [doi: 10.1365/s13291-018-0178-9](https://doi.org/10.1365/s13291-018-0178-9)
"""
function initial_condition_orszag_tang(x, t, equations::IdealGlmMhdEquations2D)
  # setup taken from Derigs et al. DMV article (2018)
  # domain must be [0, 1] x [0, 1], γ = 5/3
  rho = 1.0
  v1 = -sin(2.0*pi*x[2])
  v2 =  sin(2.0*pi*x[1])
  v3 = 0.0
  p = 1.0 / equations.gamma
  B1 = -sin(2.0*pi*x[2]) / equations.gamma
  B2 =  sin(4.0*pi*x[1]) / equations.gamma
  B3 = 0.0
  psi = 0.0
  return prim2cons(SVector(rho, v1, v2, v3, p, B1, B2, B3, psi), equations)
end


# Manually tuned CFL factors
orszag_tang_cfl(::BS3) = 0.69
orszag_tang_cfl(::RDPK3SpFSAL35) = 1.38
orszag_tang_cfl(::SSPRK43) = 1.40


function run_orszag_tang(alg; cfl=nothing, kwargs...)
  glm_scale = 1.0
  equations = IdealGlmMhdEquations2D(5/3, glm_scale)

  surface_flux = (flux_lax_friedrichs, flux_nonconservative_powell)
  volume_flux  = (flux_hindenlang_gassner, flux_nonconservative_powell)
  basis = LobattoLegendreBasis(3)
  indicator_sc = IndicatorHennemannGassner(equations, basis,
                                           alpha_max=0.5,
                                           alpha_min=0.001,
                                           alpha_smooth=true,
                                           variable=density_pressure)
  volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                   volume_flux_dg=volume_flux,
                                                   volume_flux_fv=surface_flux)
  solver = DGSEM(basis, surface_flux, volume_integral)

  coordinates_min = (0.0, 0.0)
  coordinates_max = (1.0, 1.0)
  mesh = TreeMesh(coordinates_min, coordinates_max,
                  initial_refinement_level=6,
                  n_cells_max=10_000)

  semi = SemidiscretizationHyperbolic(mesh, equations,
    initial_condition_orszag_tang, solver)

  tspan = (0.0, 0.5)
  ode = semidiscretize(semi, tspan)

  summary_callback = SummaryCallback()
  analysis_interval = 100
  analysis_callback = AnalysisCallback(semi, interval=analysis_interval)
  alive_callback = AliveCallback(; analysis_interval)
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



# Main function for this setup
function run_experiments()
  let alg = BS3()
    name = nameof(typeof(alg))
    # Controller suggested by Ranocha, Dalcin, Parsani, Ketcheson (2021)
    controller = PIDController(0.6, -0.2)

    sol_err, times_err, effective_cfl_numbers = run_orszag_tang(
      alg; adaptive=true, abstol=1.0e-4, reltol=1.0e-4, controller)

    @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
    open(joinpath(@__DIR__, "orszag_tang_$(name)_err.txt"), "w") do io
      println(io, "# nf = ",      sol_err.destats.nf)
      println(io, "# naccept = ", sol_err.destats.naccept)
      println(io, "# nreject = ", sol_err.destats.nreject)
      println(io, "# t cfl")
      writedlm(io, hcat(times_err, effective_cfl_numbers))
    end

    sol_cfl, times_cfl, effective_cfl_numbers_cfl = run_orszag_tang(
      alg; adaptive=false, cfl=orszag_tang_cfl(alg), dt=1.0 #= overwritten by CFL =#)

    @info "Finished simulation with CFL-based step size control" alg sol_cfl.destats.nf sol_cfl.destats.naccept sol_cfl.destats.nreject
    open(joinpath(@__DIR__, "orszag_tang_$(name)_cfl.txt"), "w") do io
      println(io, "# nf = ",      sol_cfl.destats.nf)
      println(io, "# naccept = ", sol_cfl.destats.naccept)
      println(io, "# nreject = ", sol_cfl.destats.nreject)
      println(io, "# t cfl")
      writedlm(io, hcat(times_cfl, effective_cfl_numbers_cfl))
    end

    # times_err, effective_cfl_numbers = let data = readdlm("orszag_tang_$(name)_err.txt", comments=true); data[:, 1], data[:, 2] end
    # times_cfl, effective_cfl_numbers_cfl = let data = readdlm("orszag_tang_$(name)_cfl.txt", comments=true); data[:, 1], data[:, 2] end
    # map(/, times_err[2:end] - times_err[1:end-1], times_cfl[2:end] - times_cfl[1:end-1])

    plt.clf()
    plt.plot(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Tol.")
    plt.plot(times_cfl[2:end], times_cfl[2:end]-times_cfl[1:end-1], label="CFL")
    plt.legend(loc="best")
    plt.xlabel("Time")
    plt.ylabel("Time step size")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_dt.png"), bbox_inches="tight")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_dt.pdf"), bbox_inches="tight")

    plt.clf()
    plt.plot(times_err, effective_cfl_numbers, label="Tol.")
    plt.plot(times_cfl, effective_cfl_numbers_cfl, label="CFL")
    plt.legend(loc="best")
    plt.xlabel("Time")
    plt.ylabel("Effective CFL number")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_cfl.png"), bbox_inches="tight")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_cfl.pdf"), bbox_inches="tight")
  end


  let alg = RDPK3SpFSAL35()
    name = nameof(typeof(alg))
    # Controller suggested by Ranocha, Dalcin, Parsani, Ketcheson (2021) is
    # already set as default in ordinaryDiffEq.jl

    sol_err, times_err, effective_cfl_numbers = run_orszag_tang(
      alg; adaptive=true, abstol=1.0e-4, reltol=1.0e-4)

    @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
    open(joinpath(@__DIR__, "orszag_tang_$(name)_err.txt"), "w") do io
      println(io, "# nf = ",      sol_err.destats.nf)
      println(io, "# naccept = ", sol_err.destats.naccept)
      println(io, "# nreject = ", sol_err.destats.nreject)
      println(io, "# t cfl")
      writedlm(io, hcat(times_err, effective_cfl_numbers))
    end

    sol_cfl, times_cfl, effective_cfl_numbers_cfl = run_orszag_tang(
      alg; adaptive=false, cfl=orszag_tang_cfl(alg), dt=1.0 #= overwritten by CFL =#)

    @info "Finished simulation with CFL-based step size control" alg sol_cfl.destats.nf sol_cfl.destats.naccept sol_cfl.destats.nreject
    open(joinpath(@__DIR__, "orszag_tang_$(name)_cfl.txt"), "w") do io
      println(io, "# nf = ",      sol_cfl.destats.nf)
      println(io, "# naccept = ", sol_cfl.destats.naccept)
      println(io, "# nreject = ", sol_cfl.destats.nreject)
      println(io, "# t cfl")
      writedlm(io, hcat(times_cfl, effective_cfl_numbers_cfl))
    end

    # times_err, effective_cfl_numbers = let data = readdlm("orszag_tang_$(name)_err.txt", comments=true); data[:, 1], data[:, 2] end
    # times_cfl, effective_cfl_numbers_cfl = let data = readdlm("orszag_tang_$(name)_cfl.txt", comments=true); data[:, 1], data[:, 2] end
    # map(/, times_err[2:end] - times_err[1:end-1], times_cfl[2:end] - times_cfl[1:end-1])

    plt.clf()
    plt.plot(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Tol.")
    plt.plot(times_cfl[2:end], times_cfl[2:end]-times_cfl[1:end-1], label="CFL")
    plt.legend(loc="best")
    plt.xlabel("Time")
    plt.ylabel("Time step size")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_dt.png"), bbox_inches="tight")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_dt.pdf"), bbox_inches="tight")

    plt.clf()
    plt.plot(times_err, effective_cfl_numbers, label="Tol.")
    plt.plot(times_cfl, effective_cfl_numbers_cfl, label="CFL")
    plt.legend(loc="best")
    plt.xlabel("Time")
    plt.ylabel("Effective CFL number")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_cfl.png"), bbox_inches="tight")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_cfl.pdf"), bbox_inches="tight")
  end


  let alg = SSPRK43()
    name = nameof(typeof(alg))
    # Controller suggested by Ranocha, Dalcin, Parsani, Ketcheson (2021)
    controller = PIDController(0.55, -0.27, 0.05)

    sol_err, times_err, effective_cfl_numbers = run_orszag_tang(
      alg; adaptive=true, abstol=1.0e-4, reltol=1.0e-4, controller)

    @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
    open(joinpath(@__DIR__, "orszag_tang_$(name)_err.txt"), "w") do io
      println(io, "# nf = ",      sol_err.destats.nf)
      println(io, "# naccept = ", sol_err.destats.naccept)
      println(io, "# nreject = ", sol_err.destats.nreject)
      println(io, "# t cfl")
      writedlm(io, hcat(times_err, effective_cfl_numbers))
    end

    sol_cfl, times_cfl, effective_cfl_numbers_cfl = run_orszag_tang(
      alg; adaptive=false, cfl=orszag_tang_cfl(alg), dt=1.0 #= overwritten by CFL =#)

    @info "Finished simulation with CFL-based step size control" alg sol_cfl.destats.nf sol_cfl.destats.naccept sol_cfl.destats.nreject
    open(joinpath(@__DIR__, "orszag_tang_$(name)_cfl.txt"), "w") do io
      println(io, "# nf = ",      sol_cfl.destats.nf)
      println(io, "# naccept = ", sol_cfl.destats.naccept)
      println(io, "# nreject = ", sol_cfl.destats.nreject)
      println(io, "# t cfl")
      writedlm(io, hcat(times_cfl, effective_cfl_numbers_cfl))
    end

    # times_err, effective_cfl_numbers = let data = readdlm("orszag_tang_$(name)_err.txt", comments=true); data[:, 1], data[:, 2] end
    # times_cfl, effective_cfl_numbers_cfl = let data = readdlm("orszag_tang_$(name)_cfl.txt", comments=true); data[:, 1], data[:, 2] end
    # map(/, times_err[2:end] - times_err[1:end-1], times_cfl[2:end] - times_cfl[1:end-1])

    plt.clf()
    plt.plot(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Tol.")
    plt.plot(times_cfl[2:end], times_cfl[2:end]-times_cfl[1:end-1], label="CFL")
    plt.legend(loc="best")
    plt.xlabel("Time")
    plt.ylabel("Time step size")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_dt.png"), bbox_inches="tight")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_dt.pdf"), bbox_inches="tight")

    plt.clf()
    plt.plot(times_err, effective_cfl_numbers, label="Tol.")
    plt.plot(times_cfl, effective_cfl_numbers_cfl, label="CFL")
    plt.legend(loc="best")
    plt.xlabel("Time")
    plt.ylabel("Effective CFL number")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_cfl.png"), bbox_inches="tight")
    plt.savefig(joinpath(@__DIR__, "orszag_tang_$(name)_cfl.pdf"), bbox_inches="tight")
  end
end

