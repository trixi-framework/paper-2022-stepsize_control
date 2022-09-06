# Compare CFL-based and error-based step size control for entropy-dissipative
# semidiscretizations with entropy projections on Gauss nodes of the
# compressible Euler equations for a Kelvin-Helmholtz instability.

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using DelimitedFiles

using OrdinaryDiffEq
using Trixi


include("../setup_pyplot.jl")
include("../monitor_stepsize.jl")


"""
    initial_condition_kelvin_helmholtz_instability(x, t, equations::CompressibleEulerEquations2D)

A version of the classical Kelvin-Helmholtz instability based on
- Andrés M. Rueda-Ramírez, Gregor J. Gassner (2021)
  A Subcell Finite Volume Positivity-Preserving Limiter for DGSEM Discretizations
  of the Euler Equations
  [arXiv: 2102.06017](https://arxiv.org/abs/2102.06017)

with low Atwood number `A = 3/7` (density ratio).
"""
function initial_condition_kelvin_helmholtz_instability_low_atwood(x, t, equations::CompressibleEulerEquations2D)
  A = 3/7 # Atwood number (density ratio)
  rho1 = 0.5 * one(A) # recover original with A = 3/7
  rho2 = rho1 * (1 + A) / (1 - A)

  # B is a discontinuous function with value 1 for -.5 <= x <= .5 and 0 elsewhere
  slope = 15
  B = 0.5 * (tanh(slope * x[2] + 7.5) - tanh(slope * x[2] - 7.5))

  rho = rho1 + B * (rho2 - rho1)  # rho ∈ [rho_1, rho_2]
  v1 = B - 0.5                    # v1  ∈ [-.5, .5]
  v2 = 0.1 * sin(2 * pi * x[1])
  p = 1.0
  return prim2cons(SVector(rho, v1, v2, p), equations)
end

"""
    initial_condition_kelvin_helmholtz_instability(x, t, equations::CompressibleEulerEquations2D)

A version of the classical Kelvin-Helmholtz instability based on
- Andrés M. Rueda-Ramírez, Gregor J. Gassner (2021)
  A Subcell Finite Volume Positivity-Preserving Limiter for DGSEM Discretizations
  of the Euler Equations
  [arXiv: 2102.06017](https://arxiv.org/abs/2102.06017)

with high Atwood number `A = 0.7` (density ratio).
"""
function initial_condition_kelvin_helmholtz_instability_high_atwood(x, t, equations::CompressibleEulerEquations2D)
  A = 0.7 # Atwood number (density ratio)
  rho1 = 0.5 * one(A) # recover original with A = 3/7
  rho2 = rho1 * (1 + A) / (1 - A)

  # B is a discontinuous function with value 1 for -.5 <= x <= .5 and 0 elsewhere
  slope = 15
  B = 0.5 * (tanh(slope * x[2] + 7.5) - tanh(slope * x[2] - 7.5))

  rho = rho1 + B * (rho2 - rho1)  # rho ∈ [rho_1, rho_2]
  v1 = B - 0.5                    # v1  ∈ [-.5, .5]
  v2 = 0.1 * sin(2 * pi * x[1])
  p = 1.0
  return prim2cons(SVector(rho, v1, v2, p), equations)
end

shortname(x) = string(x)
shortname(::typeof(initial_condition_kelvin_helmholtz_instability_low_atwood)) = "lowAtwood"
shortname(::typeof(initial_condition_kelvin_helmholtz_instability_high_atwood)) = "highAtwood"


# Manually tuned CFL factors
kelvin_helmholtz_cfl(::BS3,
  ::typeof(initial_condition_kelvin_helmholtz_instability_low_atwood)) = 0.57
# kelvin_helmholtz_cfl(::BS3,
#   ::typeof(initial_condition_kelvin_helmholtz_instability_high_atwood)) =
kelvin_helmholtz_cfl(::RDPK3SpFSAL35,
  ::typeof(initial_condition_kelvin_helmholtz_instability_low_atwood)) = 1.10
# kelvin_helmholtz_cfl(::RDPK3SpFSAL35,
#   ::typeof(initial_condition_kelvin_helmholtz_instability_high_atwood)) =
kelvin_helmholtz_cfl(::RDPK3SpFSAL49,
::typeof(initial_condition_kelvin_helmholtz_instability_low_atwood)) = 1.85
# kelvin_helmholtz_cfl(::RDPK3SpFSAL49,
#   ::typeof(initial_condition_kelvin_helmholtz_instability_high_atwood)) =
kelvin_helmholtz_cfl(::SSPRK43,
  ::typeof(initial_condition_kelvin_helmholtz_instability_low_atwood)) = 1.03
# kelvin_helmholtz_cfl(::SSPRK43,
#   ::typeof(initial_condition_kelvin_helmholtz_instability_high_atwood)) =


function run_kelvin_helmholtz(alg, initial_condition; cfl=nothing, kwargs...)
  dg = DGMulti(polydeg = 3, element_type = Quad(), approximation_type = GaussSBP(),
               surface_integral = SurfaceIntegralWeakForm(flux_lax_friedrichs),
               volume_integral = VolumeIntegralFluxDifferencing(flux_ranocha))

  equations = CompressibleEulerEquations2D(1.4)

  mesh = DGMultiMesh(dg, cells_per_dimension=(32, 32), periodicity=true)
  semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg)

  tspan = (0.0, 5.0)
  ode = semidiscretize(semi, tspan)

  summary_callback = SummaryCallback()
  alive_callback = AliveCallback(alive_interval=50)
  basic_callbacks = (summary_callback, alive_callback)

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

    let initial_condition = initial_condition_kelvin_helmholtz_instability_low_atwood
      name = string(nameof(typeof(alg))) * "_" *  shortname(initial_condition)

      # tolerances 1.0e-4, 1.0e-5, 1.0e-6 yield basically the same results
      sol_err, times_err, effective_cfl_numbers = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=true, controller,
        abstol=1.0e-4, reltol=1.0e-4)

      @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_err.txt"), "w") do io
        println(io, "# nf = ",      sol_err.destats.nf)
        println(io, "# naccept = ", sol_err.destats.naccept)
        println(io, "# nreject = ", sol_err.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_err, effective_cfl_numbers))
      end

      sol_cfl, times_cfl, effective_cfl_numbers_cfl = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=false, dt=1.0, #= overwritten by CFL =#
        cfl=kelvin_helmholtz_cfl(alg, initial_condition))

      @info "Finished simulation with CFL-based step size control" alg sol_cfl.destats.nf sol_cfl.destats.naccept sol_cfl.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.txt"), "w") do io
        println(io, "# nf = ",      sol_cfl.destats.nf)
        println(io, "# naccept = ", sol_cfl.destats.naccept)
        println(io, "# nreject = ", sol_cfl.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_cfl, effective_cfl_numbers_cfl))
      end

      # times_err, effective_cfl_numbers = let data = readdlm("kelvin_helmholtz_$(name)_err.txt", comments=true); data[:, 1], data[:, 2] end
      # times_cfl, effective_cfl_numbers_cfl = let data = readdlm("kelvin_helmholtz_$(name)_cfl.txt", comments=true); data[:, 1], data[:, 2] end
      # map(/, times_err[2:end] - times_err[1:end-1], times_cfl[2:end] - times_cfl[1:end-1])

      plt.clf()
      plt.plot(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Tol.")
      plt.plot(times_cfl[2:end], times_cfl[2:end]-times_cfl[1:end-1], label="CFL")
      plt.legend(loc="best")
      plt.xlabel("Time")
      plt.ylabel("Time step size")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_dt.png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_dt.pdf"), bbox_inches="tight")

      plt.clf()
      plt.plot(times_err, effective_cfl_numbers, label="Tol.")
      plt.plot(times_cfl, effective_cfl_numbers_cfl, label="CFL")
      plt.legend(loc="best")
      plt.xlabel("Time")
      plt.ylabel("Effective CFL number")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.pdf"), bbox_inches="tight")
    end

    let initial_condition = initial_condition_kelvin_helmholtz_instability_high_atwood
      name = string(nameof(typeof(alg))) * "_" *  shortname(initial_condition)

      # tolerances 1.0e-4, 1.0e-5, 1.0e-6 yield basically the same results
      sol_err, times_err, effective_cfl_numbers = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=true, controller,
        abstol=1.0e-4, reltol=1.0e-4)

      @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_err.txt"), "w") do io
        println(io, "# nf = ",      sol_err.destats.nf)
        println(io, "# naccept = ", sol_err.destats.naccept)
        println(io, "# nreject = ", sol_err.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_err, effective_cfl_numbers))
      end

      plt.clf()
      ax1 = plt.gca()
      plt.plot(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Time step size")
      plt.legend(loc="upper left")
      plt.xlabel("Time")
      plt.ylabel("Time step size")
      ax2 = ax1.twinx()
      plt.plot([], [])
      plt.plot(times_err, effective_cfl_numbers, label="Effective CFL number")
      plt.legend(loc="upper right", bbox_to_anchor=(1.0, 0.9))
      plt.ylabel("Effective CFL number")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name).png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name).pdf"), bbox_inches="tight")
    end
  end


  let alg = RDPK3SpFSAL35()
    # Controller suggested by Ranocha, Dalcin, Parsani, Ketcheson (2021) is
    # already set as default in ordinaryDiffEq.jl

    let initial_condition = initial_condition_kelvin_helmholtz_instability_low_atwood
      name = string(nameof(typeof(alg))) * "_" *  shortname(initial_condition)

      # tolerances 1.0e-4, 1.0e-5 yield basically the same results
      sol_err, times_err, effective_cfl_numbers = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=true,
        abstol=1.0e-4, reltol=1.0e-4)

      @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_err.txt"), "w") do io
        println(io, "# nf = ",      sol_err.destats.nf)
        println(io, "# naccept = ", sol_err.destats.naccept)
        println(io, "# nreject = ", sol_err.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_err, effective_cfl_numbers))
      end

      sol_cfl, times_cfl, effective_cfl_numbers_cfl = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=false, dt=1.0, #= overwritten by CFL =#
        cfl=kelvin_helmholtz_cfl(alg, initial_condition))

      @info "Finished simulation with CFL-based step size control" alg sol_cfl.destats.nf sol_cfl.destats.naccept sol_cfl.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.txt"), "w") do io
        println(io, "# nf = ",      sol_cfl.destats.nf)
        println(io, "# naccept = ", sol_cfl.destats.naccept)
        println(io, "# nreject = ", sol_cfl.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_cfl, effective_cfl_numbers_cfl))
      end

      # times_err, effective_cfl_numbers = let data = readdlm("kelvin_helmholtz_$(name)_err.txt", comments=true); data[:, 1], data[:, 2] end
      # times_cfl, effective_cfl_numbers_cfl = let data = readdlm("kelvin_helmholtz_$(name)_cfl.txt", comments=true); data[:, 1], data[:, 2] end
      # map(/, times_err[2:end] - times_err[1:end-1], times_cfl[2:end] - times_cfl[1:end-1])

      plt.clf()
      plt.plot(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Tol.")
      plt.plot(times_cfl[2:end], times_cfl[2:end]-times_cfl[1:end-1], label="CFL")
      plt.legend(loc="best")
      plt.xlabel("Time")
      plt.ylabel("Time step size")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_dt.png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_dt.pdf"), bbox_inches="tight")

      plt.clf()
      plt.plot(times_err, effective_cfl_numbers, label="Tol.")
      plt.plot(times_cfl, effective_cfl_numbers_cfl, label="CFL")
      plt.legend(loc="best")
      plt.xlabel("Time")
      plt.ylabel("Effective CFL number")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.pdf"), bbox_inches="tight")
    end

    let initial_condition = initial_condition_kelvin_helmholtz_instability_high_atwood
      name = string(nameof(typeof(alg))) * "_" *  shortname(initial_condition)

      # tolerances 1.0e-4, 1.0e-5 yield basically the same results
      sol_err, times_err, effective_cfl_numbers = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=true,
        abstol=1.0e-4, reltol=1.0e-4)

      @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_err.txt"), "w") do io
        println(io, "# nf = ",      sol_err.destats.nf)
        println(io, "# naccept = ", sol_err.destats.naccept)
        println(io, "# nreject = ", sol_err.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_err, effective_cfl_numbers))
      end

      plt.clf()
      ax1 = plt.gca()
      plt.plot(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Time step size")
      plt.legend(loc="upper left")
      plt.xlabel("Time")
      plt.ylabel("Time step size")
      ax2 = ax1.twinx()
      plt.plot([], [])
      plt.plot(times_err, effective_cfl_numbers, label="Effective CFL number")
      plt.legend(loc="upper right", bbox_to_anchor=(1.0, 0.9))
      plt.ylabel("Effective CFL number")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name).png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name).pdf"), bbox_inches="tight")
    end
  end


  let alg = RDPK3SpFSAL49()
    # Controller suggested by Ranocha, Dalcin, Parsani, Ketcheson (2021) is
    # already set as default in ordinaryDiffEq.jl

    let initial_condition = initial_condition_kelvin_helmholtz_instability_low_atwood
      name = string(nameof(typeof(alg))) * "_" *  shortname(initial_condition)

      sol_err, times_err, effective_cfl_numbers = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=true,
        abstol=1.0e-5, reltol=1.0e-5)

      @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_err.txt"), "w") do io
        println(io, "# nf = ",      sol_err.destats.nf)
        println(io, "# naccept = ", sol_err.destats.naccept)
        println(io, "# nreject = ", sol_err.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_err, effective_cfl_numbers))
      end

      sol_cfl, times_cfl, effective_cfl_numbers_cfl = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=false, dt=1.0, #= overwritten by CFL =#
        cfl=kelvin_helmholtz_cfl(alg, initial_condition))

      @info "Finished simulation with CFL-based step size control" alg sol_cfl.destats.nf sol_cfl.destats.naccept sol_cfl.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.txt"), "w") do io
        println(io, "# nf = ",      sol_cfl.destats.nf)
        println(io, "# naccept = ", sol_cfl.destats.naccept)
        println(io, "# nreject = ", sol_cfl.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_cfl, effective_cfl_numbers_cfl))
      end

      # times_err, effective_cfl_numbers = let data = readdlm("kelvin_helmholtz_$(name)_err.txt", comments=true); data[:, 1], data[:, 2] end
      # times_cfl, effective_cfl_numbers_cfl = let data = readdlm("kelvin_helmholtz_$(name)_cfl.txt", comments=true); data[:, 1], data[:, 2] end
      # map(/, times_err[2:end] - times_err[1:end-1], times_cfl[2:end] - times_cfl[1:end-1])

      plt.clf()
      plt.plot(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Tol.")
      plt.plot(times_cfl[2:end], times_cfl[2:end]-times_cfl[1:end-1], label="CFL")
      plt.legend(loc="best")
      plt.xlabel("Time")
      plt.ylabel("Time step size")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_dt.png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_dt.pdf"), bbox_inches="tight")

      plt.clf()
      plt.plot(times_err, effective_cfl_numbers, label="Tol.")
      plt.plot(times_cfl, effective_cfl_numbers_cfl, label="CFL")
      plt.legend(loc="best")
      plt.xlabel("Time")
      plt.ylabel("Effective CFL number")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.pdf"), bbox_inches="tight")
    end

    let initial_condition = initial_condition_kelvin_helmholtz_instability_high_atwood
      name = string(nameof(typeof(alg))) * "_" *  shortname(initial_condition)

      sol_err, times_err, effective_cfl_numbers = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=true,
        abstol=1.0e-6, reltol=1.0e-6)

      @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_err.txt"), "w") do io
        println(io, "# nf = ",      sol_err.destats.nf)
        println(io, "# naccept = ", sol_err.destats.naccept)
        println(io, "# nreject = ", sol_err.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_err, effective_cfl_numbers))
      end

      plt.clf()
      ax1 = plt.gca()
      plt.plot(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Time step size")
      plt.legend(loc="upper left")
      plt.xlabel("Time")
      plt.ylabel("Time step size")
      ax2 = ax1.twinx()
      plt.plot([], [])
      plt.plot(times_err, effective_cfl_numbers, label="Effective CFL number")
      plt.legend(loc="upper right", bbox_to_anchor=(1.0, 0.9))
      plt.ylabel("Effective CFL number")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name).png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name).pdf"), bbox_inches="tight")
    end
  end


  let alg = SSPRK43()
    # Controller suggested by Ranocha, Dalcin, Parsani, Ketcheson (2021)
    controller = PIDController(0.55, -0.27, 0.05)

    let initial_condition = initial_condition_kelvin_helmholtz_instability_low_atwood
      name = string(nameof(typeof(alg))) * "_" *  shortname(initial_condition)

      sol_err, times_err, effective_cfl_numbers = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=true, controller,
        abstol=1.0e-4, reltol=1.0e-4)

      @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_err.txt"), "w") do io
        println(io, "# nf = ",      sol_err.destats.nf)
        println(io, "# naccept = ", sol_err.destats.naccept)
        println(io, "# nreject = ", sol_err.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_err, effective_cfl_numbers))
      end

      sol_cfl, times_cfl, effective_cfl_numbers_cfl = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=false, dt=1.0, #= overwritten by CFL =#
        cfl=kelvin_helmholtz_cfl(alg, initial_condition))

      @info "Finished simulation with CFL-based step size control" alg sol_cfl.destats.nf sol_cfl.destats.naccept sol_cfl.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.txt"), "w") do io
        println(io, "# nf = ",      sol_cfl.destats.nf)
        println(io, "# naccept = ", sol_cfl.destats.naccept)
        println(io, "# nreject = ", sol_cfl.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_cfl, effective_cfl_numbers_cfl))
      end

      # times_err, effective_cfl_numbers = let data = readdlm("kelvin_helmholtz_$(name)_err.txt", comments=true); data[:, 1], data[:, 2] end
      # times_cfl, effective_cfl_numbers_cfl = let data = readdlm("kelvin_helmholtz_$(name)_cfl.txt", comments=true); data[:, 1], data[:, 2] end
      # map(/, times_err[2:end] - times_err[1:end-1], times_cfl[2:end] - times_cfl[1:end-1])

      plt.clf()
      plt.plot(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Tol.")
      plt.plot(times_cfl[2:end], times_cfl[2:end]-times_cfl[1:end-1], label="CFL")
      plt.legend(loc="best")
      plt.xlabel("Time")
      plt.ylabel("Time step size")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_dt.png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_dt.pdf"), bbox_inches="tight")

      plt.clf()
      plt.plot(times_err, effective_cfl_numbers, label="Tol.")
      plt.plot(times_cfl, effective_cfl_numbers_cfl, label="CFL")
      plt.legend(loc="best")
      plt.xlabel("Time")
      plt.ylabel("Effective CFL number")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_cfl.pdf"), bbox_inches="tight")
    end

    let initial_condition = initial_condition_kelvin_helmholtz_instability_high_atwood
      name = string(nameof(typeof(alg))) * "_" *  shortname(initial_condition)

      sol_err, times_err, effective_cfl_numbers = run_kelvin_helmholtz(
        alg, initial_condition; adaptive=true, controller,
        abstol=1.0e-4, reltol=1.0e-4)

      @info "Finished simulation with error-based step size control" alg sol_err.destats.nf sol_err.destats.naccept sol_err.destats.nreject
      open(joinpath(@__DIR__, "kelvin_helmholtz_$(name)_err.txt"), "w") do io
        println(io, "# nf = ",      sol_err.destats.nf)
        println(io, "# naccept = ", sol_err.destats.naccept)
        println(io, "# nreject = ", sol_err.destats.nreject)
        println(io, "# t cfl")
        writedlm(io, hcat(times_err, effective_cfl_numbers))
      end

      plt.clf()
      ax1 = plt.gca()
      plt.plot(times_err[2:end], times_err[2:end]-times_err[1:end-1], label="Time step size")
      plt.legend(loc="upper left")
      plt.xlabel("Time")
      plt.ylabel("Time step size")
      ax2 = ax1.twinx()
      plt.plot([], [])
      plt.plot(times_err, effective_cfl_numbers, label="Effective CFL number")
      plt.legend(loc="upper right", bbox_to_anchor=(1.0, 0.9))
      plt.ylabel("Effective CFL number")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name).png"), bbox_inches="tight")
      plt.savefig(joinpath(@__DIR__, "kelvin_helmholtz_$(name).pdf"), bbox_inches="tight")
    end
  end
end

