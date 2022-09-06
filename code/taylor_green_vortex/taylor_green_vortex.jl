# Compare CFL-based and error-based step size control for entropy-dissipative
# semidiscretizations of the compressible Euler equations of an ideal gas for
# the classical inviscid Taylor-Green vortex setup.

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using DelimitedFiles

using OrdinaryDiffEq
using Trixi
using Trixi2Vtk


include("../setup_pyplot.jl")
include("../monitor_stepsize.jl")


"""
    initial_condition_taylor_green_vortex(x, t, equations::CompressibleEulerEquations3D)

The classical inviscid Taylor-Green vortex.
"""
function initial_condition_taylor_green_vortex(x, t, equations::CompressibleEulerEquations3D)
  A  = 1.0 # magnitude of speed
  Ms = 0.1 # maximum Mach number

  rho = 1.0
  v1  =  A * sin(x[1]) * cos(x[2]) * cos(x[3])
  v2  = -A * cos(x[1]) * sin(x[2]) * cos(x[3])
  v3  = 0.0
  p   = (A / Ms)^2 * rho / equations.gamma # scaling to get Ms
  p   = p + 1.0/16.0 * A^2 * rho * (cos(2*x[1])*cos(2*x[3]) + 2*cos(2*x[2]) + 2*cos(2*x[1]) + cos(2*x[2])*cos(2*x[3]))

  return prim2cons(SVector(rho, v1, v2, v3, p), equations)
end


# Mapping based on
# - Section 5.2 of
#   Chan, Fernández, Carpenter (2019)
#   Efficient entropy stable Gauss collocation methods
# - Section 6.1 of
#   https://arξv.org/abs/2012.12040
function mapping(ξ_, η_, ζ_)
  # transform input variables from [-1, 1] onto [-π, π]
  ξ = π * ξ_
  η = π * η_
  ζ = π * ζ_

  # box lengths in physical coordinates
  Lx = Ly = Lz = 2 * π

  # boc centers in physical coordinates
  cx = cy = cz = 0.0

  y = η + Ly / 8 * (cos(3π * (ξ - cx) / Lx) *
                    cos( π * (η - cy) / Ly) *
                    cos( π * (ζ - cz) / Lz))

  x = ξ + Lx / 8 * (cos( π * (ξ - cx) / Lx) *
                    cos(4π * (y - cy) / Ly) *
                    cos( π * (ζ - cz) / Lz))

  z = ζ + Lz / 8 * (cos( π * (x - cx) / Lx) *
                    cos(2π * (y - cy) / Ly) *
                    cos( π * (ζ - cz) / Lz))

  return SVector(x, y, z)
end


# Manually tuned CFL factors
taylor_green_vortex_cfl(::BS3, ::Type{TreeMesh}) = 0.87
taylor_green_vortex_cfl(::RDPK3SpFSAL35, ::Type{TreeMesh}) = 1.71
taylor_green_vortex_cfl(::SSPRK43, ::Type{TreeMesh}) = 1.59

taylor_green_vortex_cfl(::BS3, ::Type{StructuredMesh}) = 1.35
taylor_green_vortex_cfl(::RDPK3SpFSAL35, ::Type{StructuredMesh}) = 2.67
taylor_green_vortex_cfl(::SSPRK43, ::Type{StructuredMesh}) = 2.77


function run_taylor_green_vortex(alg, mesh_type; cfl=nothing, kwargs...)
  equations = CompressibleEulerEquations3D(1.4)

  volume_flux = flux_ranocha_turbo
  solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs,
                 volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

  coordinates_min = (-1.0, -1.0, -1.0) .* pi
  coordinates_max = ( 1.0,  1.0,  1.0) .* pi
  initial_refinement_level = 3
  cells_per_dimension = 2^initial_refinement_level .* (1, 1, 1)
  if mesh_type == TreeMesh
    mesh = TreeMesh(coordinates_min, coordinates_max;
                    initial_refinement_level,
                    n_cells_max=10_000, periodicity=true)
  elseif mesh_type == StructuredMesh
    mesh = StructuredMesh(cells_per_dimension, mapping, periodicity=true)
  end

  semi = SemidiscretizationHyperbolic(mesh, equations,
    initial_condition_taylor_green_vortex, solver)

  tspan = (0.0, 10.0)
  ode = semidiscretize(semi, tspan)

  summary_callback = SummaryCallback()
  analysis_interval = 5000
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
  generate_data()
  generate_plots()
  generate_visualization_files()
end


function generate_data()
  # Algorithms and controller suggested by Ranocha, Dalcin, Parsani, Ketcheson (2021)
  parameters = [
    (BS3(), PIDController(0.60, -0.20, 0.00)),
    (RDPK3SpFSAL35(), PIDController(0.70, -0.23, 0.00)),
    (SSPRK43(), PIDController(0.55, -0.27, 0.05)),
  ]

  tolerances = round.(10.0 .^ range(-10, -2, 25), sigdigits=3)

  for mesh_type in [#=TreeMesh,=# StructuredMesh]
    for (alg, controller) in parameters
      name = nameof(typeof(alg))
      @info "Running" alg mesh_type

      rhs_evaluations = Int[]
      for tol in tolerances
        try
          sol_err, _ = run_taylor_green_vortex(alg, mesh_type,
              adaptive=true, controller=controller, abstol=tol, reltol=tol)
          if sol_err.retcode !== :Success
            error("retcode = $(sol_err.retcode)")
          end
          push!(rhs_evaluations, sol_err.destats.nf)
        catch e
          @warn "Crashed" e alg mesh_type tol
          reset_threads!()
          push!(rhs_evaluations, typemin(Int))
        end
      end

      cfl = taylor_green_vortex_cfl(alg, mesh_type)
      sol_cfl, _ = run_taylor_green_vortex(alg, mesh_type;
          adaptive=false, cfl=cfl, dt=1.0 #= overwritten by CFL =#)

      open(joinpath(@__DIR__, "taylor_green_vortex_$(name)_$(mesh_type).txt"), "w") do io
        println(io, "# CFL: nf = ", sol_cfl.destats.nf)
        println(io, "# CFL: cfl = ", cfl)
        println(io, "# tol nf")
        writedlm(io, hcat(tolerances, rhs_evaluations))
      end

      @info "Finished running" alg mesh_type
    end
  end
end


function generate_plots()
  parameters = [
    (BS3(), "\$\\mathrm{BS3(2)3_F}\$", "#E69F00", "-", "4"),
    (RDPK3SpFSAL35(), "\$\\mathrm{RDPK3(2)5_F[3S^*_+]}\$", "#56B4E9", "--", "2"),
    (SSPRK43(), "\$\\mathrm{SSP3(2)4[3S^*_+]}\$", "#009E73", "-.", "3"),
  ]

  for mesh_type in [TreeMesh, StructuredMesh]
    plt.clf()
    for (alg, label, color, linestyle, marker) in parameters
      name = nameof(typeof(alg))
      tolerances, rhs_evaluations =  let data = readdlm("taylor_green_vortex_$(name)_$(mesh_type).txt", comments=true); data[:, 1], data[:, 2] end

      data = readdlm("taylor_green_vortex_$(name)_$(mesh_type).txt")
      rhs_evaluations_cfl = data[1, end]
      cfl = data[2, end]

      plt.loglog(tolerances, rhs_evaluations; label, color, marker, linestyle="none")
      plt.plot(extrema(tolerances), rhs_evaluations_cfl * ones(2); color, linestyle,
               label="with opt. \$\\nu = $(cfl)\$")
      @info " Plotting" mesh_type alg cfl
    end

    plt.legend(loc="best")
    plt.xlabel("Tolerance")
    plt.ylabel("#RHS evaluations")
    plt.ylim((4.0e3, 9.0e4)) # plt.ylim((4.0e3, 4.5e4))

    plt.savefig(joinpath(@__DIR__, "taylor_green_vortex_$(mesh_type).png"), bbox_inches="tight")
    plt.savefig(joinpath(@__DIR__, "taylor_green_vortex_$(mesh_type).pdf"), bbox_inches="tight")
  end
end


function generate_visualization_files()
  for mesh_type in [TreeMesh, StructuredMesh]
    isdir("out") && rm("out", recursive=true)

    equations = CompressibleEulerEquations3D(1.4)

    volume_flux = flux_ranocha_turbo
    solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs,
                  volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

    coordinates_min = (-1.0, -1.0, -1.0) .* pi
    coordinates_max = ( 1.0,  1.0,  1.0) .* pi
    initial_refinement_level = 3
    cells_per_dimension = 2^initial_refinement_level .* (1, 1, 1)
    if mesh_type == TreeMesh
      mesh = TreeMesh(coordinates_min, coordinates_max;
                      initial_refinement_level,
                      n_cells_max=10_000, periodicity=true)
    elseif mesh_type == StructuredMesh
      mesh = StructuredMesh(cells_per_dimension, mapping, periodicity=true)
    end

    semi = SemidiscretizationHyperbolic(mesh, equations,
      initial_condition_taylor_green_vortex, solver)

    tspan = (0.0, 0.0)
    ode = semidiscretize(semi, tspan)

    save_solution = SaveSolutionCallback(interval=100,
                                         save_initial_solution=true,
                                         save_final_solution=false)

    solve(ode, BS3(), save_everystep=false, callback=save_solution)

    trixi2vtk("out/solution_*.h5",
              output_directory=joinpath(@__DIR__, "out_$(mesh_type)"),
              nvisnodes=19)
    isdir("out") && rm("out", recursive=true)
  end
end
