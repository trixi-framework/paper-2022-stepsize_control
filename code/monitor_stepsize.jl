
# A callback monitoring the effective CFL factor of the current
# time step size after a successful time step.
# Assumes that Trixi and OrdinaryDiffEq are in scope.
# If the next predicted time step size should be motitored, you
# need to change
# - `u_ode = integrator.uprev` to `u_ode = integrator.u`
# - `integrator.dt` to `get_proposed_dt(integrator)`
function monitor_stepsize()
  times = Vector{Float64}()
  effective_cfl_numbers = Vector{Float64}()

  cb = DiscreteCallback(
    (u, t, integrator) -> true, # condition
    (integrator) -> begin  # affect!
      t = integrator.t
      u_ode = integrator.uprev
      semi = integrator.p
      mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
      u = Trixi.wrap_array(u_ode, mesh, equations, solver, cache)

      # Avoids crashes due to positivity issues when computing the CFL factor,
      # in particular for Gauss methods
      local dt::Float64
      try
        dt = Trixi.max_dt(u, t, mesh, Trixi.have_constant_speed(equations), equations,
                          solver, cache)
      catch e
        dt = NaN
        println(e)
      end

      push!(times, t)
      push!(effective_cfl_numbers, (integrator.t - integrator.tprev) / dt)
      u_modified!(integrator, false)
    end,
    save_positions=(false,false),
    # initialize=(cb, u, t, integrator) -> cb.affect!(integrator)
  )

  return times, effective_cfl_numbers, cb
end


# Throwing exceptions disables threads, see
# https://github.com/JuliaSIMD/Polyester.jl/issues/30
import PolyesterWeave, ThreadingUtilities

function reset_threads!()
  PolyesterWeave.reset_workers!()
  for i in 1:Threads.nthreads()-1
    ThreadingUtilities.initialize_task(i)
  end
end
