# Compute spectra of plain DG and shock capturing semidiscretizations and their
# corresponding linear CFL factors for a range of RK methods.

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using LinearAlgebra

using DiffEqDevTools
using OrdinaryDiffEq
using Trixi

include("../setup_pyplot.jl")


# An indicator mocking the interface of Trixi.jl for shock capturing indicators
# with a fixed value.
struct IndicatorFixed{Cache<:NamedTuple} <: Trixi.AbstractIndicator
  cache::Cache
end

function IndicatorFixed(equations, basis; alpha_fixed)
  alpha = Vector{real(basis)}()
  cache = (; alpha, alpha_fixed)

  return IndicatorFixed{typeof(cache)}(cache)
end

function (indicator::IndicatorFixed)(u, mesh, equations, dg, cache;
                                     kwargs...)
  alpha = indicator.cache.alpha
  resize!(alpha, nelements(dg, cache))
  alpha .= indicator.cache.alpha_fixed
  return alpha
end


# Compute the spectra of plain DGSEM and shock capturing semidiscretizations
# for linear advection in 2D.
# You can call it as
# λ_dg, λ_sc = compute_spectra(angle=π/4, shock_capturing_amount=0.5); plt.clf(); plt.scatter(real(λ_dg), imag(λ_dg), label="DG"); plt.scatter(real(λ_sc), imag(λ_sc), label="SC"); plt.gca().set_aspect("equal"); plt.xlabel(L"\mathrm{Re}(\lambda)"); plt.ylabel(L"\mathrm{Im}(\lambda)"); plt.legend(loc="upper left")
# plt.savefig("spectra_dg_sc.png")
function compute_spectra(; angle=π/4, shock_capturing_amount=0.5)
  equations = LinearScalarAdvectionEquation2D(sin(angle), cos(angle))
  mesh = TreeMesh((-1.0, -1.0), (1.0, 1.0), initial_refinement_level=3,
                  n_cells_max=10^4, periodicity=true)

  # plain DG semidiscretization
  surface_flux = flux_lax_friedrichs
  solver_dg = DGSEM(polydeg=3, surface_flux=surface_flux)
  semi_dg = SemidiscretizationHyperbolic(mesh, equations,
    initial_condition_convergence_test, solver_dg)

  # shock capturing semidiscretization with fixed amount of FV
  indicator = IndicatorFixed(equations, solver_dg.basis,
                             alpha_fixed=shock_capturing_amount)
  solver_sc = DGSEM(solver_dg.basis, surface_flux,
    VolumeIntegralShockCapturingHG(indicator,
      volume_flux_dg=flux_central, volume_flux_fv=surface_flux))
  semi_sc = SemidiscretizationHyperbolic(mesh, equations,
    initial_condition_convergence_test, solver_sc)

  A_dg, b_dg = linear_structure(semi_dg)
  @assert iszero(b_dg)
  J = Matrix(A_dg)
  λ_dg = eigvals(J)

  A_sc, b_sc = linear_structure(semi_sc)
  @assert iszero(b_sc)
  J = Matrix(A_sc)
  λ_sc = eigvals(J)

  return λ_dg, λ_sc
end



# Plot the stability region of an RK method given by the Butcher tableau `tab`.
function plot_stability_region(tab::DiffEqBase.ExplicitRKTableau;
                               fig=nothing, plot_points=500, scale_factor=1.0)
  x = range(-1.45, 0.15, length=plot_points)
  y = range(-0.8, 0.8, length=plot_points)
  z = abs.(stability_region.((x' .+ im .* y) .* scale_factor, Ref(tab)))
  plt.figure(fig)
  plt.clf()
  plt.contourf(x, y, z, [0, 1], colors="gray")
  plt.xlabel(L"\mathrm{Re}(\lambda)")
  plt.ylabel(L"\mathrm{Im}(\lambda)")
  plt.gca().set_aspect("equal")
  return fig
end

# Compute the scaling factor so that all components of `factor * λ` are in the
# region of absolute stability of the RK method given by `tab`.
function compute_cfl_factor(λ, tab, scale_factor)
  tol_factor = 1.0e-7
  tol_value = 1.0e-14
  factor_min = 1.0e-8
  factor_max = 1.0e2
  factor = factor_max
  r = abs.(stability_region.(λ .* scale_factor .* factor_min, Ref(tab)))
  @assert all(<=(1 + tol_value), r)
  for _ in 1:10^4
    factor = 0.5 * (factor_min + factor_max)
    r .= abs.(stability_region.(λ .* scale_factor .* factor, Ref(tab)))
    if all(<=(1 + tol_value), r)
      factor_min = factor
    else
      factor_max = factor
    end

    factor_max - factor_min < tol_factor && break
  end
  return factor
end

# Plot the stability region of the RK method `alg` and the scaled DG spectra
# ` λ_dg, λ_sc` so that `λ_dg` is in the region of absolute stability of `alg`.
function plot_spectra_and_stability_region(alg, λ_dg, λ_sc;
                                           fig=nothing)
  A, b = deduce_Butcher_tableau(alg)
  tab = DiffEqBase.ExplicitRKTableau(A, vec(sum(A, dims=2)), b,
                                     OrdinaryDiffEq.alg_order(alg))

  # plot stability region scaled by effective number of stages
  scale_factor = length(b) - (OrdinaryDiffEq.isfsal(alg) ? 1 : 0)
  fig = plot_stability_region(tab; fig, scale_factor)

  # bisect CFL factor
  factor = compute_cfl_factor(λ_dg, tab, scale_factor)
  cfl_ratio_dg_sc = factor / compute_cfl_factor(λ_sc, tab, scale_factor)
  @info "Linear CFL factors" alg cfl_ratio_dg_sc factor

  # plot spectra on top
  kwargs = (color="black", linestyle="none", markersize=8, markeredgewidth=1.5)
  plt.plot(real(λ_dg) * factor, imag(λ_dg) * factor, label="DG", marker="x"; kwargs...)
  plt.plot(real(λ_sc) * factor, imag(λ_sc) * factor, label="SC", marker="."; kwargs...)
  plt.legend(loc="upper left")

  return fig
end



# Main function for this setup
function run_experiments()
  λ_dg, λ_sc = compute_spectra(angle=π/4, shock_capturing_amount=0.5)

  for alg in (SSPRK43(), BS3(), RDPK3SpFSAL35(), RDPK3SpFSAL49())
    plot_spectra_and_stability_region(alg, λ_dg, λ_sc; fig=plt.gcf())
    #plt.title(nameof(typeof(alg)))
    plt.savefig("spectra_DG_SC_$(nameof(typeof(alg))).png", bbox_inches="tight")
    plt.savefig("spectra_DG_SC_$(nameof(typeof(alg))).pdf", bbox_inches="tight")
  end
end
