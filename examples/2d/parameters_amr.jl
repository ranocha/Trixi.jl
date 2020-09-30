# TODO: Taal refactor, rename to
# - linear_advection_amr.jl
# - advection_amr.jl
# or something similar?

using OrdinaryDiffEq
using Trixi

###############################################################################
# semidiscretization of the linear advection equation

advectionvelocity = (1.0, 1.0)
# advectionvelocity = (0.2, -0.3)
equations = LinearScalarAdvectionEquation2D(advectionvelocity)

initial_conditions = initial_conditions_gauss

surface_flux = flux_lax_friedrichs
solver = DGSEM(3, surface_flux)

coordinates_min = (-5, -5)
coordinates_max = ( 5,  5)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=4,
                n_cells_max=30_000)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_conditions, solver)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
alive_callback = AliveCallback(analysis_interval=analysis_interval)
analysis_callback = AnalysisCallback(semi, analysis_interval=analysis_interval,
                                     extra_analysis_integrals=(entropy,))

save_solution = SaveSolutionCallback(solution_interval=100,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=:primitive)

# TODO: Taal, IO
# restart_interval = 10

amr_indicator = IndicatorThreeLevel(semi, IndicatorMax(semi),
                                    base_level=4,
                                    med_level=5, med_threshold=0.1,
                                    max_level=6, max_threshold=0.6)
amr_callback = AMRCallback(semi, amr_indicator,
                           interval=1, # TODO: Taal debug, interval=5
                           adapt_initial_conditions=true,
                           adapt_initial_conditions_only_refine=true)

stepsize_callback = StepsizeCallback(cfl=1.6)

# TODO: Taal decide, first AMR or save solution etc.
callbacks = CallbackSet(amr_callback, alive_callback);
# callbacks = CallbackSet(summary_callback, amr_callback, stepsize_callback, analysis_callback, save_solution, alive_callback); # TODO: Taal debug


###############################################################################
# run the simulation

function unstable_check(dt, u, semi, t)
  # u = Trixi.wrap_array(u_ode, semi)
  unstable_idx = findfirst(u -> abs(u) < 0.01, u)
  if unstable_idx === nothing
    return false
  end
  @warn "unstable_check" unstable_idx u[unstable_idx] t
  return true
end
# sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false), dt=stepsize_callback(ode), # TODO: Taal debug
sol = solve(ode, BS3(), unstable_check = unstable_check,
            save_everystep=false, callback=callbacks);
# summary_callback() # print the timer summary
nothing # TODO: Taal debug
