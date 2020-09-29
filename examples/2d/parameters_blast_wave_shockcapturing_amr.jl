# TODO: Taal refactor, rename to
# - euler_blast_wave_shockcapturing_amr.jl
# or something similar?

using OrdinaryDiffEq
using Trixi

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

initial_conditions = initial_conditions_blast_wave

surface_flux = flux_lax_friedrichs
volume_flux  = flux_chandrashekar
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

coordinates_min = (-2, -2)
coordinates_max = ( 2,  2)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=5, # TODO: Taal debug
                n_cells_max=10_000)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_conditions, solver)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 12.5)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
alive_callback = AliveCallback(analysis_interval=analysis_interval)
analysis_callback = AnalysisCallback(semi, analysis_interval=analysis_interval)

# indicator_amr = IndicatorMax(semi, variable=density)
# amr_indicator = IndicatorThreeLevel(semi, indicator_amr,
#                                     base_level=4,
#                                     max_level=6, max_threshold=0.5)
indicator_amr = IndicatorHennemannGassner(semi, # TODO: Taal debug
                                          alpha_max=0.5,
                                          alpha_min=0.001,
                                          alpha_smooth=true,
                                          variable=density_pressure)
amr_indicator = IndicatorThreeLevel(semi, indicator_amr,
                                    base_level=4,
                                    max_level=6, max_threshold=0.1)
amr_callback = AMRCallback(semi, amr_indicator,
                           interval=5, # TODO: Taal debug
                           adapt_initial_conditions=true, # TODO: Taal debug
                           adapt_initial_conditions_only_refine=true) # TODO: Taal debug
amr_callback(ode) # TODO: Taal debug

stepsize_callback = StepsizeCallback(cfl=0.25) # TODO: Taal debug cfl=1.0

save_solution = SaveSolutionCallback(solution_interval=10, # TODO: Taal debug, 100
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=:primitive)

# callbacks = CallbackSet(summary_callback, stepsize_callback, analysis_callback, save_solution, alive_callback)
callbacks = CallbackSet(summary_callback, amr_callback, stepsize_callback, analysis_callback, save_solution, alive_callback) # TODO: Taal debug


###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false), dt=stepsize_callback(ode),
            save_everystep=false, callback=callbacks);
# sol = solve(ode, BS3(), #adaptive=false, dt=stepsize_callback(ode), # TODO: Taal debug
#             save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary
