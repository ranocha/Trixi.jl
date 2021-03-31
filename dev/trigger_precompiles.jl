
#=
Run
```julia
using SnoopCompile
tinf = @snoopi_deep include("dev/trigger_precompiles.jl")
ttot, pcs = SnoopCompile.parcel(tinf);
SnoopCompile.write("dev/precompile", pcs)

```
or
```julia
using SnoopCompile
inf_timing = @snoopi tmin=0.001 include("test/test_examples_2d.jl") # TODO: Replace DiffEq by our method, seems to be better somehow...
pc = SnoopCompile.parcel(inf_timing)
SnoopCompile.write("dev/precompile", pc)

```
or
```julia
using SnoopCompile
SnoopCompile.@snoopc "/tmp/snoop_compile.log" include("dev/trigger_precompiles.jl")
data = SnoopCompile.read("/tmp/snoop_compile.log")
pc = SnoopCompile.parcel(reverse!(data[2]))
SnoopCompile.write("dev/precompile", pc)

```
Then, add
```julia
include("../dev/precompile/precompile_Trixi.jl")
_precompile_()
```
to Trixi/src/Trixi.jl. Remember to remove these lines before generating new
precompile statements since these will override the existing ones!

In preliminary tests, `@snoopi` seems to result in better performance.

See
- https://timholy.github.io/SnoopCompile.jl/stable/snoopi/
- https://timholy.github.io/SnoopCompile.jl/stable/snoopc/
=#

# Below is a basically a copy of elixir_advection_amr.jl


using Trixi
using OrdinaryDiffEq

###############################################################################
# semidiscretization of the linear advection equation

advectionvelocity = (1.0, 1.0)
# advectionvelocity = (0.2, -0.3)
equations = LinearScalarAdvectionEquation2D(advectionvelocity)

initial_condition = initial_condition_gauss

surface_flux = flux_lax_friedrichs
basis = LobattoLegendreBasis(3)
solver = DGSEM(basis, surface_flux)

coordinates_min = (-5, -5)
coordinates_max = ( 5,  5)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=4,
                n_cells_max=30_000)


semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

amr_controller = ControllerThreeLevel(semi, IndicatorMax(semi, variable=first),
                                      base_level=4,
                                      med_level=5, med_threshold=0.1,
                                      max_level=6, max_threshold=0.6)
amr_callback = AMRCallback(semi, amr_controller,
                           interval=5,
                           adapt_initial_condition=true,
                           adapt_initial_condition_only_refine=true)

stepsize_callback = StepsizeCallback(cfl=1.6)

save_solution = SaveSolutionCallback(interval=100,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=:primitive)

save_restart = SaveRestartCallback(interval=100,
                                   save_final_restart=true)

analysis_interval = 100
alive_callback = AliveCallback(analysis_interval=analysis_interval)
analysis_callback = AnalysisCallback(semi,
                                     interval=analysis_interval,
                                     extra_analysis_integrals=(entropy,))

# TODO: Taal decide, first AMR or save solution etc.
callbacks = CallbackSet(summary_callback, amr_callback, stepsize_callback,
                        save_restart, save_solution,
                        analysis_callback, alive_callback);


###############################################################################
# run the simulation

u_ode = copy(ode.u0)
du_ode = similar(u_ode)
Trixi.rhs!(du_ode, u_ode, semi, first(tspan))

# sol = Trixi.solve(ode, Trixi.SimpleAlgorithm2N45(),
#                   dt=stepsize_callback(ode), # solve needs some value here but it will be overwritten by the stepsize_callback
#                   save_everystep=false, callback=callbacks);
summary_callback() # print the timer summary

Trixi.print_timer(Trixi.timer())

# It seems to be better to not include this file when using @snoopi
# trixi_include(joinpath(examples_dir(), "2d", "elixir_advection_amr.jl"))


# u_ode = copy(sol.u[end])
# du_ode = similar(u_ode)
# Trixi.rhs!(du_ode, u_ode, semi, sol.t[end])
