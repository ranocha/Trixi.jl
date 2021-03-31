# julia --check-bounds=no dev/time_compilation.jl  | egrep "seconds"

using Printf
@time using Trixi
@time using OrdinaryDiffEq

@time begin
  let
    println("-"^80)
    println("Linear advection 1D")
    @printf("%-18s: ", "equations")
    @time equations = LinearScalarAdvectionEquation1D(1.0)
    @printf("%-18s: ", "solver")
    @time solver = DGSEM(3, flux_lax_friedrichs)
    @printf("%-18s: ", "mesh")
    @time mesh = TreeMesh(-1.0, 1.0, initial_refinement_level=4, n_cells_max=30_000)
    @printf("%-18s: ", "semi")
    @time semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_convergence_test, solver)
    @printf("%-18s: ", "ode")
    @time ode = semidiscretize(semi, (0.0, 1.0))
    @printf("%-18s: ", "summary_callback")
    @time summary_callback = SummaryCallback()
    @printf("%-18s: ", "analysis_callback")
    @time analysis_callback = AnalysisCallback(semi, interval=100)
    @printf("%-18s: ", "save_solution")
    @time save_solution = SaveSolutionCallback(interval=100, solution_variables=cons2prim)
    @printf("%-18s: ", "stepsize_callback")
    @time stepsize_callback = StepsizeCallback(cfl=1.6)
    @printf("%-18s: ", "callbacks")
    # @time callbacks = CallbackSet(summary_callback, analysis_callback, save_solution, stepsize_callback)
    @time callbacks = CallbackSet(summary_callback, save_solution, stepsize_callback)
    # @time callbacks = CallbackSet(stepsize_callback)
    @printf("%-18s: ", "solve")
    @time sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false), dt=1.0,  save_everystep=false, callback=callbacks)
    @printf("%-18s: ", "summary_callback()")
    @time summary_callback()
  end

  let
    println("-"^80)
    println("Compressible Euler equations 2D")
    @printf("%-18s: ", "equations")
    @time equations = CompressibleEulerEquations2D(1.4)
    @printf("%-18s: ", "solver")
    @time solver = DGSEM(4, FluxPlusDissipation(flux_ranocha, DissipationLocalLaxFriedrichs()), VolumeIntegralFluxDifferencing(flux_ranocha))
    @printf("%-18s: ", "mesh")
    @time mesh = TreeMesh((-1.0, -1.0), (1.0, 1.0), initial_refinement_level=3, n_cells_max=30_000)
    @printf("%-18s: ", "semi")
    @time semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_convergence_test, solver)
    @printf("%-18s: ", "ode")
    @time ode = semidiscretize(semi, (0.0, 1.0))
    @printf("%-18s: ", "summary_callback")
    @time summary_callback = SummaryCallback()
    @printf("%-18s: ", "analysis_callback")
    @time analysis_callback = AnalysisCallback(semi, interval=100)
    @printf("%-18s: ", "save_solution")
    @time save_solution = SaveSolutionCallback(interval=100, solution_variables=cons2prim)
    @printf("%-18s: ", "stepsize_callback")
    @time stepsize_callback = StepsizeCallback(cfl=1.6)
    @printf("%-18s: ", "callbacks")
    # @time callbacks = CallbackSet(summary_callback, analysis_callback, save_solution, stepsize_callback)
    @time callbacks = CallbackSet(summary_callback, save_solution, stepsize_callback)
    # @time callbacks = CallbackSet(stepsize_callback)
    @printf("%-18s: ", "solve")
    @time sol = solve(ode, BS3(), save_everystep=false, callback=callbacks)
    @printf("%-18s: ", "summary_callback()")
    @time summary_callback()
  end
end # @time
