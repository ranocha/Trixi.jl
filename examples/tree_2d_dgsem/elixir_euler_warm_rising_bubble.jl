
using OrdinaryDiffEq
using Trixi

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

"""
    initial_condition_warm_rising_bubble(x, t, equations::CompressibleEulerEquations2D)

The warm rising bubble test ("test case 2") taken from
- Francis X. Giraldo, M. Restelli (2008)
  A study of spectral element and discontinuous Galerkin methods for the Navier-Stokes
  equations in nonhydrostatic mesoscale atmospheric modeling: Equation sets and test cases
  [DOI: 10.1016/j.jcp.2007.12.009](https://doi.org/10.1016/j.jcp.2007.12.009)

See also
- Maciej Waruszewski, Jeremy E. Kozdon, Lucas C. Wilcox, Thomas H. Gibson,
  Francis X. Giraldo (2022)
  Entropy stable discontinuous Galerkin methods for balance laws in
  non-conservative form: Applications to the Euler equations with gravity
  [DOI: 10.1016/j.jcp.2022.111507](https://doi.org/10.1016/j.jcp.2022.111507)
"""
function initial_condition_warm_rising_bubble(x, t, equations::CompressibleEulerEquations2D)
  g = 9.81 # gravitational constant
  γ = equations.gamma

  cv = 719.0
  cp = γ * cv
  R = cp - cv

  θref = 300.0 # Kelvin
  p0 = 1.0e5   # Pascal
  xc = 500.0   # meter
  zc = 350.0   # meter
  rc = 250.0   # meter

  r = sqrt((x[1] - xc)^2 + (x[2] - zc)^2)
  θ = θref
  if r <= rc
    θ += 0.5 * (1 + cospi(r / rc))
  end

  π_exner = 1 - g / (cp * θref) * x[2]
  ϱ = p0 / (R * θ) * π_exner^(cv / R)

  ϱv1 = 0.0
  ϱv2 = 0.0

  Φ = g * x[2]
  T = θ * π_exner
  ϱe = ϱ * (cv * T + Φ)

  return SVector(ϱ, ϱv1, ϱv2, ϱe)
end

function cons2pot(u, equations::CompressibleEulerEquations2D)
  g = 9.81 # gravitational constant
  γ = equations.gamma

  cv = 719.0
  cp = γ * cv
  R = cp - cv

  rho, v1, v2, p = cons2prim(u, equations)

  Φ = constants(law).grav * z

  p = EulerGravity.pressure(law, ρ, ρu⃗, ρe, Φ)
  p0 = 1e5
  T = p / (R_d * ρ)
  θ = T * (p0 / p)^(1 - 1 / constants(law).γ)

  return (rho, rho * v1, rho * v2, θ)
end

surface_flux = flux_lax_friedrichs
volume_flux  = flux_ranocha_turbo
basis = LobattoLegendreBasis(3)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max = 0.5,
                                         alpha_min = 0.001,
                                         alpha_smooth = true,
                                         variable = density_pressure)
# volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
#                                                  volume_flux_dg = volume_flux,
#                                                  volume_flux_fv = surface_flux)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (1000.0, 1000.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 5,
                n_cells_max = 10_000,
                periodicity = (true, false))

@inline function boundary_condition_noflux(u_inner, orientation::Integer,
                                           direction, x, t,
                                           surface_flux_function,
                                           equations::CompressibleEulerEquations2D)
    rho, v1, v2, p = cons2prim(u_inner, equations)

    if orientation == 1
      v1 = zero(v1)
    else # orientation == 2
      v2 = zero(v2)
    end

    u_boundary = prim2cons(SVector(rho, v1, v2, p), equations)

    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
      flux = surface_flux_function(u_inner, u_boundary, orientation,
                                   equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation,
                                    equations)
    end

    return flux
end

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)
                      #  y_neg = boundary_condition_noflux,
                      #  y_pos = boundary_condition_noflux)

function source_terms(u, x, t, equations::CompressibleEulerEquations2D)
  g = 9.81 # gravitational constant
  rho = u[1]
  rho_v2 = u[3]

  return SVector(0.0, 0.0, -g * rho, -g * rho_v2)
end

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_warm_rising_bubble, solver;
                                    boundary_conditions, source_terms)


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 700.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 500
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution)


###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(thread = OrdinaryDiffEq.True());
            controller = PIDController(0.55, -0.27, 0.05),
            ode_default_options()..., callback = callbacks);
summary_callback() # print the timer summary


#=

julia> begin
       x = view(semi.cache.elements.node_coordinates, 1, :, :, :)
       y = view(semi.cache.elements.node_coordinates, 2, :, :, :)
       u = Trixi.wrap_array_native(sol.u[1], semi)
       rho    = view(u, 1, :, :, :)
       rho_v1 = view(u, 2, :, :, :)
       rho_v2 = view(u, 3, :, :, :)
       rho_e  = view(u, 4, :, :, :)
       g = 9.81
       γ = equations.gamma
       cv = 719.0
       cp = γ * cv
       R = cp - cv
       Φ = g * y
       p = @. (γ - 1) * (rho_e - rho * Φ - 0.5 * (rho_v1^2 + rho_v2^2) / rho)
       p0 = 1.0e5
       T = @. p / (R * rho)
       θ = @. T * (p0 / p)^(1 - 1 / γ)
       θ0 = 300.0
       θref = θ0 # TODO
       pd = ScalarPlotData2D(θ .- θref, semi)
       plot(pd)
       end


=#
