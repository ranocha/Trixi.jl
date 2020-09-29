
mutable struct StepsizeCallback{RealT}
  cfl_number::RealT
end


function Base.show(io::IO, cb::DiscreteCallback{Condition,Affect!}) where {Condition, Affect!<:StepsizeCallback}
  stepsize_callback = cb.affect!
  @unpack cfl_number = stepsize_callback
  print(io, "StepsizeCallback with CFL number ", cfl_number)
end
# TODO: Taal bikeshedding, implement a method with more information and the signature
# function Base.show(io::IO, ::MIME"text/plain", cb::DiscreteCallback{Condition,Affect!}) where {Condition, Affect!<:StepsizeCallback}
# end


function StepsizeCallback(; cfl::Real=1.0)
  # when is the callback activated
  condition = (u, t, integrator) -> true

  stepsize_callback = StepsizeCallback(cfl)

  DiscreteCallback(condition, stepsize_callback,
                   save_positions=(false,false),
                   initialize=initialize!)
end


function initialize!(cb::DiscreteCallback{Condition,Affect!}, u, t, integrator) where {Condition, Affect!<:StepsizeCallback}
  cb.affect!(integrator)
end


@inline function (stepsize_callback::StepsizeCallback)(integrator)
  # TODO: Taal decide, shall we set the time step even if the integrator is adaptive?
  if !integrator.opts.adaptive
    t = integrator.t
    u_ode = integrator.u
    semi = integrator.p
    mesh, equations, solver, cache = mesh_equations_solver_cache(semi)
    @unpack cfl_number = stepsize_callback
    u = wrap_array(u_ode, mesh, equations, solver, cache)

    dt = #=@timeit_debug timer() "calculate dt"=# cfl_number * max_dt(u, t, mesh,
                                                                  have_constant_speed(equations), equations,
                                                                  solver, cache)
    set_proposed_dt!(integrator, dt)
    integrator.opts.dtmax = dt
    integrator.dtcache = dt
  end

  u_modified!(integrator, false)
  return nothing
end


function (cb::DiscreteCallback{Condition,Affect!})(ode::ODEProblem) where {Condition, Affect!<:StepsizeCallback}
  stepsize_callback = cb.affect!
  @unpack cfl_number = stepsize_callback
  u_ode = ode.u0
  t = first(ode.tspan)
  semi = ode.p
  mesh, equations, solver, cache = mesh_equations_solver_cache(semi)
  u = wrap_array(u_ode, mesh, equations, solver, cache)

  return cfl_number * max_dt(u, t, mesh, have_constant_speed(equations), equations, solver, cache)
end


include("stepsize_dg2d.jl")
