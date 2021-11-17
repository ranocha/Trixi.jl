# TODO: This probably needs to live somewhere else.

# The following code can be used for micro-benchmarks while developing this part.
# Once everything is settled down, the new kernel should become a specialized
# method of the old kernel, of course.

#=
using BenchmarkTools, StructArrays, StrideArrays
using Revise, Trixi

redirect_stdout(devnull) do
  trixi_include("examples/tree_3d_dgsem/elixir_euler_ec.jl",
    volume_flux=flux_shima_etal, surface_flux=flux_shima_etal, tspan=(0.0, 0.0),
    polydeg=3)
end

u_ode = copy(sol.u[end]); du_ode = similar(u_ode);
# LoopVectorization.jl takes forever to compile the flux differencing kernels
# if it needs to specialize for high polynomial degrees at compile time via
# `indices(...)`, e.g., `polydeg = 7`
u = Trixi.wrap_array(u_ode, semi); du = Trixi.wrap_array(du_ode, semi);
u_prim_ode = zero(u_ode); u_prim = Trixi.wrap_array(u_prim_ode, semi); # TODO: element-wise
# u = PtrArray(pointer(u_ode), (nvariables(equations), ntuple(_ -> nnodes(solver), ndims(mesh))...,
#   nelements(solver, semi.cache)));
# du = PtrArray(pointer(du_ode), (nvariables(equations), ntuple(_ -> nnodes(solver), ndims(mesh))...,
#   nelements(solver, semi.cache)));
Trixi.rhs!(du_ode, u_ode, semi, 0.0)

u_aos_ode = reinterpret(SVector{nvariables(equations), eltype(u_ode)}, u_ode);
u_soa_ode = StructArray(u_aos_ode); du_soa_ode = copy(u_soa_ode);
# u_soa = StructArray{SVector{nvariables(equations), eltype(u_ode)}}(
#   map(u_ode -> PtrArray(pointer(u_ode), (ntuple(_ -> StaticInt(nnodes(solver)), ndims(mesh))...,
#     nelements(solver, semi.cache))), StructArrays.components(u_soa_ode)));
# du_soa = StructArray{SVector{nvariables(equations), eltype(u_ode)}}(
#   map(u_ode -> PtrArray(pointer(u_ode), (ntuple(_ -> StaticInt(nnodes(solver)), ndims(mesh))...,
#     nelements(solver, semi.cache))), StructArrays.components(du_soa_ode)));
u_soa = StructArray{SVector{nvariables(equations), eltype(u_ode)}}(
  map(u_ode -> PtrArray(pointer(u_ode), (ntuple(_ -> nnodes(solver), ndims(mesh))...,
    nelements(solver, semi.cache))), StructArrays.components(u_soa_ode)));
du_soa = StructArray{SVector{nvariables(equations), eltype(u_ode)}}(
  map(u_ode -> PtrArray(pointer(u_ode), (ntuple(_ -> nnodes(solver), ndims(mesh))...,
    nelements(solver, semi.cache))), StructArrays.components(du_soa_ode)));

@benchmark Trixi.split_form_kernel!(
  $du, $u, $(Ref(1))[], $(semi.mesh),
  $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
  $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))

@benchmark Trixi.flux_differencing_kernel!(
  $du, $u, $(Ref(1))[], $(semi.mesh),
  $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
  $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))

@benchmark Trixi.flux_differencing_kernel!(
  $du, $u_prim, $u, $(Ref(1))[], $(semi.mesh),
  $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
  $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))

@benchmark Trixi.flux_differencing_kernel!(
  $du_soa, $u_soa, $(Ref(1))[], $(semi.mesh),
  $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
  $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))


function run_benchmarks(polydegs)
  baseline = Vector{Float64}()
  turbo_on = Vector{Float64}()
  prim_var = Vector{Float64}()
  rearranged = Vector{Float64}()

  for polydeg in polydegs
    @show polydeg

    redirect_stdout(devnull) do
      trixi_include("examples/tree_3d_dgsem/elixir_euler_ec.jl",
        volume_flux=flux_shima_etal, surface_flux=flux_shima_etal,
        tspan=(0.0, 0.0), initial_refinement_level=0,
        polydeg=polydeg)
    end
    u_ode = copy(ode.u0)
    du_ode = zero(u_ode)
    u_prim_ode = zero(u_ode)

    GC.@preserve u_ode du_ode u_prim_ode begin
      u = Trixi.wrap_array(u_ode, semi)
      du = Trixi.wrap_array(du_ode, semi)
      u_prim = Trixi.wrap_array(u_prim_ode, semi)

      result = @benchmark Trixi.split_form_kernel!(
        $du, $u, $(Ref(1))[], $(semi.mesh),
        $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
        $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))
      push!(baseline, time(minimum(result)))

      result = @benchmark Trixi.flux_differencing_kernel!(
        $du, $u, $(Ref(1))[], $(semi.mesh),
        $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
        $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))
      push!(turbo_on, time(minimum(result)))

      result = @benchmark Trixi.flux_differencing_kernel!(
        $du, $u_prim, $u, $(Ref(1))[], $(semi.mesh),
        $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
        $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))
      push!(prim_var, time(minimum(result)))

      result = @benchmark Trixi.flux_differencing_kernel_rearranged!(
        $du, $u_prim, $u, $(Ref(1))[], $(semi.mesh),
        $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
        $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))
      push!(rearranged, time(minimum(result)))
    end
  end

  hcat(polydegs, baseline,
       turbo_on, turbo_on ./ baseline,
       prim_var, prim_var ./ baseline,
       rearranged, rearranged ./ baseline)
end

results = run_benchmarks(1:15)
# 15×8 Matrix{Float64}:
#   1.0     130.454     196.455  1.50593      157.305  1.20583       88.0729  0.675126
#   2.0     651.877     669.094  1.02641      468.454  0.718623     395.095   0.606087
#   3.0    2078.22     1464.4    0.704641    1130.5    0.543975    1135.5     0.54638
#   4.0    5075.0      5163.0    1.01734     3633.0    0.715862    3504.0     0.690443
#   5.0   10604.0      8733.0    0.823557    6325.4    0.596511    6803.6     0.641607
#   6.0   19774.0     14073.0    0.711692    9989.0    0.505158   12532.0     0.633762
#   7.0   33792.0     20634.0    0.610618   14936.0    0.441998   19263.0     0.570046
#   8.0   55387.0     41995.0    0.75821    29538.0    0.533302   37043.0     0.668803
#   9.0   85167.0     57706.0    0.677563   39303.0    0.461482   55207.0     0.648221
#  10.0  124718.0     77033.0    0.617657   54412.0    0.43628    82762.0     0.663593
#  11.0  178912.0     99571.0    0.556536   70343.0    0.393171  109271.0     0.610753
#  12.0  250083.0    163108.0    0.652215  111759.0    0.446888  165991.0     0.663744
#  13.0  346647.0    204391.0    0.589623  146411.0    0.422363  225718.0     0.651147
#  14.0  443217.0    253755.0    0.57253   177370.0    0.400188  312799.0     0.705747
#  15.0  579960.0    306063.0    0.527731  218419.0    0.37661   372407.0     0.642125

using Plots
fig_abs = scatter(results[:, 1], results[:, 2], label="Baseline",
  xguide="Polynomial degree", yguide="Runtime", yscale=:log10, legend=:bottomright);
  scatter!(results[:, 1], results[:, 3], label="Turbo");
  scatter!(results[:, 1], results[:, 5], label="Primitive variables");
  scatter!(results[:, 1], results[:, 7], label="Rearranged")
fig_rel = scatter(results[:, 1], results[:, 4], label="Turbo",
  xguide="Polynomial degree", yguide="Relative runtime");
  scatter!(results[:, 1], results[:, 6], label="Primitive variables");
  scatter!(results[:, 1], results[:, 8], label="Rearranged")
plot(fig_abs, fig_rel, layout=(2,1))
=#


@inline function flux_differencing_kernel!(du, u,
                                           element, mesh::TreeMesh{3},
                                           nonconservative_terms::Val{false},
                                           equations::CompressibleEulerEquations3D,
                                           volume_flux::typeof(flux_shima_etal),
                                           dg::DGSEM, cache,
                                           alpha=true)
  # true * [some floating point value] == [exactly the same floating point value]
  # This can (hopefully) be optimized away due to constant propagation.
  @unpack derivative_split, derivative_split_transpose = dg.basis

  # x direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du, derivative_split), (2, 2, 1))
    rho_ll    = u[1, i, j, k, element]
    rho_v1_ll = u[2, i, j, k, element]
    rho_v2_ll = u[3, i, j, k, element]
    rho_v3_ll = u[4, i, j, k, element]
    rho_e_ll  = u[5, i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for ii in eachnode(dg)
    # for ii in indices(u, 2,)
    # for ii in indices((u, derivative_split), (2, 2))
      rho_rr    = u[1, ii, j, k, element]
      rho_v1_rr = u[2, ii, j, k, element]
      rho_v2_rr = u[3, ii, j, k, element]
      rho_v3_rr = u[4, ii, j, k, element]
      rho_e_rr  = u[5, ii, j, k, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv1_avg = 0.5 * (p_ll * v1_rr + p_rr * v1_ll)

      f1 = rho_avg * v1_avg
      f2 = f1 * v1_avg + p_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg
      f5 = p_avg * v1_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv1_avg

      factor = alpha * derivative_split[i, ii]
      # factor = alpha * derivative_split_transpose[ii, i]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end

  # y direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
  # @turbo for k in indices((u, du), 4), j in indices((u, du, derivative_split), (3, 3, 1)), i in indices((u, du), 2)
    rho_ll    = u[1, i, j, k, element]
    rho_v1_ll = u[2, i, j, k, element]
    rho_v2_ll = u[3, i, j, k, element]
    rho_v3_ll = u[4, i, j, k, element]
    rho_e_ll  = u[5, i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for jj in eachnode(dg)
    # for jj in indices(u, 3)
    # for jj in indices((u, derivative_split), (3, 2))
      rho_rr    = u[1, i, jj, k, element]
      rho_v1_rr = u[2, i, jj, k, element]
      rho_v2_rr = u[3, i, jj, k, element]
      rho_v3_rr = u[4, i, jj, k, element]
      rho_e_rr  = u[5, i, jj, k, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv2_avg = 0.5 * (p_ll * v2_rr + p_rr * v2_ll)

      f1 = rho_avg * v2_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg + p_avg
      f4 = f1 * v3_avg
      f5 = p_avg*v2_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv2_avg

      factor = alpha * derivative_split[j, jj]
      # factor = alpha * derivative_split_transpose[jj, j]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end

  # z direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
  # @turbo for k in indices((u, du, derivative_split), (4, 4, 1)), j in indices((u, du), 3), i in indices((u, du), 2)
    rho_ll    = u[1, i, j, k, element]
    rho_v1_ll = u[2, i, j, k, element]
    rho_v2_ll = u[3, i, j, k, element]
    rho_v3_ll = u[4, i, j, k, element]
    rho_e_ll  = u[5, i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for kk in eachnode(dg)
    # for kk in indices(u, 4)
    # for kk in indices((u, derivative_split), (4, 2))
      rho_rr    = u[1, i, j, kk, element]
      rho_v1_rr = u[2, i, j, kk, element]
      rho_v2_rr = u[3, i, j, kk, element]
      rho_v3_rr = u[4, i, j, kk, element]
      rho_e_rr  = u[5, i, j, kk, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv3_avg = 0.5 * (p_ll * v3_rr + p_rr * v3_ll)

      f1 = rho_avg * v3_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg + p_avg
      f5 = p_avg*v3_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv3_avg

      factor = alpha * derivative_split[k, kk]
      # factor = alpha * derivative_split_transpose[kk, k]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end
end


# pre-compute primitive variables
@inline function flux_differencing_kernel!(du, u_prim, u_cons,
                                           element, mesh::TreeMesh{3},
                                           nonconservative_terms::Val{false},
                                           equations::CompressibleEulerEquations3D,
                                           volume_flux::typeof(flux_shima_etal),
                                           dg::DGSEM, cache,
                                           alpha=true)
  # true * [some floating point value] == [exactly the same floating point value]
  # This can (hopefully) be optimized away due to constant propagation.
  @unpack derivative_split, derivative_split_transpose = dg.basis

  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    rho    = u_cons[1, i, j, k, element]
    rho_v1 = u_cons[2, i, j, k, element]
    rho_v2 = u_cons[3, i, j, k, element]
    rho_v3 = u_cons[4, i, j, k, element]
    rho_e  = u_cons[5, i, j, k, element]

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2 + rho_v3 * v3))

    u_prim[1, i, j, k, element] = rho
    u_prim[2, i, j, k, element] = v1
    u_prim[3, i, j, k, element] = v2
    u_prim[4, i, j, k, element] = v3
    u_prim[5, i, j, k, element] = p
  end

  # x direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du, derivative_split), (2, 2, 1))
    rho_ll = u_prim[1, i, j, k, element]
    v1_ll  = u_prim[2, i, j, k, element]
    v2_ll  = u_prim[3, i, j, k, element]
    v3_ll  = u_prim[4, i, j, k, element]
    p_ll   = u_prim[5, i, j, k, element]

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for ii in eachnode(dg)
    # for ii in indices(u, 2,)
    # for ii in indices((u, derivative_split), (2, 2))
      rho_rr = u_prim[1, ii, j, k, element]
      v1_rr  = u_prim[2, ii, j, k, element]
      v2_rr  = u_prim[3, ii, j, k, element]
      v3_rr  = u_prim[4, ii, j, k, element]
      p_rr   = u_prim[5, ii, j, k, element]

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv1_avg = 0.5 * (p_ll * v1_rr + p_rr * v1_ll)

      f1 = rho_avg * v1_avg
      f2 = f1 * v1_avg + p_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg
      f5 = p_avg * v1_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv1_avg

      factor = alpha * derivative_split[i, ii]
      # factor = alpha * derivative_split_transpose[ii, i]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end

  # y direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
  # @turbo for k in indices((u, du), 4), j in indices((u, du, derivative_split), (3, 3, 1)), i in indices((u, du), 2)
    rho_ll = u_prim[1, i, j, k, element]
    v1_ll  = u_prim[2, i, j, k, element]
    v2_ll  = u_prim[3, i, j, k, element]
    v3_ll  = u_prim[4, i, j, k, element]
    p_ll   = u_prim[5, i, j, k, element]

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for jj in eachnode(dg)
    # for jj in indices(u, 3)
    # for jj in indices((u, derivative_split), (3, 2))
      rho_rr = u_prim[1, i, jj, k, element]
      v1_rr  = u_prim[2, i, jj, k, element]
      v2_rr  = u_prim[3, i, jj, k, element]
      v3_rr  = u_prim[4, i, jj, k, element]
      p_rr   = u_prim[5, i, jj, k, element]

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv2_avg = 0.5 * (p_ll * v2_rr + p_rr * v2_ll)

      f1 = rho_avg * v2_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg + p_avg
      f4 = f1 * v3_avg
      f5 = p_avg*v2_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv2_avg

      factor = alpha * derivative_split[j, jj]
      # factor = alpha * derivative_split_transpose[jj, j]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end

  # z direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
  # @turbo for k in indices((u, du, derivative_split), (4, 4, 1)), j in indices((u, du), 3), i in indices((u, du), 2)
    rho_ll = u_prim[1, i, j, k, element]
    v1_ll  = u_prim[2, i, j, k, element]
    v2_ll  = u_prim[3, i, j, k, element]
    v3_ll  = u_prim[4, i, j, k, element]
    p_ll   = u_prim[5, i, j, k, element]

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for kk in eachnode(dg)
    # for kk in indices(u, 4)
    # for kk in indices((u, derivative_split), (4, 2))
      rho_rr = u_prim[1, i, j, kk, element]
      v1_rr  = u_prim[2, i, j, kk, element]
      v2_rr  = u_prim[3, i, j, kk, element]
      v3_rr  = u_prim[4, i, j, kk, element]
      p_rr   = u_prim[5, i, j, kk, element]

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv3_avg = 0.5 * (p_ll * v3_rr + p_rr * v3_ll)

      f1 = rho_avg * v3_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg + p_avg
      f5 = p_avg*v3_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv3_avg

      factor = alpha * derivative_split[k, kk]
      # factor = alpha * derivative_split_transpose[kk, k]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end
end


# pre-compute primitive variables and uses triangular loops without turbo
@inline function flux_differencing_kernel_rearranged!(du, u_prim, u_cons,
                                                      element, mesh::TreeMesh{3},
                                                      nonconservative_terms::Val{false},
                                                      equations::CompressibleEulerEquations3D,
                                                      volume_flux::typeof(flux_shima_etal),
                                                      dg::DGSEM, cache,
                                                      alpha=true)
  # true * [some floating point value] == [exactly the same floating point value]
  # This can (hopefully) be optimized away due to constant propagation.
  @unpack derivative_split, derivative_split_transpose = dg.basis

  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    rho    = u_cons[1, i, j, k, element]
    rho_v1 = u_cons[2, i, j, k, element]
    rho_v2 = u_cons[3, i, j, k, element]
    rho_v3 = u_cons[4, i, j, k, element]
    rho_e  = u_cons[5, i, j, k, element]

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2 + rho_v3 * v3))

    u_prim[1, i, j, k, element] = rho
    u_prim[2, i, j, k, element] = v1
    u_prim[3, i, j, k, element] = v2
    u_prim[4, i, j, k, element] = v3
    u_prim[5, i, j, k, element] = p
  end

  # x direction
  for i in eachnode(dg), ii in (i+1):nnodes(dg)
    @turbo for k in eachnode(dg), j in eachnode(dg)
      rho_ll = u_prim[1, i, j, k, element]
      v1_ll  = u_prim[2, i, j, k, element]
      v2_ll  = u_prim[3, i, j, k, element]
      v3_ll  = u_prim[4, i, j, k, element]
      p_ll   = u_prim[5, i, j, k, element]

      rho_rr = u_prim[1, ii, j, k, element]
      v1_rr  = u_prim[2, ii, j, k, element]
      v2_rr  = u_prim[3, ii, j, k, element]
      v3_rr  = u_prim[4, ii, j, k, element]
      p_rr   = u_prim[5, ii, j, k, element]

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv1_avg = 0.5 * (p_ll * v1_rr + p_rr * v1_ll)

      f1 = rho_avg * v1_avg
      f2 = f1 * v1_avg + p_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg
      f5 = p_avg * v1_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv1_avg

      factor = alpha * derivative_split[i, ii]
      du[1, i, j, k, element] += factor * f1
      du[2, i, j, k, element] += factor * f2
      du[3, i, j, k, element] += factor * f3
      du[4, i, j, k, element] += factor * f4
      du[5, i, j, k, element] += factor * f5

      factor = alpha * derivative_split[ii, i]
      du[1, ii, j, k, element] += factor * f1
      du[2, ii, j, k, element] += factor * f2
      du[3, ii, j, k, element] += factor * f3
      du[4, ii, j, k, element] += factor * f4
      du[5, ii, j, k, element] += factor * f5
    end
  end

  # y direction
  for j in eachnode(dg), jj in (j+1):nnodes(dg)
    @turbo for k in eachnode(dg), i in eachnode(dg)
      rho_ll = u_prim[1, i, j, k, element]
      v1_ll  = u_prim[2, i, j, k, element]
      v2_ll  = u_prim[3, i, j, k, element]
      v3_ll  = u_prim[4, i, j, k, element]
      p_ll   = u_prim[5, i, j, k, element]

      rho_rr = u_prim[1, i, jj, k, element]
      v1_rr  = u_prim[2, i, jj, k, element]
      v2_rr  = u_prim[3, i, jj, k, element]
      v3_rr  = u_prim[4, i, jj, k, element]
      p_rr   = u_prim[5, i, jj, k, element]

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv2_avg = 0.5 * (p_ll * v2_rr + p_rr * v2_ll)

      f1 = rho_avg * v2_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg + p_avg
      f4 = f1 * v3_avg
      f5 = p_avg*v2_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv2_avg

      factor = alpha * derivative_split[j, jj]
      du[1, i, j, k, element] += factor * f1
      du[2, i, j, k, element] += factor * f2
      du[3, i, j, k, element] += factor * f3
      du[4, i, j, k, element] += factor * f4
      du[5, i, j, k, element] += factor * f5

      factor = alpha * derivative_split[jj, j]
      du[1, i, jj, k, element] += factor * f1
      du[2, i, jj, k, element] += factor * f2
      du[3, i, jj, k, element] += factor * f3
      du[4, i, jj, k, element] += factor * f4
      du[5, i, jj, k, element] += factor * f5
    end
  end

  # z direction
  for k in eachnode(dg), kk in (k+1):nnodes(dg)
    @turbo for j in eachnode(dg), i in eachnode(dg)
      rho_ll = u_prim[1, i, j, k, element]
      v1_ll  = u_prim[2, i, j, k, element]
      v2_ll  = u_prim[3, i, j, k, element]
      v3_ll  = u_prim[4, i, j, k, element]
      p_ll   = u_prim[5, i, j, k, element]

      rho_rr = u_prim[1, i, j, kk, element]
      v1_rr  = u_prim[2, i, j, kk, element]
      v2_rr  = u_prim[3, i, j, kk, element]
      v3_rr  = u_prim[4, i, j, kk, element]
      p_rr   = u_prim[5, i, j, kk, element]

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv3_avg = 0.5 * (p_ll * v3_rr + p_rr * v3_ll)

      f1 = rho_avg * v3_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg + p_avg
      f5 = p_avg*v3_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv3_avg

      factor = alpha * derivative_split[k, kk]
      du[1, i, j, k, element] += factor * f1
      du[2, i, j, k, element] += factor * f2
      du[3, i, j, k, element] += factor * f3
      du[4, i, j, k, element] += factor * f4
      du[5, i, j, k, element] += factor * f5

      factor = alpha * derivative_split[kk, k]
      du[1, i, j, kk, element] += factor * f1
      du[2, i, j, kk, element] += factor * f2
      du[3, i, j, kk, element] += factor * f3
      du[4, i, j, kk, element] += factor * f4
      du[5, i, j, kk, element] += factor * f5
    end
  end
end


# pre-compute primitive variables and uses triangular loops without turbo
# after permuting the memory layout
@inline function flux_differencing_kernel_permuted!(_du, u_cons,
                                                    element, mesh::TreeMesh{3},
                                                    nonconservative_terms::Val{false},
                                                    equations::CompressibleEulerEquations3D,
                                                    volume_flux::typeof(flux_shima_etal),
                                                    dg::DGSEM, cache,
                                                    alpha=true)
  # true * [some floating point value] == [exactly the same floating point value]
  # This can (hopefully) be optimized away due to constant propagation.
  @unpack derivative_split, derivative_split_transpose = dg.basis

  # u_prim = PtrArray(pointer(_u_prim),
  #                   (ntuple(_ -> StaticInt(nnodes(dg)), ndims(mesh))...,
  #                   StaticInt(nvariables(equations)),
  #                   nelements(dg, cache)))
  u_prim = StrideArray{eltype(u_cons)}(undef,
    (ntuple(_ -> StaticInt(nnodes(dg)), ndims(mesh))...,
    StaticInt(nvariables(equations))))

  du = StrideArray{eltype(u_cons)}(undef,
    (ntuple(_ -> StaticInt(nnodes(dg)), ndims(mesh))...,
    StaticInt(nvariables(equations))))

  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    rho    = u_cons[1, i, j, k, element]
    rho_v1 = u_cons[2, i, j, k, element]
    rho_v2 = u_cons[3, i, j, k, element]
    rho_v3 = u_cons[4, i, j, k, element]
    rho_e  = u_cons[5, i, j, k, element]

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2 + rho_v3 * v3))

    u_prim[i, j, k, 1] = rho
    u_prim[i, j, k, 2] = v1
    u_prim[i, j, k, 3] = v2
    u_prim[i, j, k, 4] = v3
    u_prim[i, j, k, 5] = p
  end

  # x direction
  # fill!(du, zero(eltype(du)))
  # for i in eachnode(dg), ii in (i+1):nnodes(dg)
  #   @turbo for k in eachnode(dg), j in eachnode(dg)
  #     rho_ll = u_prim[i, j, k, 1]
  #     v1_ll  = u_prim[i, j, k, 2]
  #     v2_ll  = u_prim[i, j, k, 3]
  #     v3_ll  = u_prim[i, j, k, 4]
  #     p_ll   = u_prim[i, j, k, 5]

  #     rho_rr = u_prim[ii, j, k, 1]
  #     v1_rr  = u_prim[ii, j, k, 2]
  #     v2_rr  = u_prim[ii, j, k, 3]
  #     v3_rr  = u_prim[ii, j, k, 4]
  #     p_rr   = u_prim[ii, j, k, 5]

  #     # Average each factor of products in flux
  #     rho_avg = 0.5 * (rho_ll + rho_rr)
  #     v1_avg  = 0.5 * ( v1_ll +  v1_rr)
  #     v2_avg  = 0.5 * ( v2_ll +  v2_rr)
  #     v3_avg  = 0.5 * ( v3_ll +  v3_rr)
  #     p_avg   = 0.5 * (  p_ll +   p_rr)
  #     kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
  #     pv1_avg = 0.5 * (p_ll * v1_rr + p_rr * v1_ll)

  #     f1 = rho_avg * v1_avg
  #     f2 = f1 * v1_avg + p_avg
  #     f3 = f1 * v2_avg
  #     f4 = f1 * v3_avg
  #     f5 = p_avg * v1_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv1_avg

  #     factor = alpha * derivative_split[i, ii]
  #     du[i, j, k, 1] += factor * f1
  #     du[i, j, k, 2] += factor * f2
  #     du[i, j, k, 3] += factor * f3
  #     du[i, j, k, 4] += factor * f4
  #     du[i, j, k, 5] += factor * f5

  #     factor = alpha * derivative_split[ii, i]
  #     du[ii, j, k, 1] += factor * f1
  #     du[ii, j, k, 2] += factor * f2
  #     du[ii, j, k, 3] += factor * f3
  #     du[ii, j, k, 4] += factor * f4
  #     du[ii, j, k, 5] += factor * f5
  #   end
  # end


  du_permuted = StrideArray{eltype(u_cons)}(undef,
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
    StaticInt(nvariables(equations))))

  u_prim_permuted = StrideArray{eltype(u_cons)}(undef,
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
    StaticInt(nvariables(equations))))

  @turbo for v in eachvariable(equations), k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    jk = j + nnodes(dg) * (k- 1)
    u_prim_permuted[jk, i, v] = u_prim[i, j, k, v]
  end
  fill!(du_permuted, zero(eltype(du_permuted)))

  for i in eachnode(dg), ii in (i+1):nnodes(dg)
    @turbo for jk in Base.OneTo(nnodes(dg)^2)
      rho_ll = u_prim_permuted[jk, i, 1]
      v1_ll  = u_prim_permuted[jk, i, 2]
      v2_ll  = u_prim_permuted[jk, i, 3]
      v3_ll  = u_prim_permuted[jk, i, 4]
      p_ll   = u_prim_permuted[jk, i, 5]

      rho_rr = u_prim_permuted[jk, ii, 1]
      v1_rr  = u_prim_permuted[jk, ii, 2]
      v2_rr  = u_prim_permuted[jk, ii, 3]
      v3_rr  = u_prim_permuted[jk, ii, 4]
      p_rr   = u_prim_permuted[jk, ii, 5]

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv1_avg = 0.5 * (p_ll * v1_rr + p_rr * v1_ll)

      f1 = rho_avg * v1_avg
      f2 = f1 * v1_avg + p_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg
      f5 = p_avg * v1_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv1_avg

      factor = alpha * derivative_split[i, ii]
      du_permuted[jk, i, 1] += factor * f1
      du_permuted[jk, i, 2] += factor * f2
      du_permuted[jk, i, 3] += factor * f3
      du_permuted[jk, i, 4] += factor * f4
      du_permuted[jk, i, 5] += factor * f5

      factor = alpha * derivative_split[ii, i]
      du_permuted[jk, ii, 1] += factor * f1
      du_permuted[jk, ii, 2] += factor * f2
      du_permuted[jk, ii, 3] += factor * f3
      du_permuted[jk, ii, 4] += factor * f4
      du_permuted[jk, ii, 5] += factor * f5
    end
  end

  @turbo for v in eachvariable(equations), k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    jk = j + nnodes(dg) * (k- 1)
    du[i, j, k, v] = du_permuted[jk, i, v]
  end


  # y direction
  for j in eachnode(dg), jj in (j+1):nnodes(dg)
    @turbo for k in eachnode(dg), i in eachnode(dg)
      rho_ll = u_prim[i, j, k, 1]
      v1_ll  = u_prim[i, j, k, 2]
      v2_ll  = u_prim[i, j, k, 3]
      v3_ll  = u_prim[i, j, k, 4]
      p_ll   = u_prim[i, j, k, 5]

      rho_rr = u_prim[i, jj, k, 1]
      v1_rr  = u_prim[i, jj, k, 2]
      v2_rr  = u_prim[i, jj, k, 3]
      v3_rr  = u_prim[i, jj, k, 4]
      p_rr   = u_prim[i, jj, k, 5]

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv2_avg = 0.5 * (p_ll * v2_rr + p_rr * v2_ll)

      f1 = rho_avg * v2_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg + p_avg
      f4 = f1 * v3_avg
      f5 = p_avg*v2_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv2_avg

      factor = alpha * derivative_split[j, jj]
      du[i, j, k, 1] += factor * f1
      du[i, j, k, 2] += factor * f2
      du[i, j, k, 3] += factor * f3
      du[i, j, k, 4] += factor * f4
      du[i, j, k, 5] += factor * f5

      factor = alpha * derivative_split[jj, j]
      du[i, jj, k, 1] += factor * f1
      du[i, jj, k, 2] += factor * f2
      du[i, jj, k, 3] += factor * f3
      du[i, jj, k, 4] += factor * f4
      du[i, jj, k, 5] += factor * f5
    end
  end


  # permuting appears to be more expensive than justified by possible gains
  # @turbo for v in eachvariable(equations), k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  #   ik = i + nnodes(dg) * (k- 1)
  #   u_prim_permuted[ik, j, v] = u_prim[i, j, k, v]
  # end
  # fill!(du_permuted, zero(eltype(du_permuted)))

  # for j in eachnode(dg), jj in (j+1):nnodes(dg)
  #   @turbo for ik in Base.OneTo(nnodes(dg)^2)
  #     rho_ll = u_prim_permuted[ik, j, 1]
  #     v1_ll  = u_prim_permuted[ik, j, 2]
  #     v2_ll  = u_prim_permuted[ik, j, 3]
  #     v3_ll  = u_prim_permuted[ik, j, 4]
  #     p_ll   = u_prim_permuted[ik, j, 5]

  #     rho_rr = u_prim_permuted[ik, jj, 1]
  #     v1_rr  = u_prim_permuted[ik, jj, 2]
  #     v2_rr  = u_prim_permuted[ik, jj, 3]
  #     v3_rr  = u_prim_permuted[ik, jj, 4]
  #     p_rr   = u_prim_permuted[ik, jj, 5]

  #     # Average each factor of products in flux
  #     rho_avg = 0.5 * (rho_ll + rho_rr)
  #     v1_avg  = 0.5 * ( v1_ll +  v1_rr)
  #     v2_avg  = 0.5 * ( v2_ll +  v2_rr)
  #     v3_avg  = 0.5 * ( v3_ll +  v3_rr)
  #     p_avg   = 0.5 * (  p_ll +   p_rr)
  #     kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
  #     pv2_avg = 0.5 * (p_ll * v2_rr + p_rr * v2_ll)

  #     f1 = rho_avg * v2_avg
  #     f2 = f1 * v1_avg
  #     f3 = f1 * v2_avg + p_avg
  #     f4 = f1 * v3_avg
  #     f5 = p_avg*v2_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv2_avg

  #     factor = alpha * derivative_split[j, jj]
  #     du_permuted[ik, j, 1] += factor * f1
  #     du_permuted[ik, j, 2] += factor * f2
  #     du_permuted[ik, j, 3] += factor * f3
  #     du_permuted[ik, j, 4] += factor * f4
  #     du_permuted[ik, j, 5] += factor * f5

  #     factor = alpha * derivative_split[jj, j]
  #     du_permuted[ik, jj, 1] += factor * f1
  #     du_permuted[ik, jj, 2] += factor * f2
  #     du_permuted[ik, jj, 3] += factor * f3
  #     du_permuted[ik, jj, 4] += factor * f4
  #     du_permuted[ik, jj, 5] += factor * f5
  #   end
  # end

  # @turbo for v in eachvariable(equations), k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  #   ik = i + nnodes(dg) * (k- 1)
  #   du[i, j, k, v] += du_permuted[ik, j, v]
  # end


  # z direction
  # for k in eachnode(dg), kk in (k+1):nnodes(dg)
  #   @turbo for j in eachnode(dg), i in eachnode(dg)
  #     rho_ll = u_prim[i, j, k, 1]
  #     v1_ll  = u_prim[i, j, k, 2]
  #     v2_ll  = u_prim[i, j, k, 3]
  #     v3_ll  = u_prim[i, j, k, 4]
  #     p_ll   = u_prim[i, j, k, 5]

  #     rho_rr = u_prim[i, j, kk, 1]
  #     v1_rr  = u_prim[i, j, kk, 2]
  #     v2_rr  = u_prim[i, j, kk, 3]
  #     v3_rr  = u_prim[i, j, kk, 4]
  #     p_rr   = u_prim[i, j, kk, 5]

  #     # Average each factor of products in flux
  #     rho_avg = 0.5 * (rho_ll + rho_rr)
  #     v1_avg  = 0.5 * ( v1_ll +  v1_rr)
  #     v2_avg  = 0.5 * ( v2_ll +  v2_rr)
  #     v3_avg  = 0.5 * ( v3_ll +  v3_rr)
  #     p_avg   = 0.5 * (  p_ll +   p_rr)
  #     kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
  #     pv3_avg = 0.5 * (p_ll * v3_rr + p_rr * v3_ll)

  #     f1 = rho_avg * v3_avg
  #     f2 = f1 * v1_avg
  #     f3 = f1 * v2_avg
  #     f4 = f1 * v3_avg + p_avg
  #     f5 = p_avg*v3_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv3_avg

  #     factor = alpha * derivative_split[k, kk]
  #     du[i, j, k, 1] += factor * f1
  #     du[i, j, k, 2] += factor * f2
  #     du[i, j, k, 3] += factor * f3
  #     du[i, j, k, 4] += factor * f4
  #     du[i, j, k, 5] += factor * f5

  #     factor = alpha * derivative_split[kk, k]
  #     du[i, j, kk, 1] += factor * f1
  #     du[i, j, kk, 2] += factor * f2
  #     du[i, j, kk, 3] += factor * f3
  #     du[i, j, kk, 4] += factor * f4
  #     du[i, j, kk, 5] += factor * f5
  #   end
  # end

  u_prim_perm = PtrArray(pointer(u_prim),
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
     StaticInt(nvariables(equations))))
  du_perm = PtrArray(pointer(du),
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
     StaticInt(nvariables(equations))))
  for k in eachnode(dg), kk in (k+1):nnodes(dg)
    @turbo for ij in Base.OneTo(nnodes(dg)^2)
      rho_ll = u_prim_perm[ij, k, 1]
      v1_ll  = u_prim_perm[ij, k, 2]
      v2_ll  = u_prim_perm[ij, k, 3]
      v3_ll  = u_prim_perm[ij, k, 4]
      p_ll   = u_prim_perm[ij, k, 5]

      rho_rr = u_prim_perm[ij, kk, 1]
      v1_rr  = u_prim_perm[ij, kk, 2]
      v2_rr  = u_prim_perm[ij, kk, 3]
      v3_rr  = u_prim_perm[ij, kk, 4]
      p_rr   = u_prim_perm[ij, kk, 5]

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv3_avg = 0.5 * (p_ll * v3_rr + p_rr * v3_ll)

      f1 = rho_avg * v3_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg + p_avg
      f5 = p_avg*v3_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv3_avg

      factor = alpha * derivative_split[k, kk]
      du_perm[ij, k, 1] += factor * f1
      du_perm[ij, k, 2] += factor * f2
      du_perm[ij, k, 3] += factor * f3
      du_perm[ij, k, 4] += factor * f4
      du_perm[ij, k, 5] += factor * f5

      factor = alpha * derivative_split[kk, k]
      du_perm[ij, kk, 1] += factor * f1
      du_perm[ij, kk, 2] += factor * f2
      du_perm[ij, kk, 3] += factor * f3
      du_perm[ij, kk, 4] += factor * f4
      du_perm[ij, kk, 5] += factor * f5
    end
  end

  @turbo for v in eachvariable(equations), k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  # @turbo for v in indices((_du, du), (1, 4)),
  #            k in indices((_du, du), (4, 3)),
  #            j in indices((_du, du), (3, 2)),
  #            i in indices((_du, du), (2, 1))
    _du[v, i, j, k, element] += du[i, j, k, v]
  end
end


@inline function flux_differencing_kernel!(du::StructArray, u::StructArray,
                                           element, mesh::TreeMesh{3},
                                           nonconservative_terms::Val{false},
                                           equations::CompressibleEulerEquations3D,
                                           volume_flux::typeof(flux_shima_etal),
                                           dg::DGSEM, cache,
                                           alpha=true)
  # true * [some floating point value] == [exactly the same floating point value]
  # This can (hopefully) be optimized away due to constant propagation.
  @unpack derivative_split = dg.basis

  u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e = StructArrays.components(u)
  du_rho, du_rho_v1, du_rho_v2, du_rho_v3, du_rho_e = StructArrays.components(du)

  # x direction
  # @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  @turbo for k in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e,
                           du_rho, du_rho_v1, du_rho_v2, du_rho_v3, du_rho_e), 3),
             j in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e,
                           du_rho, du_rho_v1, du_rho_v2, du_rho_v3, du_rho_e), 2),
             i in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e,
                           du_rho, du_rho_v1, du_rho_v2, du_rho_v3, du_rho_e), 1)
    rho_ll    = u_rho[   i, j, k, element]
    rho_v1_ll = u_rho_v1[i, j, k, element]
    rho_v2_ll = u_rho_v2[i, j, k, element]
    rho_v3_ll = u_rho_v3[i, j, k, element]
    rho_e_ll  = u_rho_e[ i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du_rho))
    du_2 = zero(eltype(du_rho_v1))
    du_3 = zero(eltype(du_rho_v2))
    du_4 = zero(eltype(du_rho_v3))
    du_5 = zero(eltype(du_rho_e))

    # for ii in eachnode(dg)
    for ii in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e), 1)
      rho_rr    = u_rho[   ii, j, k, element]
      rho_v1_rr = u_rho_v1[ii, j, k, element]
      rho_v2_rr = u_rho_v2[ii, j, k, element]
      rho_v3_rr = u_rho_v3[ii, j, k, element]
      rho_e_rr  = u_rho_e[ ii, j, k, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv1_avg = 0.5 * (p_ll * v1_rr + p_rr * v1_ll)

      f1 = rho_avg * v1_avg
      f2 = f1 * v1_avg + p_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg
      f5 = p_avg * v1_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv1_avg

      factor = alpha * derivative_split[i, ii]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du_rho[   i, j, k, element] += du_1
    du_rho_v1[i, j, k, element] += du_2
    du_rho_v2[i, j, k, element] += du_3
    du_rho_v3[i, j, k, element] += du_4
    du_rho_e[ i, j, k, element] += du_5
  end

  # y direction
  # @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  @turbo for k in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e,
                           du_rho, du_rho_v1, du_rho_v2, du_rho_v3, du_rho_e), 3),
             j in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e,
                           du_rho, du_rho_v1, du_rho_v2, du_rho_v3, du_rho_e), 2),
             i in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e,
                           du_rho, du_rho_v1, du_rho_v2, du_rho_v3, du_rho_e), 1)
    rho_ll    = u_rho[   i, j, k, element]
    rho_v1_ll = u_rho_v1[i, j, k, element]
    rho_v2_ll = u_rho_v2[i, j, k, element]
    rho_v3_ll = u_rho_v3[i, j, k, element]
    rho_e_ll  = u_rho_e[ i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du_rho))
    du_2 = zero(eltype(du_rho_v1))
    du_3 = zero(eltype(du_rho_v2))
    du_4 = zero(eltype(du_rho_v3))
    du_5 = zero(eltype(du_rho_e))

    # for jj in eachnode(dg)
    for jj in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e), 2)
      rho_rr    = u_rho[   i, jj, k, element]
      rho_v1_rr = u_rho_v1[i, jj, k, element]
      rho_v2_rr = u_rho_v2[i, jj, k, element]
      rho_v3_rr = u_rho_v3[i, jj, k, element]
      rho_e_rr  = u_rho_e[ i, jj, k, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv2_avg = 0.5 * (p_ll * v2_rr + p_rr * v2_ll)

      f1 = rho_avg * v2_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg + p_avg
      f4 = f1 * v3_avg
      f5 = p_avg*v2_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv2_avg

      factor = alpha * derivative_split[j, jj]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du_rho[   i, j, k, element] += du_1
    du_rho_v1[i, j, k, element] += du_2
    du_rho_v2[i, j, k, element] += du_3
    du_rho_v3[i, j, k, element] += du_4
    du_rho_e[ i, j, k, element] += du_5
  end

  # z direction
  # @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  @turbo for k in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e,
                           du_rho, du_rho_v1, du_rho_v2, du_rho_v3, du_rho_e), 3),
             j in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e,
                           du_rho, du_rho_v1, du_rho_v2, du_rho_v3, du_rho_e), 2),
             i in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e,
                           du_rho, du_rho_v1, du_rho_v2, du_rho_v3, du_rho_e), 1)
    rho_ll    = u_rho[   i, j, k, element]
    rho_v1_ll = u_rho_v1[i, j, k, element]
    rho_v2_ll = u_rho_v2[i, j, k, element]
    rho_v3_ll = u_rho_v3[i, j, k, element]
    rho_e_ll  = u_rho_e[ i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du_rho))
    du_2 = zero(eltype(du_rho_v1))
    du_3 = zero(eltype(du_rho_v2))
    du_4 = zero(eltype(du_rho_v3))
    du_5 = zero(eltype(du_rho_e))

    # for kk in eachnode(dg)
    for kk in indices((u_rho, u_rho_v1, u_rho_v2, u_rho_v3, u_rho_e), 3)
      rho_rr    = u_rho[   i, j, kk, element]
      rho_v1_rr = u_rho_v1[i, j, kk, element]
      rho_v2_rr = u_rho_v2[i, j, kk, element]
      rho_v3_rr = u_rho_v3[i, j, kk, element]
      rho_e_rr  = u_rho_e[ i, j, kk, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      kin_avg = 0.5 * (v1_ll * v1_rr + v2_ll * v2_rr + v3_ll * v3_rr)
      pv3_avg = 0.5 * (p_ll * v3_rr + p_rr * v3_ll)

      f1 = rho_avg * v3_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg + p_avg
      f5 = p_avg*v3_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv3_avg

      factor = alpha * derivative_split[k, kk]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du_rho[   i, j, k, element] += du_1
    du_rho_v1[i, j, k, element] += du_2
    du_rho_v2[i, j, k, element] += du_3
    du_rho_v3[i, j, k, element] += du_4
    du_rho_e[ i, j, k, element] += du_5
  end
end






################################################################################
# flux_ranocha
################################################################################


# plain version
@inline function flux_differencing_kernel!(du, u,
                                           element, mesh::TreeMesh{3},
                                           nonconservative_terms::Val{false},
                                           equations::CompressibleEulerEquations3D,
                                           volume_flux::typeof(flux_ranocha),
                                           dg::DGSEM, cache,
                                           alpha=true)
  # true * [some floating point value] == [exactly the same floating point value]
  # This can (hopefully) be optimized away due to constant propagation.
  @unpack derivative_split, derivative_split_transpose = dg.basis

  # x direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    rho_ll    = u[1, i, j, k, element]
    rho_v1_ll = u[2, i, j, k, element]
    rho_v2_ll = u[3, i, j, k, element]
    rho_v3_ll = u[4, i, j, k, element]
    rho_e_ll  = u[5, i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for ii in eachnode(dg)
      rho_rr    = u[1, ii, j, k, element]
      rho_v1_rr = u[2, ii, j, k, element]
      rho_v2_rr = u[3, ii, j, k, element]
      rho_v3_rr = u[4, ii, j, k, element]
      rho_e_rr  = u[5, ii, j, k, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Compute the necessary mean values
      rho_mean = ln_mean(rho_ll, rho_rr)
      # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
      # in exact arithmetic since
      #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
      #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
      inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
      v1_avg = 0.5 * (v1_ll + v1_rr)
      v2_avg = 0.5 * (v2_ll + v2_rr)
      v3_avg = 0.5 * (v3_ll + v3_rr)
      p_avg  = 0.5 * (p_ll + p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      f1 = rho_mean * v1_avg
      f2 = f1 * v1_avg + p_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg
      f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v1_rr + p_rr*v1_ll)

      factor = alpha * derivative_split[i, ii]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end

  # y direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
  # @turbo for k in indices((u, du), 4), j in indices((u, du, derivative_split), (3, 3, 1)), i in indices((u, du), 2)
    rho_ll    = u[1, i, j, k, element]
    rho_v1_ll = u[2, i, j, k, element]
    rho_v2_ll = u[3, i, j, k, element]
    rho_v3_ll = u[4, i, j, k, element]
    rho_e_ll  = u[5, i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for jj in eachnode(dg)
      rho_rr    = u[1, i, jj, k, element]
      rho_v1_rr = u[2, i, jj, k, element]
      rho_v2_rr = u[3, i, jj, k, element]
      rho_v3_rr = u[4, i, jj, k, element]
      rho_e_rr  = u[5, i, jj, k, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Compute the necessary mean values
      rho_mean = ln_mean(rho_ll, rho_rr)
      # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
      # in exact arithmetic since
      #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
      #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
      inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
      v1_avg = 0.5 * (v1_ll + v1_rr)
      v2_avg = 0.5 * (v2_ll + v2_rr)
      v3_avg = 0.5 * (v3_ll + v3_rr)
      p_avg  = 0.5 * (p_ll + p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      f1 = rho_mean * v2_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg + p_avg
      f4 = f1 * v3_avg
      f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v2_rr + p_rr*v2_ll)

      factor = alpha * derivative_split[j, jj]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end

  # z direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
  # @turbo for k in indices((u, du, derivative_split), (4, 4, 1)), j in indices((u, du), 3), i in indices((u, du), 2)
    rho_ll    = u[1, i, j, k, element]
    rho_v1_ll = u[2, i, j, k, element]
    rho_v2_ll = u[3, i, j, k, element]
    rho_v3_ll = u[4, i, j, k, element]
    rho_e_ll  = u[5, i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for kk in eachnode(dg)
      rho_rr    = u[1, i, j, kk, element]
      rho_v1_rr = u[2, i, j, kk, element]
      rho_v2_rr = u[3, i, j, kk, element]
      rho_v3_rr = u[4, i, j, kk, element]
      rho_e_rr  = u[5, i, j, kk, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Compute the necessary mean values
      rho_mean = ln_mean(rho_ll, rho_rr)
      # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
      # in exact arithmetic since
      #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
      #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
      inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
      v1_avg = 0.5 * (v1_ll + v1_rr)
      v2_avg = 0.5 * (v2_ll + v2_rr)
      v3_avg = 0.5 * (v3_ll + v3_rr)
      p_avg  = 0.5 * (p_ll + p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      f1 = rho_mean * v3_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg + p_avg
      f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v3_rr + p_rr*v3_ll)

      factor = alpha * derivative_split[k, kk]
      # factor = alpha * derivative_split_transpose[kk, k]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end
end

# make LV aware of logmean, doesn't seem to help much
@inline function flux_differencing_kernel_inline_logmean!(du, u,
                                           element, mesh::TreeMesh{3},
                                           nonconservative_terms::Val{false},
                                           equations::CompressibleEulerEquations3D,
                                           volume_flux::typeof(flux_ranocha),
                                           dg::DGSEM, cache,
                                           alpha=true)
  # true * [some floating point value] == [exactly the same floating point value]
  # This can (hopefully) be optimized away due to constant propagation.
  @unpack derivative_split, derivative_split_transpose = dg.basis

  # x direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    rho_ll    = u[1, i, j, k, element]
    rho_v1_ll = u[2, i, j, k, element]
    rho_v2_ll = u[3, i, j, k, element]
    rho_v3_ll = u[4, i, j, k, element]
    rho_e_ll  = u[5, i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for ii in eachnode(dg)
      rho_rr    = u[1, ii, j, k, element]
      rho_v1_rr = u[2, ii, j, k, element]
      rho_v2_rr = u[3, ii, j, k, element]
      rho_v3_rr = u[4, ii, j, k, element]
      rho_e_rr  = u[5, ii, j, k, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Compute the necessary mean values
      # rho_mean = ln_mean(rho_ll, rho_rr)
      x = rho_ll
      y = rho_rr
      f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
      # (x + y) / @evalpoly(f2, 2, 2/3, 2/5, 2/7)
      special_path1 = (x + y) / (2 + f2*(0.666666666666667 + f2*(0.4 + 0.285714285714286*f2)))
      regular_path1 = (y - x) / log(y / x)
      rho_mean = ifelse(f2 < 1.0e-4, special_path1, regular_path1)

      # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
      # in exact arithmetic since
      #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
      #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
      inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
      # using the explicit form below doesn't seem to help
      # x2 = rho_ll * p_rr
      # y2 = rho_rr * p_ll
      # f22 = (x2 * (x2 - 2 * y2) + y2 * y2) / (x2 * (x2 + 2 * y2) + y2 * y2) # f2 = f^2
      # # @evalpoly(f2, 2, 2/3, 2/5, 2/7) / (x + y)
      # special_path2 = (2 + f2*(0.666666666666667 + f2*(0.4 + 0.285714285714286*f2))) / (x2 + y2)
      # regular_path2 = log(y2 / x2) / (y2 - x2)
      # inv_rho_p_mean = p_ll * p_rr * ifelse(f22 < 1.0e-4, special_path2, regular_path2)

      v1_avg = 0.5 * (v1_ll + v1_rr)
      v2_avg = 0.5 * (v2_ll + v2_rr)
      v3_avg = 0.5 * (v3_ll + v3_rr)
      p_avg  = 0.5 * (p_ll + p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      f1 = rho_mean * v1_avg
      f2 = f1 * v1_avg + p_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg
      f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v1_rr + p_rr*v1_ll)

      factor = alpha * derivative_split[i, ii]
      # factor = alpha * derivative_split_transpose[ii, i]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end

  # y direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
  # @turbo for k in indices((u, du), 4), j in indices((u, du, derivative_split), (3, 3, 1)), i in indices((u, du), 2)
    rho_ll    = u[1, i, j, k, element]
    rho_v1_ll = u[2, i, j, k, element]
    rho_v2_ll = u[3, i, j, k, element]
    rho_v3_ll = u[4, i, j, k, element]
    rho_e_ll  = u[5, i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for jj in eachnode(dg)
    # for jj in indices(u, 3)
    # for jj in indices((u, derivative_split), (3, 2))
      rho_rr    = u[1, i, jj, k, element]
      rho_v1_rr = u[2, i, jj, k, element]
      rho_v2_rr = u[3, i, jj, k, element]
      rho_v3_rr = u[4, i, jj, k, element]
      rho_e_rr  = u[5, i, jj, k, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Compute the necessary mean values
      rho_mean = ln_mean(rho_ll, rho_rr)
      # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
      # in exact arithmetic since
      #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
      #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
      inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
      v1_avg = 0.5 * (v1_ll + v1_rr)
      v2_avg = 0.5 * (v2_ll + v2_rr)
      v3_avg = 0.5 * (v3_ll + v3_rr)
      p_avg  = 0.5 * (p_ll + p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      f1 = rho_mean * v2_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg + p_avg
      f4 = f1 * v3_avg
      f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v2_rr + p_rr*v2_ll)

      factor = alpha * derivative_split[j, jj]
      # factor = alpha * derivative_split_transpose[jj, j]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end

  # z direction
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  # @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
  # @turbo for k in indices((u, du, derivative_split), (4, 4, 1)), j in indices((u, du), 3), i in indices((u, du), 2)
    rho_ll    = u[1, i, j, k, element]
    rho_v1_ll = u[2, i, j, k, element]
    rho_v2_ll = u[3, i, j, k, element]
    rho_v3_ll = u[4, i, j, k, element]
    rho_e_ll  = u[5, i, j, k, element]

    v1_ll = rho_v1_ll / rho_ll
    v2_ll = rho_v2_ll / rho_ll
    v3_ll = rho_v3_ll / rho_ll
    p_ll  = (equations.gamma - 1) * (
      rho_e_ll - 0.5 * (rho_v1_ll * v1_ll + rho_v2_ll * v2_ll + rho_v3_ll * v3_ll))

    du_1 = zero(eltype(du))
    du_2 = zero(eltype(du))
    du_3 = zero(eltype(du))
    du_4 = zero(eltype(du))
    du_5 = zero(eltype(du))

    for kk in eachnode(dg)
    # for kk in indices(u, 4)
    # for kk in indices((u, derivative_split), (4, 2))
      rho_rr    = u[1, i, j, kk, element]
      rho_v1_rr = u[2, i, j, kk, element]
      rho_v2_rr = u[3, i, j, kk, element]
      rho_v3_rr = u[4, i, j, kk, element]
      rho_e_rr  = u[5, i, j, kk, element]

      v1_rr = rho_v1_rr / rho_rr
      v2_rr = rho_v2_rr / rho_rr
      v3_rr = rho_v3_rr / rho_rr
      p_rr  = (equations.gamma - 1) * (
        rho_e_rr - 0.5 * (rho_v1_rr * v1_rr + rho_v2_rr * v2_rr + rho_v3_rr * v3_rr))

      # Compute the necessary mean values
      rho_mean = ln_mean(rho_ll, rho_rr)
      # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
      # in exact arithmetic since
      #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
      #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
      inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
      v1_avg = 0.5 * (v1_ll + v1_rr)
      v2_avg = 0.5 * (v2_ll + v2_rr)
      v3_avg = 0.5 * (v3_ll + v3_rr)
      p_avg  = 0.5 * (p_ll + p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      f1 = rho_mean * v3_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg + p_avg
      f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v3_rr + p_rr*v3_ll)

      factor = alpha * derivative_split[k, kk]
      # factor = alpha * derivative_split_transpose[kk, k]
      du_1 += factor * f1
      du_2 += factor * f2
      du_3 += factor * f3
      du_4 += factor * f4
      du_5 += factor * f5
    end

    du[1, i, j, k, element] += du_1
    du[2, i, j, k, element] += du_2
    du[3, i, j, k, element] += du_3
    du[4, i, j, k, element] += du_4
    du[5, i, j, k, element] += du_5
  end
end


# pre-compute primitive variables and uses triangular loops without turbo
# after permuting the memory layout
@inline function flux_differencing_kernel_permuted!(_du, u_cons,
                                                    element, mesh::TreeMesh{3},
                                                    nonconservative_terms::Val{false},
                                                    equations::CompressibleEulerEquations3D,
                                                    volume_flux::typeof(flux_ranocha),
                                                    dg::DGSEM, cache,
                                                    alpha=true)
  # true * [some floating point value] == [exactly the same floating point value]
  # This can (hopefully) be optimized away due to constant propagation.
  @unpack derivative_split, derivative_split_transpose = dg.basis


  u_prim = StrideArray{eltype(u_cons)}(undef,
    (ntuple(_ -> StaticInt(nnodes(dg)), ndims(mesh))...,
    StaticInt(nvariables(equations))))

  du = StrideArray{eltype(u_cons)}(undef,
    (ntuple(_ -> StaticInt(nnodes(dg)), ndims(mesh))...,
    StaticInt(nvariables(equations))))
  # fill!(du, zero(eltype(du))) # TODO: Can be removed

  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    rho    = u_cons[1, i, j, k, element]
    rho_v1 = u_cons[2, i, j, k, element]
    rho_v2 = u_cons[3, i, j, k, element]
    rho_v3 = u_cons[4, i, j, k, element]
    rho_e  = u_cons[5, i, j, k, element]

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2 + rho_v3 * v3))

    u_prim[i, j, k, 1] = rho
    u_prim[i, j, k, 2] = v1
    u_prim[i, j, k, 3] = v2
    u_prim[i, j, k, 4] = v3
    u_prim[i, j, k, 5] = p
  end


  # x direction
  du_permuted = StrideArray{eltype(u_cons)}(undef,
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
    StaticInt(nvariables(equations))))

  u_prim_permuted = StrideArray{eltype(u_cons)}(undef,
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
    StaticInt(nvariables(equations))))

  @turbo for v in eachvariable(equations), k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    jk = j + nnodes(dg) * (k- 1)
    u_prim_permuted[jk, i, v] = u_prim[i, j, k, v]
  end
  fill!(du_permuted, zero(eltype(du_permuted)))

  for i in eachnode(dg), ii in (i+1):nnodes(dg)
    @turbo for jk in Base.OneTo(nnodes(dg)^2)
    # @turbo for jk in indices((du_permuted, u_prim_permuted), 1)
      rho_ll = u_prim_permuted[jk, i, 1]
      v1_ll  = u_prim_permuted[jk, i, 2]
      v2_ll  = u_prim_permuted[jk, i, 3]
      v3_ll  = u_prim_permuted[jk, i, 4]
      p_ll   = u_prim_permuted[jk, i, 5]

      rho_rr = u_prim_permuted[jk, ii, 1]
      v1_rr  = u_prim_permuted[jk, ii, 2]
      v2_rr  = u_prim_permuted[jk, ii, 3]
      v3_rr  = u_prim_permuted[jk, ii, 4]
      p_rr   = u_prim_permuted[jk, ii, 5]

      # Compute the necessary mean values
      rho_mean = ln_mean(rho_ll, rho_rr)
      # x1 = rho_ll
      # y1 = rho_rr
      # rho_mean = (y1 - x1) / log(y1 / x1) # TODO: logmean
      # z1 = (x1 * (x1 - 2 * y1) + y1 * y1) / (x1 * (x1 + 2 * y1) + y1 * y1)
      # special_path1 = (x1 + y1) / (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1)))
      # regular_path1 = (y1 - x1) / log(y1 / x1)
      # rho_mean = ifelse(z1 < 1.0e-4, special_path1, regular_path1)

      # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
      # in exact arithmetic since
      #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
      #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
      inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
      # x2 = rho_ll * p_rr
      # y2 = rho_rr * p_ll
      # inv_rho_p_mean = p_ll * p_rr * log(y2 / x2) / (y2 - x2) # TODO: logmean
      # z2 = (x2 * (x2 - 2 * y2) + y2 * y2) / (x2 * (x2 + 2 * y2) + y2 * y2)
      # special_path2 = (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1))) / (x2 + y2)
      # regular_path2 = log(y2 / x2) / (y2 - x2)
      # inv_rho_p_mean = p_ll * p_rr * ifelse(z2 < 1.0e-4, special_path2, regular_path2)

      v1_avg = 0.5 * (v1_ll + v1_rr)
      v2_avg = 0.5 * (v2_ll + v2_rr)
      v3_avg = 0.5 * (v3_ll + v3_rr)
      p_avg  = 0.5 * (p_ll + p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      f1 = rho_mean * v1_avg
      f2 = f1 * v1_avg + p_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg
      f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v1_rr + p_rr*v1_ll)

      factor = alpha * derivative_split[i, ii]
      du_permuted[jk, i, 1] += factor * f1
      du_permuted[jk, i, 2] += factor * f2
      du_permuted[jk, i, 3] += factor * f3
      du_permuted[jk, i, 4] += factor * f4
      du_permuted[jk, i, 5] += factor * f5

      factor = alpha * derivative_split[ii, i]
      du_permuted[jk, ii, 1] += factor * f1
      du_permuted[jk, ii, 2] += factor * f2
      du_permuted[jk, ii, 3] += factor * f3
      du_permuted[jk, ii, 4] += factor * f4
      du_permuted[jk, ii, 5] += factor * f5
    end
  end

  @turbo for v in eachvariable(equations), k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    jk = j + nnodes(dg) * (k- 1)
    du[i, j, k, v] = du_permuted[jk, i, v]
  end


  # y direction
  for j in eachnode(dg), jj in (j+1):nnodes(dg)
    @turbo for k in eachnode(dg), i in eachnode(dg)
    # @turbo for k in eachnode(dg), i in eachnode(dg)
      rho_ll = u_prim[i, j, k, 1]
      v1_ll  = u_prim[i, j, k, 2]
      v2_ll  = u_prim[i, j, k, 3]
      v3_ll  = u_prim[i, j, k, 4]
      p_ll   = u_prim[i, j, k, 5]

      rho_rr = u_prim[i, jj, k, 1]
      v1_rr  = u_prim[i, jj, k, 2]
      v2_rr  = u_prim[i, jj, k, 3]
      v3_rr  = u_prim[i, jj, k, 4]
      p_rr   = u_prim[i, jj, k, 5]

      # Compute the necessary mean values
      rho_mean = ln_mean(rho_ll, rho_rr)
      # x1 = rho_ll
      # y1 = rho_rr
      # rho_mean = (y1 - x1) / log(y1 / x1) # TODO: logmean
      # z1 = (x1 * (x1 - 2 * y1) + y1 * y1) / (x1 * (x1 + 2 * y1) + y1 * y1)
      # special_path1 = (x1 + y1) / (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1)))
      # regular_path1 = (y1 - x1) / log(y1 / x1)
      # rho_mean = ifelse(z1 < 1.0e-4, special_path1, regular_path1)

      # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
      # in exact arithmetic since
      #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
      #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
      inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
      # x2 = rho_ll * p_rr
      # y2 = rho_rr * p_ll
      # inv_rho_p_mean = p_ll * p_rr * log(y2 / x2) / (y2 - x2) # TODO: logmean
      # z2 = (x2 * (x2 - 2 * y2) + y2 * y2) / (x2 * (x2 + 2 * y2) + y2 * y2)
      # special_path2 = (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1))) / (x2 + y2)
      # regular_path2 = log(y2 / x2) / (y2 - x2)
      # inv_rho_p_mean = p_ll * p_rr * ifelse(z2 < 1.0e-4, special_path2, regular_path2)

      v1_avg = 0.5 * (v1_ll + v1_rr)
      v2_avg = 0.5 * (v2_ll + v2_rr)
      v3_avg = 0.5 * (v3_ll + v3_rr)
      p_avg  = 0.5 * (p_ll + p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      f1 = rho_mean * v2_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg + p_avg
      f4 = f1 * v3_avg
      f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v2_rr + p_rr*v2_ll)

      factor = alpha * derivative_split[j, jj]
      du[i, j, k, 1] += factor * f1
      du[i, j, k, 2] += factor * f2
      du[i, j, k, 3] += factor * f3
      du[i, j, k, 4] += factor * f4
      du[i, j, k, 5] += factor * f5

      factor = alpha * derivative_split[jj, j]
      du[i, jj, k, 1] += factor * f1
      du[i, jj, k, 2] += factor * f2
      du[i, jj, k, 3] += factor * f3
      du[i, jj, k, 4] += factor * f4
      du[i, jj, k, 5] += factor * f5
    end
  end


  # z direction
  u_prim_perm = PtrArray(pointer(u_prim),
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
     StaticInt(nvariables(equations))))
  du_perm = PtrArray(pointer(du),
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
     StaticInt(nvariables(equations))))
  for k in eachnode(dg), kk in (k+1):nnodes(dg)
    @turbo for ij in Base.OneTo(nnodes(dg)^2)
    # @turbo for ij in Base.OneTo(nnodes(dg)^2)
      rho_ll = u_prim_perm[ij, k, 1]
      v1_ll  = u_prim_perm[ij, k, 2]
      v2_ll  = u_prim_perm[ij, k, 3]
      v3_ll  = u_prim_perm[ij, k, 4]
      p_ll   = u_prim_perm[ij, k, 5]

      rho_rr = u_prim_perm[ij, kk, 1]
      v1_rr  = u_prim_perm[ij, kk, 2]
      v2_rr  = u_prim_perm[ij, kk, 3]
      v3_rr  = u_prim_perm[ij, kk, 4]
      p_rr   = u_prim_perm[ij, kk, 5]

      # Compute the necessary mean values
      rho_mean = ln_mean(rho_ll, rho_rr)
      # x1 = rho_ll
      # y1 = rho_rr
      # rho_mean = (y1 - x1) / log(y1 / x1) # TODO: logmean
      # z1 = (x1 * (x1 - 2 * y1) + y1 * y1) / (x1 * (x1 + 2 * y1) + y1 * y1)
      # special_path1 = (x1 + y1) / (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1)))
      # regular_path1 = (y1 - x1) / log(y1 / x1)
      # rho_mean = ifelse(z1 < 1.0e-4, special_path1, regular_path1)

      # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
      # in exact arithmetic since
      #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
      #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
      inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
      # x2 = rho_ll * p_rr
      # y2 = rho_rr * p_ll
      # inv_rho_p_mean = p_ll * p_rr * log(y2 / x2) / (y2 - x2) # TODO: logmean
      # z2 = (x2 * (x2 - 2 * y2) + y2 * y2) / (x2 * (x2 + 2 * y2) + y2 * y2)
      # special_path2 = (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1))) / (x2 + y2)
      # regular_path2 = log(y2 / x2) / (y2 - x2)
      # inv_rho_p_mean = p_ll * p_rr * ifelse(z2 < 1.0e-4, special_path2, regular_path2)

      v1_avg = 0.5 * (v1_ll + v1_rr)
      v2_avg = 0.5 * (v2_ll + v2_rr)
      v3_avg = 0.5 * (v3_ll + v3_rr)
      p_avg  = 0.5 * (p_ll + p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      f1 = rho_mean * v3_avg
      f2 = f1 * v1_avg
      f3 = f1 * v2_avg
      f4 = f1 * v3_avg + p_avg
      f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v3_rr + p_rr*v3_ll)

      factor = alpha * derivative_split[k, kk]
      du_perm[ij, k, 1] += factor * f1
      du_perm[ij, k, 2] += factor * f2
      du_perm[ij, k, 3] += factor * f3
      du_perm[ij, k, 4] += factor * f4
      du_perm[ij, k, 5] += factor * f5

      factor = alpha * derivative_split[kk, k]
      du_perm[ij, kk, 1] += factor * f1
      du_perm[ij, kk, 2] += factor * f2
      du_perm[ij, kk, 3] += factor * f3
      du_perm[ij, kk, 4] += factor * f4
      du_perm[ij, kk, 5] += factor * f5
    end
  end


  @turbo for v in eachvariable(equations), k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    _du[v, i, j, k, element] += du[i, j, k, v]
  end
end

