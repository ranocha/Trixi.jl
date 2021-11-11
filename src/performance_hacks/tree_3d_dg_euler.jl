
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
# if it needs to specialize for high polynomial degrees at compile time,
# e.g., `polydeg = 7`
# u = Trixi.wrap_array(u_ode, semi); du = Trixi.wrap_array(du_ode, semi);
u = PtrArray(pointer(u_ode), (nvariables(equations), ntuple(_ -> nnodes(solver), ndims(mesh))...,
  nelements(solver, semi.cache)));
du = PtrArray(pointer(du_ode), (nvariables(equations), ntuple(_ -> nnodes(solver), ndims(mesh))...,
  nelements(solver, semi.cache)));
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
  $du_soa, $u_soa, $(Ref(1))[], $(semi.mesh),
  $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
  $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))

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
