
# pre-compute primitive variables and uses triangular loops without turbo
# after permuting the memory layout
# @inline function split_form_kernel!(_du, u_cons,
@inline function flux_differencing_kernel_permuted!(_du, u_cons,
                                                    element, mesh::Union{StructuredMesh{3}, P4estMesh{3}},
                                                    nonconservative_terms::Val{false},
                                                    equations::CompressibleEulerEquations3D,
                                                    volume_flux::typeof(flux_shima_etal),
                                                    dg::DGSEM, cache,
                                                    alpha=true)
  # true * [some floating point value] == [exactly the same floating point value]
  # This can (hopefully) be optimized away due to constant propagation.
  @unpack derivative_split, derivative_split_transpose = dg.basis
  @unpack contravariant_vectors = cache.elements

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

  @turbo for v in eachvariable(equations),
             k in eachnode(dg),
             j in eachnode(dg),
             i in eachnode(dg)
    jk = j + nnodes(dg) * (k- 1)
    u_prim_permuted[jk, i, v] = u_prim[i, j, k, v]
  end
  fill!(du_permuted, zero(eltype(du_permuted)))

  contravariant_vectors_x = StrideArray{eltype(contravariant_vectors)}(undef,
    (StaticInt(nnodes(dg)), StaticInt(nnodes(dg)^2),
    StaticInt(ndims(mesh))))
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    jk = j + nnodes(dg) * (k- 1)
    contravariant_vectors_x[i, jk, 1] = contravariant_vectors[1, 1, i, j, k, element]
    contravariant_vectors_x[i, jk, 2] = contravariant_vectors[2, 1, i, j, k, element]
    contravariant_vectors_x[i, jk, 3] = contravariant_vectors[3, 1, i, j, k, element]
  end

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

      normal_direction_1 = 0.5 * (
        contravariant_vectors_x[i, jk, 1] + contravariant_vectors_x[ii, jk, 1])
      normal_direction_2 = 0.5 * (
        contravariant_vectors_x[i, jk, 2] + contravariant_vectors_x[ii, jk, 2])
      normal_direction_3 = 0.5 * (
        contravariant_vectors_x[i, jk, 3] + contravariant_vectors_x[ii, jk, 3])

      v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 + v3_ll * normal_direction_3
      v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 + v3_rr * normal_direction_3

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      v_dot_n_avg = 0.5 * (v_dot_n_ll + v_dot_n_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      # Calculate fluxes depending on normal_direction
      f1 = rho_avg * v_dot_n_avg
      f2 = f1 * v1_avg + p_avg * normal_direction_1
      f3 = f1 * v2_avg + p_avg * normal_direction_2
      f4 = f1 * v3_avg + p_avg * normal_direction_3
      f5 = ( f1 * velocity_square_avg + p_avg * v_dot_n_avg * equations.inv_gamma_minus_one
            + 0.5 * (p_ll * v_dot_n_rr + p_rr * v_dot_n_ll) )

      # Add scaled fluxes to RHS
      factor_i = alpha * derivative_split[i, ii]
      du_permuted[jk, i, 1] += factor_i * f1
      du_permuted[jk, i, 2] += factor_i * f2
      du_permuted[jk, i, 3] += factor_i * f3
      du_permuted[jk, i, 4] += factor_i * f4
      du_permuted[jk, i, 5] += factor_i * f5

      factor_ii = alpha * derivative_split[ii, i]
      du_permuted[jk, ii, 1] += factor_ii * f1
      du_permuted[jk, ii, 2] += factor_ii * f2
      du_permuted[jk, ii, 3] += factor_ii * f3
      du_permuted[jk, ii, 4] += factor_ii * f4
      du_permuted[jk, ii, 5] += factor_ii * f5
    end
  end

  @turbo for v in eachvariable(equations),
             k in eachnode(dg),
             j in eachnode(dg),
             i in eachnode(dg)
    jk = j + nnodes(dg) * (k- 1)
    du[i, j, k, v] = du_permuted[jk, i, v]
  end


  # y direction
  contravariant_vectors_y = StrideArray{eltype(contravariant_vectors)}(undef,
    (StaticInt(nnodes(dg)), StaticInt(nnodes(dg)), StaticInt(nnodes(dg)),
    StaticInt(ndims(mesh))))
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    contravariant_vectors_y[i, j, k, 1] = contravariant_vectors[1, 2, i, j, k, element]
    contravariant_vectors_y[i, j, k, 2] = contravariant_vectors[2, 2, i, j, k, element]
    contravariant_vectors_y[i, j, k, 3] = contravariant_vectors[3, 2, i, j, k, element]
  end

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

      normal_direction_1 = 0.5 * (
        contravariant_vectors_y[i, j, k, 1] + contravariant_vectors_y[i, jj, k, 1])
      normal_direction_2 = 0.5 * (
        contravariant_vectors_y[i, j, k, 2] + contravariant_vectors_y[i, jj, k, 2])
      normal_direction_3 = 0.5 * (
        contravariant_vectors_y[i, j, k, 3] + contravariant_vectors_y[i, jj, k, 3])

      v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 + v3_ll * normal_direction_3
      v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 + v3_rr * normal_direction_3

      # Average each factor of products in flux
      rho_avg = 0.5 * (rho_ll + rho_rr)
      v1_avg  = 0.5 * ( v1_ll +  v1_rr)
      v2_avg  = 0.5 * ( v2_ll +  v2_rr)
      v3_avg  = 0.5 * ( v3_ll +  v3_rr)
      v_dot_n_avg = 0.5 * (v_dot_n_ll + v_dot_n_rr)
      p_avg   = 0.5 * (  p_ll +   p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      # Calculate fluxes depending on normal_direction
      f1 = rho_avg * v_dot_n_avg
      f2 = f1 * v1_avg + p_avg * normal_direction_1
      f3 = f1 * v2_avg + p_avg * normal_direction_2
      f4 = f1 * v3_avg + p_avg * normal_direction_3
      f5 = ( f1 * velocity_square_avg + p_avg * v_dot_n_avg * equations.inv_gamma_minus_one
            + 0.5 * (p_ll * v_dot_n_rr + p_rr * v_dot_n_ll) )

      # Add scaled fluxes to RHS
      factor_j = alpha * derivative_split[j, jj]
      du[i, j, k, 1] += factor_j * f1
      du[i, j, k, 2] += factor_j * f2
      du[i, j, k, 3] += factor_j * f3
      du[i, j, k, 4] += factor_j * f4
      du[i, j, k, 5] += factor_j * f5

      factor_jj = alpha * derivative_split[jj, j]
      du[i, jj, k, 1] += factor_jj * f1
      du[i, jj, k, 2] += factor_jj * f2
      du[i, jj, k, 3] += factor_jj * f3
      du[i, jj, k, 4] += factor_jj * f4
      du[i, jj, k, 5] += factor_jj * f5
    end
  end


  # z direction
  GC.@preserve u_prim begin
    u_prim_perm = PtrArray(pointer(u_prim),
      (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
      StaticInt(nvariables(equations))))

    du_perm = PtrArray(pointer(du),
      (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
      StaticInt(nvariables(equations))))

    contravariant_vectors_z = StrideArray{eltype(contravariant_vectors)}(undef,
      (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
      StaticInt(ndims(mesh))))
    @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
      ij = i + nnodes(dg) * (j- 1)
      contravariant_vectors_z[ij, k, 1] = contravariant_vectors[1, 3, i, j, k, element]
      contravariant_vectors_z[ij, k, 2] = contravariant_vectors[2, 3, i, j, k, element]
      contravariant_vectors_z[ij, k, 3] = contravariant_vectors[3, 3, i, j, k, element]
    end

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

        normal_direction_1 = 0.5 * (
          contravariant_vectors_z[ij, k, 1] + contravariant_vectors_z[ij, kk, 1])
        normal_direction_2 = 0.5 * (
          contravariant_vectors_z[ij, k, 2] + contravariant_vectors_z[ij, kk, 2])
        normal_direction_3 = 0.5 * (
          contravariant_vectors_z[ij, k, 3] + contravariant_vectors_z[ij, kk, 3])

        v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 + v3_ll * normal_direction_3
        v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 + v3_rr * normal_direction_3

        # Average each factor of products in flux
        rho_avg = 0.5 * (rho_ll + rho_rr)
        v1_avg  = 0.5 * ( v1_ll +  v1_rr)
        v2_avg  = 0.5 * ( v2_ll +  v2_rr)
        v3_avg  = 0.5 * ( v3_ll +  v3_rr)
        v_dot_n_avg = 0.5 * (v_dot_n_ll + v_dot_n_rr)
        p_avg   = 0.5 * (  p_ll +   p_rr)
        velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

        # Calculate fluxes depending on normal_direction
        f1 = rho_avg * v_dot_n_avg
        f2 = f1 * v1_avg + p_avg * normal_direction_1
        f3 = f1 * v2_avg + p_avg * normal_direction_2
        f4 = f1 * v3_avg + p_avg * normal_direction_3
        f5 = ( f1 * velocity_square_avg + p_avg * v_dot_n_avg * equations.inv_gamma_minus_one
              + 0.5 * (p_ll * v_dot_n_rr + p_rr * v_dot_n_ll) )

        # Add scaled fluxes to RHS
        factor_k = alpha * derivative_split[k, kk]
        du_perm[ij, k, 1] += factor_k * f1
        du_perm[ij, k, 2] += factor_k * f2
        du_perm[ij, k, 3] += factor_k * f3
        du_perm[ij, k, 4] += factor_k * f4
        du_perm[ij, k, 5] += factor_k * f5

        factor_kk = alpha * derivative_split[kk, k]
        du_perm[ij, kk, 1] += factor_kk * f1
        du_perm[ij, kk, 2] += factor_kk * f2
        du_perm[ij, kk, 3] += factor_kk * f3
        du_perm[ij, kk, 4] += factor_kk * f4
        du_perm[ij, kk, 5] += factor_kk * f5
      end
    end
  end # GC.@preserve u_prim begin

  @turbo for v in eachvariable(equations),
             k in eachnode(dg),
             j in eachnode(dg),
             i in eachnode(dg)
    _du[v, i, j, k, element] += du[i, j, k, v]
  end
end

