
@inline function flux_differencing_kernel!(du::AbstractArray{<:Any,5}, u,
                                           element, mesh::TreeMesh{3},
                                           nonconservative_terms::Val{false},
                                           equations::CompressibleEulerEquations3D,
                                           volume_flux::typeof(flux_shima_etal),
                                           dg::DGSEM, cache,
                                           alpha=true)
  # true * [some floating point value] == [exactly the same floating point value]
  # This can (hopefully) be optimized away due to constant propagation.
  @unpack derivative_split = dg.basis

  # x direction
  # @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
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

    # for ii in eachnode(dg)
    for ii in indices(u, 2)
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
      kin_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      pv1_avg = 0.5 * (p_ll * v1_rr + p_rr * v1_ll)
      f1 = rho_avg * v1_avg
      f2 = rho_avg * v1_avg * v1_avg + p_avg
      f3 = rho_avg * v1_avg * v2_avg
      f4 = rho_avg * v1_avg * v3_avg
      f5 = p_avg * v1_avg * equations.inv_gamma_minus_one + rho_avg * v1_avg * kin_avg + pv1_avg

      factor = alpha * derivative_split[i, ii]
      du[1, i, j, k, element] += factor * f1
      du[2, i, j, k, element] += factor * f2
      du[3, i, j, k, element] += factor * f3
      du[4, i, j, k, element] += factor * f4
      du[5, i, j, k, element] += factor * f5
    end
  end

  # y direction
  # @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
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

    # for jj in eachnode(dg)
    for jj in indices(u, 3)
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
      kin_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      pv2_avg = 0.5 * (p_ll*v2_rr + p_rr*v2_ll)
      f1 = rho_avg * v2_avg
      f2 = rho_avg * v2_avg * v1_avg
      f3 = rho_avg * v2_avg * v2_avg + p_avg
      f4 = rho_avg * v2_avg * v3_avg
      f5 = p_avg*v2_avg * equations.inv_gamma_minus_one + rho_avg*v2_avg*kin_avg + pv2_avg

      factor = alpha * derivative_split[j, jj]
      du[1, i, j, k, element] += factor * f1
      du[2, i, j, k, element] += factor * f2
      du[3, i, j, k, element] += factor * f3
      du[4, i, j, k, element] += factor * f4
      du[5, i, j, k, element] += factor * f5
    end
  end

  # z direction
  # @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
  @turbo for k in indices((u, du), 4), j in indices((u, du), 3), i in indices((u, du), 2)
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

    # for kk in eachnode(dg)
    for kk in indices(u, 4)
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
      kin_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      pv3_avg = 0.5 * (p_ll*v3_rr + p_rr*v3_ll)
      f1 = rho_avg * v3_avg
      f2 = rho_avg * v3_avg * v1_avg
      f3 = rho_avg * v3_avg * v2_avg
      f4 = rho_avg * v3_avg * v3_avg + p_avg
      f5 = p_avg*v3_avg * equations.inv_gamma_minus_one + rho_avg*v3_avg*kin_avg + pv3_avg

      factor = alpha * derivative_split[k, kk]
      du[1, i, j, k, element] += factor * f1
      du[2, i, j, k, element] += factor * f2
      du[3, i, j, k, element] += factor * f3
      du[4, i, j, k, element] += factor * f4
      du[5, i, j, k, element] += factor * f5
    end
  end
end
