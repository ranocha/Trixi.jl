
#=

function run_benchmarks_shima_etal(polydegs)
  baseline = Vector{Float64}()
  permuted = Vector{Float64}()

  for polydeg in polydegs
    @show polydeg

    redirect_stdout(devnull) do
      trixi_include(joinpath(examples_dir(), "structured_3d_dgsem", "elixir_euler_ec.jl"),
        volume_flux=flux_shima_etal, surface_flux=flux_shima_etal,
        tspan=(0.0, 0.0), cells_per_dimension=(1, 1, 1),
        polydeg=polydeg)
    end
    u_ode = copy(ode.u0)
    du_ode = zero(u_ode)

    GC.@preserve u_ode du_ode begin
      u = Trixi.wrap_array(u_ode, semi)
      du = Trixi.wrap_array(du_ode, semi)

      result = @benchmark Trixi.split_form_kernel!(
        $du, $u, $(Ref(1))[], $(semi.mesh),
        $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
        $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))
      push!(baseline, time(minimum(result)))

      result = @benchmark Trixi.flux_differencing_kernel_permuted!(
        $du, $u, $(Ref(1))[], $(semi.mesh),
        $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
        $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))
      push!(permuted, time(minimum(result)))
    end
  end

  hcat(polydegs, baseline,
       permuted, permuted ./ baseline)
end

results_shima_etal = run_benchmarks_shima_etal(1:15)
# 15×4 Matrix{Float64}:
#   1.0     220.77     160.966  0.729112
#   2.0    1020.1      511.594  0.501513
#   3.0    3125.75    1073.9    0.343566
#   4.0    7580.75    3198.25   0.421891
#   5.0   15747.0     5709.17   0.362556
#   6.0   29120.0    11502.0    0.394986
#   7.0   49876.0    19443.0    0.389827
#   8.0   80628.0    30683.0    0.38055
#   9.0  124761.0    45374.0    0.363687
#  10.0  187370.0    70199.0    0.374654
#  11.0  257747.0    89271.0    0.346351
#  12.0  368058.0   140885.0    0.382779
#  13.0  488652.0   181817.0    0.372079
#  14.0  652503.0   254041.0    0.389333
#  15.0  887587.0   532083.0    0.599471

let results = results_shima_etal
  fig_abs = plot(xguide="Polynomial degree", yguide="Runtime",
                 yscale=:log10, legend=:bottomright)
    scatter!(results[:, 1], results[:, 3], label="Permuted")
    scatter!(results[:, 1], results[:, 2], label="Baseline")
  fig_rel = plot(xguide="Polynomial degree", yguide="Relative runtime")
    scatter!(results[:, 1], results[:, 4], label="Permuted")
    plot!(yrange=(0.2, 1.1))
  plot(fig_abs, fig_rel, layout=(2,1), dpi=300, title="Shima et al. (curved)")
end


=#


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
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
     StaticInt(ndims(mesh))))
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    jk = j + nnodes(dg) * (k- 1)
    contravariant_vectors_x[jk, i, 1] = contravariant_vectors[1, 1, i, j, k, element]
    contravariant_vectors_x[jk, i, 2] = contravariant_vectors[2, 1, i, j, k, element]
    contravariant_vectors_x[jk, i, 3] = contravariant_vectors[3, 1, i, j, k, element]
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
        contravariant_vectors_x[jk, i, 1] + contravariant_vectors_x[jk, ii, 1])
      normal_direction_2 = 0.5 * (
        contravariant_vectors_x[jk, i, 2] + contravariant_vectors_x[jk, ii, 2])
      normal_direction_3 = 0.5 * (
        contravariant_vectors_x[jk, i, 3] + contravariant_vectors_x[jk, ii, 3])

      v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 + v3_ll * normal_direction_3
      v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 + v3_rr * normal_direction_3

      # Compute necessary mean values
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

      # Compute necessary mean values
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

        # Compute necessary mean values
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



#=

function run_benchmarks_ranocha(polydegs)
  baseline = Vector{Float64}()
  permuted_log = Vector{Float64}()

  for polydeg in polydegs
    @show polydeg

    redirect_stdout(devnull) do
      trixi_include(joinpath(examples_dir(), "structured_3d_dgsem", "elixir_euler_ec.jl"),
        volume_flux=flux_ranocha, surface_flux=flux_ranocha,
        tspan=(0.0, 0.0), cells_per_dimension=(1, 1, 1),
        polydeg=polydeg)
    end
    u_ode = copy(ode.u0)
    du_ode = zero(u_ode)

    GC.@preserve u_ode du_ode begin
      u = Trixi.wrap_array(u_ode, semi)
      du = Trixi.wrap_array(du_ode, semi)

      result = @benchmark Trixi.split_form_kernel!(
        $du, $u, $(Ref(1))[], $(semi.mesh),
        $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
        $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))
      push!(baseline, time(minimum(result)))

      result = @benchmark Trixi.flux_differencing_kernel_permuted_log!(
        $du, $u, $(Ref(1))[], $(semi.mesh),
        $(Trixi.have_nonconservative_terms(semi.equations)), $(semi.equations),
        $(semi.solver.volume_integral.volume_flux), $(semi.solver), $(semi.cache))
      push!(permuted_log, time(minimum(result)))
    end
  end

  hcat(polydegs, baseline,
       permuted_log, permuted_log ./ baseline)
end

results_ranocha = run_benchmarks_ranocha(1:15)
# 7×4 Matrix{Float64}:
#  1.0    349.098    259.504  0.743357
#  2.0   1705.8      865.895  0.507618
#  3.0   5406.67    1999.4    0.369803
#  4.0  13268.0     6403.6    0.482635
#  5.0  29053.0    10749.0    0.369979
#  6.0  52635.0    20070.0    0.381305
#  7.0  89433.0    34274.0    0.383237

let results = results_ranocha
  fig_abs = plot(xguide="Polynomial degree", yguide="Runtime",
                 yscale=:log10, legend=:bottomright);
    scatter!(results[:, 1], results[:, 3], label="Permuted, log")
    scatter!(results[:, 1], results[:, 2], label="Baseline")
  fig_rel = plot(xguide="Polynomial degree", yguide="Relative runtime")
    scatter!(results[:, 1], results[:, 4], label="Permuted, log")
    plot!(yrange=(0.2, 1.1))
  plot(fig_abs, fig_rel, layout=(2,1), dpi=300, title="EC flux (curved)")
end


=#



# pre-compute primitive variables, uses triangular loops without turbo
# after permuting the memory layout, and pre-computes logs
# @inline function split_form_kernel!(_du, u_cons,
@inline function flux_differencing_kernel_permuted_log!(_du, u_cons,
                                                    element, mesh::Union{StructuredMesh{3}, P4estMesh{3}},
                                                    nonconservative_terms::Val{false},
                                                    equations::CompressibleEulerEquations3D,
                                                    volume_flux::typeof(flux_ranocha),
                                                    dg::DGSEM, cache,
                                                    alpha=true)
  # true * [some floating point value] == [exactly the same floating point value]
  # This can (hopefully) be optimized away due to constant propagation.
  @unpack derivative_split, derivative_split_transpose = dg.basis
  @unpack contravariant_vectors = cache.elements


  u_prim = StrideArray{eltype(u_cons)}(undef,
    (ntuple(_ -> StaticInt(nnodes(dg)), ndims(mesh))...,
     StaticInt(nvariables(equations) + 2))) # We also compute logs

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
    u_prim[i, j, k, 6] = log(rho)
    u_prim[i, j, k, 7] = log(p)
  end


  # x direction
  du_permuted = StrideArray{eltype(u_cons)}(undef,
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
     StaticInt(nvariables(equations))))

  u_prim_permuted = StrideArray{eltype(u_cons)}(undef,
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
     StaticInt(nvariables(equations) + 2)))

  @turbo for v in indices(u_prim, 4),
             k in eachnode(dg),
             j in eachnode(dg),
             i in eachnode(dg)
    jk = j + nnodes(dg) * (k- 1)
    u_prim_permuted[jk, i, v] = u_prim[i, j, k, v]
  end
  fill!(du_permuted, zero(eltype(du_permuted)))

  contravariant_vectors_x = StrideArray{eltype(contravariant_vectors)}(undef,
    (StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
     StaticInt(ndims(mesh))))
  @turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
    jk = j + nnodes(dg) * (k- 1)
    contravariant_vectors_x[jk, i, 1] = contravariant_vectors[1, 1, i, j, k, element]
    contravariant_vectors_x[jk, i, 2] = contravariant_vectors[2, 1, i, j, k, element]
    contravariant_vectors_x[jk, i, 3] = contravariant_vectors[3, 1, i, j, k, element]
  end

  for i in eachnode(dg), ii in (i+1):nnodes(dg)
    @turbo for jk in Base.OneTo(nnodes(dg)^2)
      rho_ll = u_prim_permuted[jk, i, 1]
      v1_ll  = u_prim_permuted[jk, i, 2]
      v2_ll  = u_prim_permuted[jk, i, 3]
      v3_ll  = u_prim_permuted[jk, i, 4]
      p_ll   = u_prim_permuted[jk, i, 5]
      log_rho_ll = u_prim_permuted[jk, i, 6]
      log_p_ll   = u_prim_permuted[jk, i, 7]

      rho_rr = u_prim_permuted[jk, ii, 1]
      v1_rr  = u_prim_permuted[jk, ii, 2]
      v2_rr  = u_prim_permuted[jk, ii, 3]
      v3_rr  = u_prim_permuted[jk, ii, 4]
      p_rr   = u_prim_permuted[jk, ii, 5]
      log_rho_rr = u_prim_permuted[jk, ii, 6]
      log_p_rr   = u_prim_permuted[jk, ii, 7]

      normal_direction_1 = 0.5 * (
        contravariant_vectors_x[jk, i, 1] + contravariant_vectors_x[jk, ii, 1])
      normal_direction_2 = 0.5 * (
        contravariant_vectors_x[jk, i, 2] + contravariant_vectors_x[jk, ii, 2])
      normal_direction_3 = 0.5 * (
        contravariant_vectors_x[jk, i, 3] + contravariant_vectors_x[jk, ii, 3])

      v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 + v3_ll * normal_direction_3
      v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 + v3_rr * normal_direction_3

      # Compute necessary mean values
      # rho_mean = ln_mean(rho_ll, rho_rr)
      x1 = rho_ll
      log_x1 = log_rho_ll
      y1 = rho_rr
      log_y1 = log_rho_rr
      x1_plus_y1 = x1 + y1
      y1_minus_x1 = y1 - x1
      z1 = y1_minus_x1^2 / x1_plus_y1^2
      special_path1 = x1_plus_y1 / (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1)))
      regular_path1 = y1_minus_x1 / (log_y1 - log_x1)
      rho_mean = ifelse(z1 < 1.0e-4, special_path1, regular_path1)

      # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
      # in exact arithmetic since
      #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
      #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
      # inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
      x2 = rho_ll * p_rr
      log_x2 = log_rho_ll + log_p_rr
      y2 = rho_rr * p_ll
      log_y2 = log_rho_rr + log_p_ll
      x2_plus_y2 = x2 + y2
      y2_minus_x2 = y2 - x2
      z2 = y2_minus_x2^2 / x2_plus_y2^2
      special_path2 = (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1))) / x2_plus_y2
      regular_path2 = (log_y2 - log_x2) / y2_minus_x2
      inv_rho_p_mean = p_ll * p_rr * ifelse(z2 < 1.0e-4, special_path2, regular_path2)

      v1_avg = 0.5 * (v1_ll + v1_rr)
      v2_avg = 0.5 * (v2_ll + v2_rr)
      v3_avg = 0.5 * (v3_ll + v3_rr)
      p_avg  = 0.5 * (p_ll + p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      # Calculate fluxes depending on normal_direction
      f1 = rho_mean * 0.5 * (v_dot_n_ll + v_dot_n_rr)
      f2 = f1 * v1_avg + p_avg * normal_direction_1
      f3 = f1 * v2_avg + p_avg * normal_direction_2
      f4 = f1 * v3_avg + p_avg * normal_direction_3
      f5 = ( f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one )
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
      log_rho_ll = u_prim[i, j, k, 6]
      log_p_ll   = u_prim[i, j, k, 7]

      rho_rr = u_prim[i, jj, k, 1]
      v1_rr  = u_prim[i, jj, k, 2]
      v2_rr  = u_prim[i, jj, k, 3]
      v3_rr  = u_prim[i, jj, k, 4]
      p_rr   = u_prim[i, jj, k, 5]
      log_rho_rr = u_prim[i, jj, k, 6]
      log_p_rr   = u_prim[i, jj, k, 7]

      normal_direction_1 = 0.5 * (
        contravariant_vectors_y[i, j, k, 1] + contravariant_vectors_y[i, jj, k, 1])
      normal_direction_2 = 0.5 * (
        contravariant_vectors_y[i, j, k, 2] + contravariant_vectors_y[i, jj, k, 2])
      normal_direction_3 = 0.5 * (
        contravariant_vectors_y[i, j, k, 3] + contravariant_vectors_y[i, jj, k, 3])

      v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 + v3_ll * normal_direction_3
      v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 + v3_rr * normal_direction_3

      # Compute necessary mean values
      # rho_mean = ln_mean(rho_ll, rho_rr)
      x1 = rho_ll
      log_x1 = log_rho_ll
      y1 = rho_rr
      log_y1 = log_rho_rr
      x1_plus_y1 = x1 + y1
      y1_minus_x1 = y1 - x1
      z1 = y1_minus_x1^2 / x1_plus_y1^2
      special_path1 = x1_plus_y1 / (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1)))
      regular_path1 = y1_minus_x1 / (log_y1 - log_x1)
      rho_mean = ifelse(z1 < 1.0e-4, special_path1, regular_path1)

      # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
      # in exact arithmetic since
      #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
      #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
      # inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
      x2 = rho_ll * p_rr
      log_x2 = log_rho_ll + log_p_rr
      y2 = rho_rr * p_ll
      log_y2 = log_rho_rr + log_p_ll
      x2_plus_y2 = x2 + y2
      y2_minus_x2 = y2 - x2
      z2 = y2_minus_x2^2 / x2_plus_y2^2
      special_path2 = (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1))) / x2_plus_y2
      regular_path2 = (log_y2 - log_x2) / y2_minus_x2
      inv_rho_p_mean = p_ll * p_rr * ifelse(z2 < 1.0e-4, special_path2, regular_path2)

      v1_avg = 0.5 * (v1_ll + v1_rr)
      v2_avg = 0.5 * (v2_ll + v2_rr)
      v3_avg = 0.5 * (v3_ll + v3_rr)
      p_avg  = 0.5 * (p_ll + p_rr)
      velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

      # Calculate fluxes depending on normal_direction
      f1 = rho_mean * 0.5 * (v_dot_n_ll + v_dot_n_rr)
      f2 = f1 * v1_avg + p_avg * normal_direction_1
      f3 = f1 * v2_avg + p_avg * normal_direction_2
      f4 = f1 * v3_avg + p_avg * normal_direction_3
      f5 = ( f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one )
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
       StaticInt(nvariables(equations) + 2)))

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
        log_rho_ll = u_prim_perm[ij, k, 6]
        log_p_ll   = u_prim_perm[ij, k, 7]

        rho_rr = u_prim_perm[ij, kk, 1]
        v1_rr  = u_prim_perm[ij, kk, 2]
        v2_rr  = u_prim_perm[ij, kk, 3]
        v3_rr  = u_prim_perm[ij, kk, 4]
        p_rr   = u_prim_perm[ij, kk, 5]
        log_rho_rr = u_prim_perm[ij, kk, 6]
        log_p_rr   = u_prim_perm[ij, kk, 7]

        normal_direction_1 = 0.5 * (
          contravariant_vectors_z[ij, k, 1] + contravariant_vectors_z[ij, kk, 1])
        normal_direction_2 = 0.5 * (
          contravariant_vectors_z[ij, k, 2] + contravariant_vectors_z[ij, kk, 2])
        normal_direction_3 = 0.5 * (
          contravariant_vectors_z[ij, k, 3] + contravariant_vectors_z[ij, kk, 3])

        v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 + v3_ll * normal_direction_3
        v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 + v3_rr * normal_direction_3

        # Compute necessary mean values
        # rho_mean = ln_mean(rho_ll, rho_rr)
        x1 = rho_ll
        log_x1 = log_rho_ll
        y1 = rho_rr
        log_y1 = log_rho_rr
        x1_plus_y1 = x1 + y1
        y1_minus_x1 = y1 - x1
        z1 = y1_minus_x1^2 / x1_plus_y1^2
        special_path1 = x1_plus_y1 / (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1)))
        regular_path1 = y1_minus_x1 / (log_y1 - log_x1)
        rho_mean = ifelse(z1 < 1.0e-4, special_path1, regular_path1)

        # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
        # in exact arithmetic since
        #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
        #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
        # inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
        x2 = rho_ll * p_rr
        log_x2 = log_rho_ll + log_p_rr
        y2 = rho_rr * p_ll
        log_y2 = log_rho_rr + log_p_ll
        x2_plus_y2 = x2 + y2
        y2_minus_x2 = y2 - x2
        z2 = y2_minus_x2^2 / x2_plus_y2^2
        special_path2 = (2 + z1*(0.666666666666667 + z1*(0.4 + 0.285714285714286*z1))) / x2_plus_y2
        regular_path2 = (log_y2 - log_x2) / y2_minus_x2
        inv_rho_p_mean = p_ll * p_rr * ifelse(z2 < 1.0e-4, special_path2, regular_path2)

        v1_avg = 0.5 * (v1_ll + v1_rr)
        v2_avg = 0.5 * (v2_ll + v2_rr)
        v3_avg = 0.5 * (v3_ll + v3_rr)
        p_avg  = 0.5 * (p_ll + p_rr)
        velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

        # Calculate fluxes depending on normal_direction
        f1 = rho_mean * 0.5 * (v_dot_n_ll + v_dot_n_rr)
        f2 = f1 * v1_avg + p_avg * normal_direction_1
        f3 = f1 * v2_avg + p_avg * normal_direction_2
        f4 = f1 * v3_avg + p_avg * normal_direction_3
        f5 = ( f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one )
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
