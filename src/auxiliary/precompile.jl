
#=
The code contained in this file is inspired by an analysis performed
using SnoopCompile, although most of it is written manually.

This kind of analysis was performed using the following code.
```julia
using SnoopCompile
tinf = @snoopi_deep begin
  include("examples/1d/elixir_advection_extended.jl")
  include("examples/2d/elixir_euler_ec.jl")
  include("examples/3d/elixir_mhd_alfven_wave_mortar.jl")
  trixi_include(default_example())
end
ttot, pcs = SnoopCompile.parcel(tinf); @show ttot; SnoopCompile.write("dev/precompile", pcs)

```
After running the code above, SnoopCompile has generated the file
`dev/precompile/precompile_Trixi.jl`, which contains precompile statements
in the function `_precompile_`.
More information on this process can be found in the docs of SnoopCompile,
in particular at https://timholy.github.io/SnoopCompile.jl/stable/.

This kind of analysis helps finding type-unstable parts of the code, e.g.
`init_interfaces` etc. in https://github.com/trixi-framework/Trixi.jl/pull/307.
Moreover, it allows to generate precompile statements which reduce the latency
by caching type inference results.
The latency can be measured by running
```bash
julia --threads=1 -e '@time using Trixi; @time trixi_include(default_example())'
```
Alternatively, a more detailed overview of compilations times can be obtained by
```bash
julia --threads=1 --check-bounds=no dev/time_compilation.jl  | egrep "seconds"
```

We add `@assert` to the precompile statements below to make sure that we don't include
failing precompile statements, cf. https://timholy.github.io/SnoopCompile.jl/stable/snoopi/.
If any assertions below fail, it is generally safe to just disable the failing call
to precompile (or to not include precompile.jl at all).
To still get the same latency reductions, you should consider to adapt the failing precompile
statements in accordance with the changes in Trixi's source code. Please, feel free to ping
the core developers of Trixi to get help with that.
=#

# manually generated precompile statements
function _precompile_manual_()
  ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

  # Constructors: mesh
  for RealT in (Int, Float64,)
    # 1D
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Trixi.Type)),NamedTuple{(:initial_refinement_level, :n_cells_max),Tuple{Int,Int}},Type{TreeMesh},RealT,RealT})
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Trixi.Type)),NamedTuple{(:initial_refinement_level, :n_cells_max),Tuple{Int,Int}},Type{TreeMesh},Tuple{RealT},Tuple{RealT}})
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Trixi.Type)),NamedTuple{(:initial_refinement_level, :n_cells_max, :periodicity),Tuple{Int,Int,Bool}},Type{TreeMesh},RealT,RealT})

    # 2D
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Trixi.Type)),NamedTuple{(:initial_refinement_level, :n_cells_max),Tuple{Int,Int}},Type{TreeMesh},Tuple{RealT,RealT},Tuple{RealT,RealT}})
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Trixi.Type)),NamedTuple{(:initial_refinement_level, :n_cells_max, :periodicity),Tuple{Int,Int,Bool}},Type{TreeMesh},Tuple{RealT,RealT},Tuple{RealT,RealT}})
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:initial_refinement_level, :refinement_patches, :n_cells_max), Tuple{Int, Tuple{NamedTuple{(:type, :coordinates_min, :coordinates_max), Tuple{String, Tuple{RealT, RealT}, Tuple{RealT, RealT}}}}, Int}},Type{TreeMesh},Tuple{RealT, RealT},Tuple{RealT, RealT}})

    # 3D
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Trixi.Type)),NamedTuple{(:initial_refinement_level, :n_cells_max),Tuple{Int,Int}},Type{TreeMesh},Tuple{RealT,RealT,RealT},Tuple{RealT,RealT,RealT}})
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Trixi.Type)),NamedTuple{(:initial_refinement_level, :n_cells_max, :periodicity),Tuple{Int,Int,Bool}},Type{TreeMesh},Tuple{RealT,RealT,RealT},Tuple{RealT,RealT,RealT}})
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:initial_refinement_level, :refinement_patches, :n_cells_max), Tuple{Int, Tuple{NamedTuple{(:type, :coordinates_min, :coordinates_max), Tuple{String, Tuple{RealT, RealT, RealT}, Tuple{RealT, RealT, RealT}}}}, Int}},Type{TreeMesh},Tuple{RealT, RealT, RealT},Tuple{RealT, RealT, RealT}})
  end
  for TreeType in (SerialTree, ParallelTree), NDIMS in 1:3
    @assert Base.precompile(Tuple{typeof(initialize!),TreeMesh{NDIMS,TreeType{NDIMS}},Int,Tuple{},Tuple{}})
    for RealT in (Float64,)
      @assert Base.precompile(Tuple{typeof(initialize!),TreeMesh{3, SerialTree{3}},Int64,Tuple{NamedTuple{(:type, :coordinates_min, :coordinates_max), Tuple{String, NTuple{NDIMS, RealT}, NTuple{NDIMS, RealT}}}},Tuple{}})
    end
    @assert Base.precompile(Tuple{typeof(save_mesh_file),TreeMesh{NDIMS,TreeType{NDIMS}},String})
  end

  # Constructors: linear advection
  for RealT in (Float64,)
    @assert Base.precompile(Tuple{Type{LinearScalarAdvectionEquation1D},RealT})
    @assert Base.precompile(Tuple{Type{LinearScalarAdvectionEquation2D},RealT,RealT})
    @assert Base.precompile(Tuple{Type{LinearScalarAdvectionEquation2D},Tuple{RealT,RealT}})
    @assert Base.precompile(Tuple{Type{LinearScalarAdvectionEquation3D},RealT,RealT,RealT})
    @assert Base.precompile(Tuple{Type{LinearScalarAdvectionEquation3D},Tuple{RealT,RealT,RealT}})
  end

  # Constructors: hyperbolic diffusion
  for RealT in (Float64,)
    @assert Base.precompile(Tuple{Type{HyperbolicDiffusionEquations1D},})
    @assert Base.precompile(Tuple{Type{HyperbolicDiffusionEquations2D},})
    @assert Base.precompile(Tuple{Type{HyperbolicDiffusionEquations3D},})
  end

  # Constructors: Euler
  for RealT in (Float64,)
    @assert Base.precompile(Tuple{Type{CompressibleEulerEquations1D},RealT})
    @assert Base.precompile(Tuple{Type{CompressibleEulerEquations2D},RealT})
    @assert Base.precompile(Tuple{Type{CompressibleEulerEquations3D},RealT})
  end

  # Constructors: MHD
  for RealT in (Float64,)
    @assert Base.precompile(Tuple{Type{IdealGlmMhdEquations1D},RealT})
    @assert Base.precompile(Tuple{Type{IdealGlmMhdEquations2D},RealT})
    @assert Base.precompile(Tuple{Type{IdealGlmMhdEquations3D},RealT})
  end

  # Constructors: LBM
  for RealT in (Float64,)
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:Ma, :Re), Tuple{RealT, RealT}},Type{LatticeBoltzmannEquations2D}})
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:Ma, :Re), Tuple{RealT, Int}},Type{LatticeBoltzmannEquations2D}})
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:Ma, :Re), Tuple{RealT, RealT}},Type{LatticeBoltzmannEquations3D}})
    @assert Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:Ma, :Re), Tuple{RealT, Int}},Type{LatticeBoltzmannEquations3D}})
  end

  # Constructors of the basis are inherently type-unstable since we pass integers
  # and use their values as parameters of static arrays.
  # Nevertheless, we can still precompile methods used to construct the bases.
  Base.precompile(Tuple{Type{LobattoLegendreBasis},Int})
  for RealT in (Float64,)
    Base.precompile(Tuple{Type{LobattoLegendreBasis},RealT,Int})
    @assert Base.precompile(Tuple{typeof(calc_dhat),Vector{RealT},Vector{RealT}})
    @assert Base.precompile(Tuple{typeof(calc_dsplit),Vector{RealT},Vector{RealT}})
    @assert Base.precompile(Tuple{typeof(polynomial_derivative_matrix),Vector{RealT}})
    @assert Base.precompile(Tuple{typeof(polynomial_interpolation_matrix),Vector{RealT},Vector{RealT}})
    @assert Base.precompile(Tuple{typeof(barycentric_weights),Vector{RealT}})
    @assert Base.precompile(Tuple{typeof(calc_lhat),RealT,Vector{RealT},Vector{RealT}})
    @assert Base.precompile(Tuple{typeof(lagrange_interpolating_polynomials),RealT,Vector{RealT},Vector{RealT}})
    @assert Base.precompile(Tuple{typeof(calc_q_and_l),Int,RealT})
    @assert Base.precompile(Tuple{typeof(legendre_polynomial_and_derivative),Int,RealT})
    @assert Base.precompile(Tuple{typeof(vandermonde_legendre),Vector{RealT}})
  end
  @assert Base.precompile(Tuple{typeof(gauss_lobatto_nodes_weights),Int})
  @assert Base.precompile(Tuple{typeof(gauss_nodes_weights),Int})
  @assert Base.precompile(Tuple{typeof(calc_forward_upper),Int})
  @assert Base.precompile(Tuple{typeof(calc_forward_lower),Int})
  @assert Base.precompile(Tuple{typeof(calc_reverse_upper),Int,Val{:gauss}})
  @assert Base.precompile(Tuple{typeof(calc_reverse_lower),Int,Val{:gauss}})
  @assert Base.precompile(Tuple{typeof(calc_reverse_upper),Int,Val{:gauss_lobatto}})
  @assert Base.precompile(Tuple{typeof(calc_reverse_lower),Int,Val{:gauss_lobatto}})

  # init_elements, interfaces, etc.
  for RealT in (Float64,), polydeg in 1:7, nvars in 1:5
    uEltype = RealT
    nnodes_ = polydeg + 1

    # 1D, serial
    @assert Base.precompile(Tuple{typeof(init_boundaries),Array{Int,1},TreeMesh{1,SerialTree{1}},ElementContainer1D{RealT,uEltype,nvars,polydeg}})
    @assert Base.precompile(Tuple{typeof(init_interfaces),Array{Int,1},TreeMesh{1,SerialTree{1}},ElementContainer1D{RealT,uEltype,nvars,polydeg}})

    # 2D, serial
    @assert Base.precompile(Tuple{typeof(init_boundaries),Array{Int,1},TreeMesh{2,SerialTree{2}},ElementContainer2D{RealT,uEltype,nvars,polydeg}})
    @assert Base.precompile(Tuple{typeof(init_interfaces),Array{Int,1},TreeMesh{2,SerialTree{2}},ElementContainer2D{RealT,uEltype,nvars,polydeg}})
    @assert Base.precompile(Tuple{typeof(init_mortars!),L2MortarContainer2D{uEltype,nvars,polydeg},ElementContainer2D{RealT,uEltype,nvars,polydeg},TreeMesh{2, SerialTree{2}}})

    # 2D, parallel
    @assert Base.precompile(Tuple{typeof(init_boundaries),Array{Int,1},TreeMesh{2,ParallelTree{2}},ElementContainer2D{RealT,uEltype,nvars,polydeg}})
    @assert Base.precompile(Tuple{typeof(init_interfaces),Array{Int,1},TreeMesh{2,ParallelTree{2}},ElementContainer2D{RealT,uEltype,nvars,polydeg}})
    @assert Base.precompile(Tuple{typeof(init_mortars!),L2MortarContainer2D{uEltype,nvars,polydeg},ElementContainer2D{RealT,uEltype,nvars,polydeg},TreeMesh{2, ParallelTree{2}}})
    @assert Base.precompile(Tuple{typeof(init_mpi_interfaces),Array{Int,1},TreeMesh{2,ParallelTree{2}},ElementContainer2D{RealT,uEltype,nvars,polydeg}})

    # 3D, serial
    @assert Base.precompile(Tuple{typeof(init_boundaries),Array{Int,1},TreeMesh{3,SerialTree{3}},ElementContainer3D{RealT,uEltype,nvars,polydeg}})
    @assert Base.precompile(Tuple{typeof(init_interfaces),Array{Int,1},TreeMesh{3,SerialTree{3}},ElementContainer3D{RealT,uEltype,nvars,polydeg}})
    @assert Base.precompile(Tuple{typeof(init_mortars!),L2MortarContainer3D{uEltype,nvars,polydeg},ElementContainer3D{RealT,uEltype,nvars,polydeg},TreeMesh{3, SerialTree{3}}})
  end

  # Show/display methods: general
  Base.precompile(Tuple{typeof(summary_box),IO,String,Vector{Pair{String, Any}}})

  # Show/display methods: equations
  Base.precompile(Tuple{typeof(show),typeof(stdout),MIME{Symbol("text/plain")},AbstractEquations})
  Base.precompile(Tuple{typeof(show),IOContext{typeof(stdout)},MIME{Symbol("text/plain")},AbstractEquations})

  # Show/display methods: meshes
  for NDIMS in 1:3
    # serial
    @assert Base.precompile(Tuple{typeof(show),typeof(stdout),TreeMesh{NDIMS, SerialTree{NDIMS}}})
    @assert Base.precompile(Tuple{typeof(show),IOContext{typeof(stdout)},TreeMesh{NDIMS, SerialTree{NDIMS}}})
    @assert Base.precompile(Tuple{typeof(show),IOContext{typeof(stdout)},MIME{Symbol("text/plain")},TreeMesh{NDIMS, SerialTree{NDIMS}}})
    @assert Base.precompile(Tuple{typeof(summary_line),IOContext{typeof(stdout)},String,TreeMesh{NDIMS, SerialTree{NDIMS}}})
    # parallel
    @assert Base.precompile(Tuple{typeof(show),typeof(stdout),TreeMesh{NDIMS, ParallelTree{NDIMS}}})
    @assert Base.precompile(Tuple{typeof(show),IOContext{typeof(stdout)},MIME"text/plain",TreeMesh{NDIMS, ParallelTree{NDIMS}}})
    @assert Base.precompile(Tuple{typeof(summary_line),IOContext{typeof(stdout)},String,TreeMesh{NDIMS, ParallelTree{NDIMS}}})
  end

  # Show/display methods: DG
  Base.precompile(Tuple{typeof(show),typeof(stdout),DG})
  Base.precompile(Tuple{typeof(show),IOContext{typeof(stdout)},MIME{Symbol("text/plain")},DG})
  Base.precompile(Tuple{typeof(show),typeof(stdout),SemidiscretizationHyperbolic})
  Base.precompile(Tuple{typeof(show),IOContext{typeof(stdout)},MIME{Symbol("text/plain")},SemidiscretizationHyperbolic})

  # multiply_dimensionwise! as used in the analysis callback
  for RealT in (Float64,)
    @assert Base.precompile(Tuple{typeof(multiply_dimensionwise!),Array{RealT, 2},Matrix{RealT},SubArray{RealT, 2, Array{RealT, 3}, Tuple{Base.Slice{Base.OneTo{Int}}, Base.Slice{Base.OneTo{Int}}, Int}, true}})
    @assert Base.precompile(Tuple{typeof(multiply_dimensionwise!),Array{RealT, 3},Matrix{RealT},SubArray{RealT, 3, Array{RealT, 4}, Tuple{Base.Slice{Base.OneTo{Int}}, Base.Slice{Base.OneTo{Int}}, Base.Slice{Base.OneTo{Int}}, Int}, true},Array{RealT, 3}})
    @assert Base.precompile(Tuple{typeof(multiply_dimensionwise!),Array{RealT, 4},Matrix{RealT},SubArray{RealT, 4, Array{RealT, 5}, Tuple{Base.Slice{Base.OneTo{Int}}, Base.Slice{Base.OneTo{Int}}, Base.Slice{Base.OneTo{Int}}, Base.Slice{Base.OneTo{Int}}, Int}, true},Array{RealT, 4},Array{RealT, 4}})
  end

  # infrastructure, special elixirs
  @assert Base.precompile(Tuple{typeof(trixi_include),String})
  @assert Base.precompile(Tuple{typeof(default_example)})

  return nothing
end
