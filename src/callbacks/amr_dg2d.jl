
# TODO: Taal dimension agnostic
function (amr_callback::AMRCallback)(u_ode::AbstractVector, mesh::TreeMesh,
                                     equations, dg::DG, cache;
                                     only_refine=false, only_coarsen=false)
  @unpack indicator, adaptor = amr_callback

  u = wrap_array(u_ode, mesh, equations, dg, cache)
  lambda = @timeit timer() "indicator" indicator(u, mesh, equations, dg, cache)

  leaf_cell_ids = leaf_cells(mesh.tree)
  @assert length(lambda) == length(leaf_cell_ids) ("Indicator (length = $(length(lambda))) and leaf cell (length = $(length(leaf_cell_ids))) arrays have different length")

  to_refine  = Int[]
  to_coarsen = Int[]
  for element in eachelement(dg, cache)
    indicator_value = lambda[element]
    if indicator_value > 0
      push!(to_refine, leaf_cell_ids[element])
    elseif indicator_value < 0
      push!(to_coarsen, leaf_cell_ids[element])
    end
  end

  # TODO: Taal debug
  pointer_before__refine = pointer(u_ode)
  length_before__refine  = length(u_ode)
  let u = wrap_array(u_ode, mesh, equations, dg, cache)
    unstable_idx = findfirst(u -> abs(u) < 0.01, u)
    if unstable_idx !== nothing
      @info "AMR, before refine!" unstable_idx u[unstable_idx] pointer_before__refine length_before__refine
    end
    # unstable_idx_ode = findfirst(u -> abs(u) < 0.01, u_ode)
    # if unstable_idx !== nothing
    #   @info "AMR, before refine!" unstable_idx_ode u_ode[unstable_idx_ode]
    # end
  end

  @timeit timer() "refine" if !only_coarsen && !isempty(to_refine)
    # refine mesh
    refined_original_cells = @timeit timer() "mesh" refine!(mesh.tree, to_refine)

    # refine solver
    @timeit timer() "solver" refine!(u_ode, adaptor, mesh, equations, dg, cache, refined_original_cells)
    # TODO: Taal implement, passive solvers?
    # if !isempty(passive_solvers)
    #   @timeit timer() "passive solvers" for ps in passive_solvers
    #     refine!(ps, mesh, refined_original_cells)
    #   end
    # end
  else
    # If there is nothing to refine, create empty array for later use
    refined_original_cells = Int[]
  end

  # TODO: Taal debug
  pointer_before_coarsen = pointer(u_ode)
  length_before_coarsen  = length(u_ode)
  let u = wrap_array(u_ode, mesh, equations, dg, cache)
    unstable_idx = findfirst(u -> abs(u) < 0.01, u)
    if unstable_idx !== nothing
      @info "AMR, before coarsen!" unstable_idx u[unstable_idx] pointer_before__refine pointer_before_coarsen length_before__refine length_before_coarsen
    end
    # unstable_idx_ode = findfirst(u -> abs(u) < 0.01, u_ode)
    # if unstable_idx !== nothing
    #   @info "AMR, before coarsen!" unstable_idx_ode u_ode[unstable_idx_ode]
    # end
  end

  @timeit timer() "coarsen" if !only_refine && !isempty(to_coarsen)
    # Since the cells may have been shifted due to refinement, first we need to
    # translate the old cell ids to the new cell ids
    if !isempty(to_coarsen)
      to_coarsen = original2refined(to_coarsen, refined_original_cells, mesh)
    end

    # Next, determine the parent cells from which the fine cells are to be
    # removed, since these are needed for the coarsen! function. However, since
    # we only want to coarsen if *all* child cells are marked for coarsening,
    # we count the coarsening indicators for each parent cell and only coarsen
    # if all children are marked as such (i.e., where the count is 2^ndims). At
    # the same time, check if a cell is marked for coarsening even though it is
    # *not* a leaf cell -> this can only happen if it was refined due to 2:1
    # smoothing during the preceding refinement operation.
    parents_to_coarsen = zeros(Int, length(mesh.tree))
    for cell_id in to_coarsen
      # If cell has no parent, it cannot be coarsened
      if !has_parent(mesh.tree, cell_id)
        continue
      end

      # If cell is not leaf (anymore), it cannot be coarsened
      if !is_leaf(mesh.tree, cell_id)
        continue
      end

      # Increase count for parent cell
      parent_id = mesh.tree.parent_ids[cell_id]
      parents_to_coarsen[parent_id] += 1
    end

    # Extract only those parent cells for which all children should be coarsened
    to_coarsen = collect(1:length(parents_to_coarsen))[parents_to_coarsen .== 2^ndims(mesh)]

    # Finally, coarsen mesh
    coarsened_original_cells = @timeit timer() "mesh" coarsen!(mesh.tree, to_coarsen)

    # Convert coarsened parent cell ids to the list of child cell ids that have
    # been removed, since this is the information that is expected by the solver
    removed_child_cells = zeros(Int, n_children_per_cell(mesh.tree) * length(coarsened_original_cells))
    for (index, coarse_cell_id) in enumerate(coarsened_original_cells)
      for child in 1:n_children_per_cell(mesh.tree)
        removed_child_cells[n_children_per_cell(mesh.tree) * (index-1) + child] = coarse_cell_id + child
      end
    end

    # coarsen solver
    @timeit timer() "solver" coarsen!(u_ode, adaptor, mesh, equations, dg, cache, removed_child_cells)
    # TODO: Taal implement, passive solvers?
    # if !isempty(passive_solvers)
    #   @timeit timer() "passive solvers" for ps in passive_solvers
    #     coarsen!(ps, mesh, removed_child_cells)
    #   end
    # end
  else
    # If there is nothing to coarsen, create empty array for later use
    coarsened_original_cells = Int[]
  end

  # TODO: Taal debug
  pointer_after__coarsen = pointer(u_ode)
  length_after__coarsen  = length(u_ode)
  let u = wrap_array(u_ode, mesh, equations, dg, cache)
    unstable_idx = findfirst(u -> abs(u) < 0.01, u)
    if unstable_idx !== nothing
      @info "AMR, after coarsen!" unstable_idx u[unstable_idx] pointer_before__refine pointer_before_coarsen pointer_after__coarsen length_before__refine length_before_coarsen length_after__coarsen isempty(to_refine) isempty(to_coarsen)
    end
    # unstable_idx_ode = findfirst(u -> abs(u) < 0.01, u_ode)
    # if unstable_idx !== nothing
    #   @info "AMR, after coarsen!" unstable_idx_ode u_ode[unstable_idx_ode]
    # end
  end

  # TODO: Taal debug
  let u = wrap_array(u_ode, mesh, equations, dg, cache)
    unstable_idx = findfirst(u -> abs(u) < 0.01, u)
    if unstable_idx !== nothing
      u[unstable_idx] = 1
      unstable_idx = findfirst(u -> abs(u) < 0.01, u)
    end
    if unstable_idx !== nothing
      @info "AMR, after setting" unstable_idx u[unstable_idx]
    end
  end

 # Return true if there were any cells coarsened or refined, otherwise false
 has_changed = !isempty(refined_original_cells) || !isempty(coarsened_original_cells)
 if has_changed # TODO: Taal decide, where shall we set this?
  # don't set it to has_changed since there can be changes from earlier calls
  mesh.unsaved_changes = true
 end

 return has_changed
end

# function original2refined(original_cell_ids, refined_original_cells, mesh) in src/amr/amr.jl


function refine!(u_ode::AbstractVector, adaptor, mesh::TreeMesh{2}, equations, dg::DGSEM, cache, cells_to_refine)
  # Return early if there is nothing to do
  if isempty(cells_to_refine)
    return
  end

  # Determine for each existing element whether it needs to be refined
  needs_refinement = falses(nelements(dg, cache))

  # The "Ref(...)" is such that we can vectorize the search but not the array that is searched
  elements_to_refine = searchsortedfirst.(Ref(cache.elements.cell_ids[1:nelements(dg, cache)]),
                                          cells_to_refine)
  needs_refinement[elements_to_refine] .= true

  # Retain current solution data
  old_n_elements = nelements(dg, cache)
  old_u_ode = copy(u_ode)
  old_u     = wrap_array(old_u_ode, mesh, equations, dg, cache)

  # Get new list of leaf cells
  leaf_cell_ids = leaf_cells(mesh.tree)

  # Initialize new elements container
  elements = init_elements(leaf_cell_ids, mesh,
                           real(dg), nvariables(equations), polydeg(dg))
  copy!(cache.elements, elements)
  @assert nelements(dg, cache) == nelements(elements) # TODO: Taal debug
  @assert nelements(dg, cache) > old_n_elements # TODO: Taal debug

  resize!(u_ode, nvariables(equations) * nnodes(dg)^ndims(mesh) * nelements(dg, cache))
  u = wrap_array(u_ode, mesh, equations, dg, cache)

  # Loop over all elements in old container and either copy them or refine them
  element_id = 1
  for old_element_id in 1:old_n_elements
    if needs_refinement[old_element_id]
      # Refine element and store solution directly in new data structure
      # TODO: Taal debug
      # refine_element!(u, element_id, old_u, old_element_id,
      #                 adaptor, equations, dg)
      for element in element_id:element_id+3, j in eachnode(dg), i in eachnode(dg), v in eachvariable(equations)
        u[v, i, j, element] = old_u[v, i, j, old_element_id]
      end
      element_id += 2^ndims(mesh)
    else
      # Copy old element data to new element container
      # @views u[:, .., element_id] .= old_u[:, .., old_element_id] # TODO: Taal debug
      for j in eachnode(dg), i in eachnode(dg), v in eachvariable(equations)
        u[v, i, j, element_id] = old_u[v, i, j, old_element_id]
      end
      element_id += 1
    end
  end
  @assert element_id == nelements(dg, cache) + 1 || element_id == nelements(dg, cache) + 2^ndims(mesh) "element_id = $element_id, nelements(dg, cache) = $(nelements(dg, cache))" # TODO: Taal debug

  # TODO: Taal performance, allow initializing the stuff in place, making use of resize!
  # Initialize new interfaces container
  interfaces = init_interfaces(leaf_cell_ids, mesh, elements,
                               real(dg), nvariables(equations), polydeg(dg))
  copy!(cache.interfaces, interfaces)

  # Initialize boundaries
  boundaries, _ = init_boundaries(leaf_cell_ids, mesh, elements,
                                  real(dg), nvariables(equations), polydeg(dg))
  copy!(cache.boundaries, boundaries)

  # Initialize new mortar containers
  mortars = init_mortars(leaf_cell_ids, mesh, elements,
                         real(dg), nvariables(equations), polydeg(dg), dg.mortar)
  copy!(cache.mortars, mortars)

  # Sanity check
  if isperiodic(mesh.tree) && nmortars(mortars) == 0
    @assert ninterfaces(interfaces) == 2 * nelements(dg, cache) ("For 2D and periodic domains and conforming elements, the number of interfaces must be twice the number of elements")
  end

  return nothing
end


# TODO: Taal compare performance of different implementations
# Refine solution data u for an element, using L2 projection (interpolation)
function refine_element!(u::AbstractArray{<:Any,4}, element_id, old_u, old_element_id,
                         adaptor::LobattoLegendreAdaptorL2, equations, dg)
  @unpack forward_upper, forward_lower = adaptor

  # Store new element ids
  lower_left_id  = element_id
  lower_right_id = element_id + 1
  upper_left_id  = element_id + 2
  upper_right_id = element_id + 3

  # TODO: Taal debug
  @assert old_element_id >= 1
  @assert size(old_u, 1) == nvariables(equations)
  @assert size(old_u, 2) == nnodes(dg)
  @assert size(old_u, 3) == nnodes(dg)
  @assert size(old_u, 4) >= old_element_id
  @assert     element_id >= 1
  @assert size(    u, 1) == nvariables(equations)
  @assert size(    u, 2) == nnodes(dg)
  @assert size(    u, 3) == nnodes(dg)
  @assert size(    u, 4) >= element_id + 3

  # Interpolate to lower left element
  for j in eachnode(dg), i in eachnode(dg)
    acc = zero(get_node_vars(u, equations, dg, i, j, element_id))
    for l in eachnode(dg), k in eachnode(dg)
      acc += get_node_vars(old_u, equations, dg, k, l, old_element_id) * forward_lower[i, k] * forward_lower[j, l]
    end
    set_node_vars!(u, acc, equations, dg, i, j, lower_left_id)
  end

  # Interpolate to lower right element
  for j in eachnode(dg), i in eachnode(dg)
    acc = zero(get_node_vars(u, equations, dg, i, j, element_id))
    for l in eachnode(dg), k in eachnode(dg)
      acc += get_node_vars(old_u, equations, dg, k, l, old_element_id) * forward_upper[i, k] * forward_lower[j, l]
    end
    set_node_vars!(u, acc, equations, dg, i, j, lower_right_id)
  end

  # Interpolate to upper left element
  for j in eachnode(dg), i in eachnode(dg)
    acc = zero(get_node_vars(u, equations, dg, i, j, element_id))
    for l in eachnode(dg), k in eachnode(dg)
      acc += get_node_vars(old_u, equations, dg, k, l, old_element_id) * forward_lower[i, k] * forward_upper[j, l]
    end
    set_node_vars!(u, acc, equations, dg, i, j, upper_left_id)
  end

  # Interpolate to upper right element
  for j in eachnode(dg), i in eachnode(dg)
    acc = zero(get_node_vars(u, equations, dg, i, j, element_id))
    for l in eachnode(dg), k in eachnode(dg)
      acc += get_node_vars(old_u, equations, dg, k, l, old_element_id) * forward_upper[i, k] * forward_upper[j, l]
    end
    set_node_vars!(u, acc, equations, dg, i, j, upper_right_id)
  end

  return nothing
end



# Coarsen elements in the DG solver based on a list of cell_ids that should be removed
function coarsen!(u_ode::AbstractVector, adaptor, mesh::TreeMesh{2}, equations, dg::DGSEM, cache, child_cells_to_coarsen)
  # Return early if there is nothing to do
  if isempty(child_cells_to_coarsen)
    return
  end

  # Determine for each old element whether it needs to be removed
  to_be_removed = falses(nelements(dg, cache))
  # The "Ref(...)" is such that we can vectorize the search but not the array that is searched
  elements_to_remove = searchsortedfirst.(Ref(cache.elements.cell_ids[1:nelements(dg, cache)]),
                                          child_cells_to_coarsen)
  to_be_removed[elements_to_remove] .= true

  # Retain current solution data
  old_n_elements = nelements(dg, cache)
  old_u_ode = copy(u_ode)
  old_u     = wrap_array(old_u_ode, mesh, equations, dg, cache)

  # Get new list of leaf cells
  leaf_cell_ids = leaf_cells(mesh.tree)

  # TODO: Taal debug
  pointer_before_coarsen = pointer(u_ode)
  length_before_coarsen  = length(u_ode)
  let u = wrap_array(u_ode, mesh, equations, dg, cache)
    unstable_idx = findfirst(u -> abs(u) < 0.01, u)
    if unstable_idx !== nothing
      @info "AMR, inside before coarsen!" unstable_idx u[unstable_idx] pointer_before_coarsen length_before_coarsen findfirst(u -> abs(u) < 0.01, old_u)
    end
    # unstable_idx_ode = findfirst(u -> abs(u) < 0.01, u_ode)
    # if unstable_idx !== nothing
    #   @info "AMR, inside before coarsen!" unstable_idx_ode u_ode[unstable_idx_ode]
    # end

    unstable_idx = findfirst(u -> abs(u) < 0.01, u_ode)
    if unstable_idx !== nothing
      @info "AMR, inside before coarsen! for u_ode" unstable_idx u_ode[unstable_idx]
    end

    # unstable_idx = findfirst(u -> abs(u) < 0.01, old_u_ode)
    # if unstable_idx !== nothing
    #   @info "AMR, inside before coarsen! for old_u_ode" unstable_idx old_u_ode[unstable_idx]
    # end

    unstable_idx = findfirst(u -> abs(u) < 0.01, old_u)
    if unstable_idx !== nothing
      @info "AMR, inside before coarsen! for old_u" unstable_idx old_u[unstable_idx]
    end
  end

  # Initialize new elements container
  elements = init_elements(leaf_cell_ids, mesh,
                           real(dg), nvariables(equations), polydeg(dg))
  copy!(cache.elements, elements)
  @assert nelements(dg, cache) == nelements(elements) # TODO: Taal debug
  @assert nelements(dg, cache) < old_n_elements # TODO: Taal debug

  resize!(u_ode, nvariables(equations) * nnodes(dg)^ndims(mesh) * nelements(dg, cache))
  u = wrap_array(u_ode, mesh, equations, dg, cache)

  # TODO: Taal debug
  pointer_medium_coarsen = pointer(u_ode)
  length_medium_coarsen  = length(u_ode)
  let u = wrap_array(u_ode, mesh, equations, dg, cache)
    unstable_idx = findfirst(u -> abs(u) < 0.01, u)
    if unstable_idx !== nothing
      @info "AMR, inside medium coarsen!" unstable_idx u[unstable_idx] pointer_before_coarsen pointer_medium_coarsen length_before_coarsen length_medium_coarsen findfirst(u -> abs(u) < 0.01, old_u)
    end

    unstable_idx = findfirst(u -> abs(u) < 0.01, old_u)
    if unstable_idx !== nothing
      @info "AMR, inside medium coarsen! for old_u" unstable_idx old_u[unstable_idx]
    end
  end

  # Loop over all elements in old container and either copy them or coarsen them
  skip = 0
  element_id = 1
  for old_element_id in 1:old_n_elements
    # If skip is non-zero, we just coarsened 2^ndims elements and need to omit the following elements
    if skip > 0
      skip -= 1
      continue
    end

    if to_be_removed[old_element_id]
      # If an element is to be removed, sanity check if the following elements
      # are also marked - otherwise there would be an error in the way the
      # cells/elements are sorted
      @assert all(to_be_removed[old_element_id:(old_element_id+2^ndims(mesh)-1)]) "bad cell/element order"

      # Coarsen elements and store solution directly in new data structure
      # TODO: Taal debug
      # coarsen_elements!(u, element_id, old_u, old_element_id,
      #                   adaptor, equations, dg)
      # u[:, :, :, element_id] .= 0.5
      for j in eachnode(dg), i in eachnode(dg), v in eachvariable(equations)
        u[v, i, j, element_id] = old_u[v, i, j, old_element_id]
      end
      element_id += 1
      skip = 2^ndims(mesh) - 1
    else
      # Copy old element data to new element container
      # TODO: Taal debug
      # @views u[:, .., element_id] .= old_u[:, .., old_element_id]
      for j in eachnode(dg), i in eachnode(dg), v in eachvariable(equations)
        u[v, i, j, element_id] = old_u[v, i, j, old_element_id]
      end
      element_id += 1
    end
  end
  @assert element_id == nelements(dg, cache) + 1 "element_id = $element_id, nelements(dg, cache) = $(nelements(dg, cache))" # TODO: Taal debug

  # TODO: Taal debug
  pointer_after__coarsen = pointer(u_ode)
  length_after__coarsen  = length(u_ode)
  let u = wrap_array(u_ode, mesh, equations, dg, cache)
    unstable_idx = findfirst(u -> abs(u) < 0.01, u)
    if unstable_idx !== nothing
      @info "AMR, inside after coarsen!" unstable_idx u[unstable_idx] pointer_before_coarsen pointer_medium_coarsen pointer_after__coarsen length_before_coarsen length_medium_coarsen length_after__coarsen findfirst(u -> abs(u) < 0.01, old_u)
    end
  end

  # TODO: Taal performance, allow initializing the stuff in place, making use of resize!
  # Initialize new interfaces container
  interfaces = init_interfaces(leaf_cell_ids, mesh, elements,
                               real(dg), nvariables(equations), polydeg(dg))
  copy!(cache.interfaces, interfaces)

  # Initialize boundaries
  boundaries, _ = init_boundaries(leaf_cell_ids, mesh, elements,
                                  real(dg), nvariables(equations), polydeg(dg))
  copy!(cache.boundaries, boundaries)

  # Initialize new mortar containers
  mortars = init_mortars(leaf_cell_ids, mesh, elements,
                         real(dg), nvariables(equations), polydeg(dg), dg.mortar)
  copy!(cache.mortars, mortars)

  # Sanity check
  if isperiodic(mesh.tree) && nmortars(mortars) == 0
    @assert ninterfaces(interfaces) == 2 * nelements(dg, cache) ("For 2D and periodic domains and conforming elements, the number of interfaces must be twice the number of elements")
  end

  return nothing
end


# TODO: Taal compare performance of different implementations
# Coarsen solution data u for four elements, using L2 projection
function coarsen_elements!(u::AbstractArray{<:Any,4}, element_id, old_u, old_element_id,
                           adaptor::LobattoLegendreAdaptorL2, equations, dg)
  @unpack reverse_upper, reverse_lower = adaptor

  # Store old element ids
  lower_left_id  = old_element_id
  lower_right_id = old_element_id + 1
  upper_left_id  = old_element_id + 2
  upper_right_id = old_element_id + 3

  # TODO: Taal debug
  @assert old_element_id >= 1
  @assert size(old_u, 1) == nvariables(equations)
  @assert size(old_u, 2) == nnodes(dg)
  @assert size(old_u, 3) == nnodes(dg)
  @assert size(old_u, 4) >= old_element_id + 3
  @assert     element_id >= 1
  @assert size(    u, 1) == nvariables(equations)
  @assert size(    u, 2) == nnodes(dg)
  @assert size(    u, 3) == nnodes(dg)
  @assert size(    u, 4) >= element_id

  for j in eachnode(dg), i in eachnode(dg)
    acc = zero(get_node_vars(u, equations, dg, i, j, element_id))

    # Project from lower left element
    for l in eachnode(dg), k in eachnode(dg)
      acc += get_node_vars(old_u, equations, dg, k, l, lower_left_id) * reverse_lower[i, k] * reverse_lower[j, l]
    end

    # Project from lower right element
    for l in eachnode(dg), k in eachnode(dg)
      acc += get_node_vars(old_u, equations, dg, k, l, lower_right_id) * reverse_upper[i, k] * reverse_lower[j, l]
    end

    # Project from upper left element
    for l in eachnode(dg), k in eachnode(dg)
      acc += get_node_vars(old_u, equations, dg, k, l, upper_left_id) * reverse_lower[i, k] * reverse_upper[j, l]
    end

    # Project from upper right element
    for l in eachnode(dg), k in eachnode(dg)
      acc += get_node_vars(old_u, equations, dg, k, l, upper_right_id) * reverse_upper[i, k] * reverse_upper[j, l]
    end

    # Update value
    set_node_vars!(u, acc, equations, dg, i, j, element_id)
  end
end


function indicator_cache(mesh::TreeMesh{2}, equations, dg::DG, cache)

  indicator_value = Vector{real(dg)}(undef, nelements(dg, cache))
  return (; indicator_value)
end

# TODO: Taal refactor, merge the two loops of IndicatorThreeLevel and IndicatorLöhner?
#       But that would remove the simplest possibility to write that stuff to a file...
#       We could of course implement some additional logic and workarounds, but is it worth the effort?
function (indicator::IndicatorThreeLevel)(u::AbstractArray{<:Any,4},
                                          mesh::TreeMesh{2}, equations, dg::DG, cache)

  @unpack indicator_value = indicator.cache
  resize!(indicator_value, nelements(dg, cache))

  alpha = indicator.indicator(u, equations, dg, cache)

  #=Threads.@threads=# for element in eachelement(dg, cache)
    cell_id = cache.elements.cell_ids[element]
    current_level = mesh.tree.levels[cell_id]

    # set target level
    target_level = current_level
    if alpha[element] > indicator.max_threshold
      target_level = indicator.max_level
    elseif alpha[element] > indicator.med_threshold
      if indicator.med_level > 0
        target_level = indicator.med_level
        # otherwise, target_level = current_level
        # set med_level = -1 to implicitly use med_level = current_level
      end
    else
      target_level = indicator.base_level
    end

    # compare target level with actual level to set indicator
    if current_level < target_level
      indicator_value[element] = 1 # refine!
    elseif current_level > target_level
      indicator_value[element] = -1 # coarsen!
    else
      indicator_value[element] = 0 # we're good
    end
  end

  return indicator_value
end


# this method is used when the indicator is constructed as for shock-capturing volume integrals
function create_cache(::Type{IndicatorLöhner}, equations::AbstractEquations{2}, basis::LobattoLegendreBasis)

  alpha = Vector{real(basis)}()

  A = Array{real(basis), ndims(equations)}
  indicator_threaded = [A(undef, nnodes(basis), nnodes(basis)) for _ in 1:Threads.nthreads()]

  return (; alpha, indicator_threaded)
end

# this method is used when the indicator is constructed as for AMR
function create_cache(typ::Type{IndicatorLöhner}, mesh, equations::AbstractEquations{2}, dg::DGSEM, cache)
  create_cache(typ, equations, dg.basis)
end


function (löhner::IndicatorLöhner)(u::AbstractArray{<:Any,4}, equations, dg::DGSEM, cache)
  @assert nnodes(dg) >= 3 "IndicatorLöhner only works for nnodes >= 3 (polydeg > 1)"
  @unpack alpha, indicator_threaded = löhner.cache
  resize!(alpha, nelements(dg, cache))

  #=Threads.@threads=# for element in eachelement(dg, cache)
    indicator = indicator_threaded[Threads.threadid()]

    # Calculate indicator variables at Gauss-Lobatto nodes
    for j in eachnode(dg), i in eachnode(dg)
      u_local = get_node_vars(u, equations, dg, i, j, element)
      indicator[i, j] = löhner.variable(u_local, equations)
    end

    estimate = zero(real(dg))
    for j in eachnode(dg), i in 2:nnodes(dg)-1
      # x direction
      u0 = indicator[i,   j]
      up = indicator[i+1, j]
      um = indicator[i-1, j]
      estimate = max(estimate, löhner(um, u0, up))
    end

    for j in 2:nnodes(dg)-1, i in eachnode(dg)
      # y direction
      u0 = indicator[i, j, ]
      up = indicator[i, j+1]
      um = indicator[i, j-1]
      estimate = max(estimate, löhner(um, u0, up))
    end

    # use the maximum as DG element indicator
    alpha[element] = estimate
  end

  return alpha
end

# dirty Löhner estimate, direction by direction, assuming constant nodes
@inline function (löhner::IndicatorLöhner)(um::Real, u0::Real, up::Real)
  num = abs(up - 2 * u0 + um)
  den = abs(up - u0) + abs(u0-um) + löhner.f_wave * (abs(up) + 2 * abs(u0) + abs(um))
  return num / den
end



# this method is used when the indicator is constructed as for shock-capturing volume integrals
function create_cache(::Type{IndicatorMax}, equations::AbstractEquations{2}, basis::LobattoLegendreBasis)

  alpha = Vector{real(basis)}()

  A = Array{real(basis), ndims(equations)}
  indicator_threaded = [A(undef, nnodes(basis), nnodes(basis)) for _ in 1:Threads.nthreads()]

  return (; alpha, indicator_threaded)
end

# this method is used when the indicator is constructed as for AMR
function create_cache(typ::Type{IndicatorMax}, mesh, equations::AbstractEquations{2}, dg::DGSEM, cache)
  cache = create_cache(typ, equations, dg.basis)
end


function (indicator_max::IndicatorMax)(u::AbstractArray{<:Any,4}, equations, dg::DGSEM, cache)
  @unpack alpha, indicator_threaded = indicator_max.cache
  resize!(alpha, nelements(dg, cache))

  #=Threads.@threads=# for element in eachelement(dg, cache)
    indicator = indicator_threaded[Threads.threadid()]

    # Calculate indicator variables at Gauss-Lobatto nodes
    for j in eachnode(dg), i in eachnode(dg)
      u_local = get_node_vars(u, equations, dg, i, j, element)
      indicator[i, j] = indicator_max.variable(u_local, equations)
    end

    alpha[element] = maximum(indicator)
  end

  return alpha
end
