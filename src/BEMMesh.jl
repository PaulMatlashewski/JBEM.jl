include("BEMModel.jl")
using Base.Test

################################################################################
# Connectivity dictionary key functions
#
"""
    inverse_key(key::AbstractArray)

Given the key for the mesh connectivity d -> d', return the key for the
connectivity d' -> d
"""
function inverse_key(key::AbstractString)
    if key == "V2V"
        return "V2V"
    elseif key == "V2E"
        return "E2V"
    elseif key == "V2F"
        return "F2V"
    elseif key == "E2V"
        return "V2E"
    elseif key == "E2E"
        return "E2E"
    elseif key == "E2F"
        return "F2E"
    elseif key == "F2V"
        return "V2F"
    elseif key == "F2E"
        return "E2F"
    elseif key == "F2F"
        return "F2F"
    elseif key == "E2N"
        return "N2E"
    elseif key == "N2E"
        return "E2N"
    else
        throw(KeyError(key))
    end
end
################################################################################

################################################################################
# Connectivity access functions
#
"""
    get_indices(connectivity::MeshConnectivity, i::Int64)

Get the indices of the desired mesh entity. The argument `connectivity`
specifies a mesh connectivity d -> d'. The argument `i` specifies the mesh
entitiy (d, i). This function returns an array of `indices` that specifies
which entities are connected to (d, i). That is, it allows one to determine the
set of entities connected to (d, i) defined by { (d',j) : j ∈ `indices` }.
"""
function get_indices(connectivity::MeshConnectivity, i::Int64)
    start_index = connectivity.offsets[i]
    end_index = connectivity.offsets[i + 1] - 1
    indices = connectivity.indices[start_index:end_index]
    return indices
end

"""
    get_vertices(connectivity::MeshConnectivity, geometry::MeshGeometry,
                 i::Int64)

Get the vertices connected to the desired mesh entity.
"""
function get_vertices(connectivity::MeshConnectivity, geometry::MeshGeometry,
                      i::Int64)
    indices = get_indices(connectivity, i)
    vertices = geometry.vertices[indices]
    return vertices
end

"""
    update_connectivity!(connectivity::MeshConnectivity, index::Int64,
                         rep_number::Int64)

Update the mesh connectivity fields between specified mesh entities. If the mesh
connectivity we wish to update is d -> d', the `index` argument specifies the
mesh entitiy (d, index) while the argument `rep_number` specifies the mesh
entitiy (d', rep_number).
"""
function update_connectivity!(connectivity::MeshConnectivity, index::Int64,
                              rep_number::Int64)
    # The index in the indices array to add the replacement number
    rep_index = connectivity.offsets[index]
    if rep_index > length(connectivity.indices)
        # Add the replacement number to the end
        push!(connectivity.indices, rep_number)
    else
        # Insert the replacement number at replacement index. This looks strange
        # but this is how it's done in Julia...
        splice!(connectivity.indices, rep_index, [rep_number, connectivity.indices[rep_index]])
    end
    # Offset all the entity offsets after index to take into account the fact
    # we just added an index.
    for i in index + 1:length(connectivity.offsets)
        connectivity.offsets[i] += 1
    end
end

################################################################################
# Computational mesh functions and algorithms based on A. Logg paper
#
"""
    connectivity_transpose(connectivity::MeshConnectivity)

Given a mesh connectivity d -> d', compute the connectivity d' -> d for d' < d.
"""
function connectivity_transpose!(topology::MeshTopology, key::AbstractString)
    # Get the transpose key
    transpose_key = inverse_key(key)

    # Starting Connectivity
    connectivity = topology.connectivity[key]

    # Number of dimension d' entities
    ndprime = maximum(connectivity.indices)
    # Number of dimension d entities
    nd = length(connectivity.offsets) - 1

    # Initialize connectivity indices
    new_indices = Int64[] # A priori unknown index array length

    # Initialize the offset values with 1. This will be used to place the
    # cell indices in the correct position of the array. The algorithm
    # increments the appropriate offset elements after placing an index
    new_offsets = ones(Int64, ndprime + 1) # Known number of d' entities

    topology.connectivity[transpose_key] = MeshConnectivity(new_indices,
                                                            new_offsets)
    # The idea here is to use the indices of d' connected to the jth
    # d-entity as an index into the current offset array. This value in the
    # current offset array gives the location to place the d-entity number in
    # the indices array. After placing the d-entity index in the indices array,
    # increment the elements in the offsets array at positions greater than the
    # d'-entity number by 1. After looping through all the cells, this will give
    # the correct connectivity information.
    for j in 1:nd
        # Get the d' indices for the jth cell
        indices = get_indices(connectivity, j)
        for index in indices
            update_connectivity!(topology.connectivity[transpose_key], index, j)
        end
    end
end

"""
    connectivity_build!(topology::MeshTopology)
"""
# TO DO: Implament this from the A. Logg. Paper. This may me needed for 3D
# Mesh computation
function connectivity_build!(topology::MeshTopology)
end

"""
    connectivity_intersection(topology::MeshTopology, key_start::AbstractString,
                              key_end::AbstractString)

Compute the connectivity d -> d' from d -> d" and d" -> d' for d >= d'.
key_start is the key for d -> d" and key_end is the key for d" -> d'.
"""
# TO DO: Implament this from the A. Logg. Paper. This may me needed for 3D
# Mesh computation
function connectivity_intersection!(topology::MeshTopology,
                                    key_start::AbstractString,
                                    key_end::AbstractString)
end

"""
    generate_connectivity!(mesh::Mesh, key::AbstractString)

Build the mesh connectivity corresponding to the key parameter.
"""
# TO DO: Might have to make this general based on the A. Logg. algorithm.
# For now, only special cases for key are implemented. Currently only
# to find E2V for 2D meshes.
function generate_connectivity!(mesh::Mesh{Dimension2}, key::AbstractString)
end
################################################################################

################################################################################
# Miscellaneous mesh functions
#
"""
    generate_normals!(mesh::Mesh)

Compute and store the unit normal vectors of the mesh. Each mesh cell has an
associated unit normal. This normal is used during integral evaluations of the
normal derivative of the free space Green's function as well as in determining
if discontinuous elements are required for flux on the boundary.

Currently only implemented for 2D linear elements.
"""
# TO DO: It would be nice if we didn't have to count on the
# vertices being given in order. This is minor, so fix up other things first.
function generate_normals!(mesh::Mesh{Dimension2})
    # Number of mesh cells
    nd = length(mesh.topology.connectivity["E2V"].offsets) - 1

    # Treat each mesh segment as a directed line segment. So each segment
    # can be written as [dx, dy] = [pt2.x - pt1.x, pt2.y - pt1.y]. The
    # outward pointing normal is given by [dy, -dx]. This requires that the
    # vertices of the 2D mesh are given in counterclockwise order.
    for i in 1:nd
        vertices = get_vertices(mesh.topology.connectivity["E2V"],
                                mesh.geometry, i)
        dx = vertices[2][1] - vertices[1][1]
        dy = vertices[2][2] - vertices[1][2]
        push!(mesh.geometry.normals, Point(dy, -1.0*dx))
    end
end

# TO DO: 3D version. Some linear algebra stuff to look up (cross products, etc.)
function generate_normals!(mesh::Mesh{Dimension3})
end
################################################################################

################################################################################
# Polynomial interpolation node functions
#
"""
    node_stepping(vertices::Array{Point}, include_verts::Array{Bool},
                  order:Int64)

Compute shape function interpolation nodes on boundary edge.
# Arguments
* `vertices`: the vertices defining the endpoints of the edge
* `include_verts`: an array of booleans indicating if a vertex
should be included as a node. include_verts[i] = true means that vertices[i]
will be included as a vertex
* `order`: the order of the polynomial on the edge. The convention is
different from mathematical convention: 1 = constant, 2 = linear, 3 = quadratic,
etc. This means the number of nodes is equal to order.
"""
function node_stepping(vertices::Array{Point}, include_verts::Array{Bool},
                       order::Int64)
    @test length(vertices) == 2 # The number of vertices for an edge is always 2
    @test length(include_verts) == 2 # Must match the number of vertices.

    nodes = Array{Point}(order)
    delta = Point(vertices[2] - vertices[1])/(order - 1) # step size
    for i in 1:order
        nodes[i] = vertices[1] + delta*(i - 1)
    end

    if include_verts[1] == false # remove the first node
        shift!(nodes)
    end
    if include_verts[2] == false # remove the last node
        pop!(nodes)
    end
    return nodes
end

"""
    function make_edge_nodes!(mesh::Mesh, index::Int64, order::Int64)

Build the node geometry and E2N connectivity for edge (2, `index`). The number
of nodes added is determined by the order of the polynomial on the edge as
indicated by `order`.
"""
##### TO DO: Can this be tidied up. A lots going on here with some repetition.
#####        Can there be some more functions?
function make_edge_nodes!(mesh::Mesh, index::Int64, order::Int64)
    # Get the vertices of the current edge
    vertices = get_vertices(mesh.topology.connectivity["E2V"], mesh.geometry,
                            index)
    # Get the indices of these vertices
    vertex_indices = get_indices(mesh.topology.connectivity["E2V"], index)
    @test length(vertices) == 2 # Edges should always have 2 vertices
    @test length(vertex_indices) == 2

    if order == 1 # Special discontinuous case for constant elements
        push!(mesh.geometry.nodes, mean(vertices)) # Node in the centre of edge
        global_number = length(mesh.geometry.nodes) # Total number of mesh nodes
        # Edge 'index' is connected to the node at index 'global_nodes'
        push!(mesh.topology.connectivity["E2N"].indices, global_number)
        # Increment all other edge offset indices by 1
        mesh.topology.connectivity["E2N"].offsets[index + 1:end] += 1
    else
        # Current number of nodes on edge
        local_nodes = mesh.topology.connectivity["E2N"].offsets[index + 1] -
                      mesh.topology.connectivity["E2N"].offsets[index]
        if local_nodes == 0
            # Make all nodes including nodes coinciding with vertices
            include_verts = [true, true] # Include both vertices as nodes
            nodes = node_stepping(vertices, include_verts, order)
            for i in 1:length(nodes)
                # Add nodes for the current edge (given by index)
                push!(mesh.geometry.nodes, nodes[i]) # Add node to geometry
                # Index of recently added node
                global_number = length(mesh.geometry.nodes)
                update_connectivity!(mesh.topology.connectivity["E2N"], index,
                                     global_number)
                # Since local_nodes == 0, the edges connected to vertices[1] and
                # vertices[2] share nodes[1] and nodes[end] respectively. So
                # we add more copies of global_number to the indices array
                # for these connected edges and modify offsets appropriately
                if i == 1
                    edge_indices = get_indices(mesh.topology.connectivity["V2E"],
                                               vertex_indices[1])
                    for edge in edge_indices
                        if edge == index
                            # This is the current edge we're operating on
                            # so we skip over this. We're only concerned with
                            # the connected edges here.
                        else
                            update_connectivity!(mesh.topology.connectivity["E2N"],
                                                 edge, global_number)
                        end
                    end
                end
                if i == length(nodes)
                    edge_indices = get_indices(mesh.topology.connectivity["V2E"],
                                               vertex_indices[end])
                    for edge in edge_indices
                        if edge == index
                            # This is the current edge we're operating on
                            # so we skip over this. We're only concerned with
                            # the connected edges here.
                        else
                            update_connectivity!(mesh.topology.connectivity["E2N"],
                                                 edge, global_number)
                        end
                    end
                end
            end
        elseif local_nodes == 1 # One of the edge vertices already has a node.
            # To determine this vertex, we first find the edges connected to
            # each vertex.
            i = 1
            include_verts = [true, true] # initialize flags
            for vertex_index in vertex_indices
                edge_indices = get_indices(mesh.topology.connectivity["V2E"],
                                           vertex_index)
                # Since local_nodes == 1, at least one of the edges connected
                # to one of the vertices must already have nodes. At the same
                # time, it MUST be the case that the other vertex CANNOT be
                # connected to edges with nodes (otherwise local_nodes == 2).
                # So the strategy here is to inspect the edges connected to each
                # vertex to find which of the vertices is connected to an edge
                # with nodes.
                for edge in edge_indices
                    if edge == index
                        # This is the current edge we're operating on
                        # so we skip over this. We're only concerned with
                        # the connected edges here.
                    else
                        edge_nodes = mesh.topology.connectivity["E2N"].offsets[edge + 1] -
                                     mesh.topology.connectivity["E2N"].offsets[edge]
                        if edge_nodes == order # This means the edge has nodes
                            include_verts[i] = false
                        end
                    end
                end
                i += 1
            end
            # As stated earlier, only ONE of the vertices can have a connected
            # edge with nodes. We test that here.
            @test (include_verts[1] == true) $ (include_verts[2] == true)

            nodes = node_stepping(vertices, include_verts, order)
            for i in 1:length(nodes)
                # Add nodes for the current edge
                push!(mesh.geometry.nodes, nodes[i]) # Add node to geometry
                # Index of recently added node
                global_number = length(mesh.geometry.nodes)
                update_connectivity!(mesh.topology.connectivity["E2N"], index,
                                     global_number)

                # Update the connected edges to the vertex without the initial node
                if (i == 1) && (include_verts[1] == true)
                    edge_indices = get_indices(mesh.topology.connectivity["V2E"],
                                               vertex_indices[1])
                    for edge in edge_indices
                        if edge == index
                            # This is the current edge we're operating on
                            # so we skip over this. We're only concerned with
                            # the connected edges here.
                        else
                            update_connectivity!(mesh.topology.connectivity["E2N"],
                                                 edge, global_number)
                        end
                    end
                end
                if (i == length(nodes)) && (include_verts[2] == true)
                    edge_indices = get_indices(mesh.topology.connectivity["V2E"],
                                               vertex_indices[end])
                    for edge in edge_indices
                        if edge == index
                        else
                            update_connectivity!(mesh.topology.connectivity["E2N"],
                                                 edge, global_number)

                        end
                    end
                end
            end
        elseif local_nodes == 2
            # Only make the internal nodes. No need to worry about updating
            # connected edges in this case
            include_verts = [false, false]
            nodes = node_stepping(vertices, include_verts, order)
            for i in 1:length(nodes)
                push!(mesh.geometry.nodes, nodes[i])
                global_number = length(mesh.geometry.nodes)
                update_connectivity!(mesh.topology.connectivity["E2N"], index,
                                     global_number)
            end
        else
            throw(ErrorException("Too many nodes on edge"))
        end
    end
end

"""
    generate_nodes!(mesh::Mesh, order::Int64)

Add the polynomial interpolation nodes to the mesh. This function is called once
during the matrix assembly procedure to calculate the geometry and topology of
the nodes with respect to the input mesh.
"""
function generate_nodes!(mesh::Mesh, order::Int64)
    # For this BEM implementation, the interpolation nodes will always be on
    # the mesh edges, for both 2D and 3D meshes.

    # Check if the edge topology is known:
    if haskey(mesh.topology.connectivity, "E2V") == false # Edges unknown
        # Generate the edge -> vertex connectivity. Note that we start with
        # E2V for 2d meshes (since the edge is of maximal dimension D for
        # 2d meshes). So this will only be called for 3d meshes.
        # TO DO: generate_connectivity is not implemented yet.
        generate_connectivity!(mesh, "E2V") # Does not exist yet
    end
    # We need to know which edges are connected to a given vertex. So we need
    # to find the vertex -> edge connectivity. This is simply the transpose of
    # E2V
    if haskey(mesh.topology.connectivity, "V2E") == false
        connectivity_transpose!(mesh.topology, "E2V")
    end

    nedges = length(mesh.topology.connectivity["E2V"].offsets) - 1 # Number of edges

    # The nodes have both a topological and geometric aspect. The geometric
    # aspect is simply their cartesian coordinates, similar to the mesh
    # vertices. Their topological aspect is how they are connected to the mesh.
    # In particular, we are interested in how they are connected to the edges
    # and cells of the mesh.

    # Since the number of nodes for edge k = offsets[k+1] - offsets[k], we give
    # a starting representation of zero nodes for each edge
    offsets = ones(Int64, nedges + 1) # Offsets for edge to node connectivity
    indices = Int64[] # Indices for edge to node connectivity
    mesh.topology.connectivity["E2N"] = MeshConnectivity(indices, offsets)

    # Each edge must have a number of nodes equal to order. Note that edges
    # will share nodes whenever they meet at a vertex.
    for i in 1:nedges
        # Build nodes for each edge
        make_edge_nodes!(mesh, i, order)
    end
    ###### TO DO: This only works for dimension 2!!!! Dimension 3 requires a
    ######        general connectivity generator. We can use transpose here
    ######        because the edges are the cells in dimension 2. So to get
    ######        the node to cell connectivity, we simply take the transpose of
    ######        the E2N
    # Make the node -> cell connectivity

    # TO DO: calculating connectivity transpose after making the nodes takes
    # some time. There's a potential performance increase by calculating the
    # node to cell connectivity while the nodes are being generated in
    # make_edge_nodes!. Performance is already fairly good, so this is a minor
    # issue.
    connectivity_transpose!(mesh.topology, "E2N") # Only for dim 2
end
################################################################################
"""
    discretize_mesh!(mesh::Mesh, order::Int64)

Compute the Lagrange interpolation node points (and the associated connectivity)
and the mesh cell outward unit normals
"""
function discretize_mesh!(mesh::Mesh, order::Int64)
    # Assemble global node system. Skip if the nodes are already in place.
    if length(mesh.geometry.nodes) == 0
        generate_nodes!(mesh, order)
    end
    # Compute mesh normals
    if length(mesh.geometry.normals) == 0
        generate_normals!(mesh)
    end
end

"""
    initialize_matrix(mesh::Mesh, mesh_bcs::MeshBoundaryCondition)

Initialize the system matrix. Returns an NxN matrix with dummy values where N
equals the number of nodes plus the number of split collocation points.
"""
function initialize_matrix(mesh::Mesh{Dimension2},
                           mesh_bcs::MeshBoundaryCondition)
    # The minimum size of the matrix is equal to the number of nodes in the mesh
    N = length(mesh.geometry.nodes)

    # We add a collocation equation whenever there are two degrees of freedom at
    # a mesh node. This happens at corner vertices whose connected edges both
    # have a Dirichlet boundary condition
    index = 1
    for vertex in mesh.geometry.vertices
        edge_indices = get_indices(mesh.topology.connectivity["V2E"], index)
        @test length(edge_indices) == 2 # Must be true for 2D mesh
        if isa(mesh_bcs.bcs[mesh_bcs.markers[edge_indices[1]]], Dirichlet) &
           isa(mesh_bcs.bcs[mesh_bcs.markers[edge_indices[2]]], Dirichlet)
           # Both boundary connected edges have a Dirichlet condition
           N += 1
       end
       index += 1
   end
   A = Array{Float64}(N,N) # Initialize the matrix
   return A
end

"""
    collocate!(A::Array{Float64}, index::Int64, coll_point::Point,
               geometry::MeshGeometry, bcs::MeshBoundaryCondition)

Modify row `index` of matrix `A` with a collocation equation corresponding to
the impulse source located at `coll_point`.
"""
function collocate!(A::Array{Float64}, index::Int64, coll_point::Point,
                    geometry::MeshGeometry, bcs::MeshBoundaryCondition)
end


"""
    assemble_system(mesh::Mesh, mesh_bcs::MeshBoundaryCondition)

Assemble the BEM system Ax = b. This is based on the direct BEM where the system
of equations is generated by selecting collocation points at every element
interpolation node. Extra degrees of freedom are treated with the multipoint
method or with discontinuous elements by moving collocation points inside the
mesh cells whenever required.
"""
function assemble_system(mesh::Mesh{Dimension2},
                         mesh_bcs::MeshBoundaryCondition, order::Int64)
    # Discretize the mesh using Lagrangian elements specified by order. For now,
    # we only use geometrically linear boundary elements although our function
    # can be any order. Future implementation will allow for geometrically
    # curved boundary elements with B-splines (up to cubic). With a B-spline
    # boundary, we can also consider using a functional basis that matches the
    # B-spline basis. Although I don't see why the Lagrangian basis wouldn't
    # work on B-splines as well. This should be tested once both are
    # implemented.
    discretize_mesh!(mesh, order)

    A = initialize_matrix(mesh, mesh_bcs)

    # Use the collocation method to generate the system of equations. The
    # collocation points are chosen to be at the nodes of the mesh. We
    # introduce mulitnodes wherever there is a geometric discontinuity (corner).
    # If the mulitnode method results in redundant equations (as in the case of
    # Dirichlet boundary conditions specified on adjacent cells), the
    # collocation points are shifted inside each boundary cell to ensure
    # unique equations.
    index = 1
    for node in mesh.geometry.nodes
        # Get the edges connected to the current node.
        edge_indices = get_indices(mesh.topology.connectivity["N2E"], index)
        # Check if there are multiple Dirichlet conditions on the connected
        # edges. If so, we need to split and move the collocation point.
        if length(edge_indices) == 1 # This is an internal node
            collocate!(A, index, node, mesh.geometry, mesh_bcs)
        else
            # If there is not one edge attached to a node, there must be two
            # when working with a 2d mesh.
            @test length(edge_indices) == 2
            for edge in edge_indices
                # For the current node, find the nearest adjacent node on the
                # current edge
                node_indices = get_indices(mesh.topology.connectivity["E2N"], edge)
                for i in 1:length(node_indices)

                    if node_indices[i] == node
                        # Skip the current node. We want the closest node
                    else
                        if i == 1
                            # Take this as the minimum distance automatically
                            delta = mesh.geometry.nodes[i] - node
                        else
                            # The minimum distance, delta, must have been
                            # initialized by this point
                            @test isdefined(:delta)
                            test_delta = mesh.geometry.nodes[i] - node
                            if norm(test_delta) < norm(delta)
                                delta = test_delta
                            end
                        end
                    end
                end
                # Place a collocation point some fractional distance d between
                # the current node and the nearest internal node on the current
                # edge.
                d = 0.25 # Test with 0.25 for now
                coll_point = node + d*delta
                collocate!(A, index, coll_point, mesh.geometry, mesh_bcs)
            end
        end
        index += 1
    end
    return A
end
