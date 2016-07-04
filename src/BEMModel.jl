include("BEMTypes.jl")

function assemble_square_model()
    # Assemble mesh
    dim = Dimension2()

    vertices = Point[]
    push!(vertices, Point(0.0, 0.0))
    push!(vertices, Point(1.0, 0.0))
    push!(vertices, Point(1.0, 1.0))
    push!(vertices, Point(0.0, 1.0))

    # Initialize empty nodes array. Nodes are added during discretization.
    nodes = Point[]
    # Initialize empty normal array. Normal vectors are added during
    # discretization.
    normals = Point[]
    geometry = MeshGeometry(vertices, nodes, normals)

    # Mesh connectivity dictionary
    connectivity = Dict{AbstractString, MeshConnectivity}()
    # Edge to vertex connectivity
    indices = [1, 2, 2, 3, 3, 4, 4, 1]
    offsets = [1, 3, 5, 7, 9]
    connectivity["E2V"] = MeshConnectivity(indices, offsets)
    topology = MeshTopology(dim, connectivity)

    mesh = Mesh(topology, geometry)

    # Assemble boundary condition
    ice(x,y) = 0.0
    fire(x,y) = 100.0
    insulated(x,y) = 0.0
    ice_boundary = Dirichlet(ice)
    fire_boundary = Dirichlet(fire)
    insulated_boundary = Neumann(insulated)
    bcs = BoundaryCondition[ice_boundary, fire_boundary, insulated_boundary]
    markers = [1, 3, 2, 3]
    mesh_bc = MeshBoundaryCondition(bcs, markers)

    return mesh, mesh_bc
end

function assemble_cube_model()
    # Assemble mesh
    dim = Dimension3()

    vertices = Point[]
    push!(vertices, Point(0.0, 0.0, 0.0))
    push!(vertices, Point(1.0, 0.0, 0.0))
    push!(vertices, Point(1.0, 0.0, 1.0))
    push!(vertices, Point(0.0, 0.0, 1.0))
    push!(vertices, Point(1.0, 1.0, 0.0))
    push!(vertices, Point(1.0, 1.0, 1.0))
    push!(vertices, Point(0.0, 1.0, 0.0))
    push!(vertices, Point(0.0, 1.0, 0.0))

    # Initialize empty nodes array. Nodes are added during discretization.
    nodes = Point[]
    # Initialize empty normal array. Normal vectors are added during
    # discretization.
    normals = Point[]
    geometry = MeshGeometry(vertices, nodes, normals)

    # Mesh connectivity dictionary
    connectivity = Dict{AbstractString, MeshConnectivity}()
    # Face to vertex connectivity
    indices = [1, 3, 4, 1, 2, 3, 2, 5, 3, 5, 6, 3, 5, 7, 6, 8, 7, 5, 8, 4, 7, 1,
               8, 4, 4, 3, 7, 3, 6, 7, 1, 2, 5, 1, 5, 8]
    offsets = [1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37]
    F2V = MeshConnectivity(indices, offsets)
    connectivity["F2V"] = F2V
    topology = MeshTopology(dim, connectivity)

    mesh = Mesh(topology, geometry)

    # Assemble boundary condition
    ice(x,y,z) = 0.0
    fire(x,y,z) = 100.0
    insulated(x,y,z) = 0.0
    ice_boundary = Dirichlet(ice)
    fire_boundary = Dirichlet(fire)
    insulated_boundary = Neumann(insulated)
    bcs = BoundaryCondition[ice_boundary, fire_boundary, insulated_boundary]
    markers = [3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1, 1]
    mesh_bc = MeshBoundaryCondition(bcs, markers)

    return mesh, mesh_bc
end
