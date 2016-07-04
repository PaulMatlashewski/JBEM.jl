# The mesh data structure is based on the paper:
# Efficient Representation of Computational Meshes. A. Logg [1]
# This describes the general data structure as well as the algorithms to compute
# mesh connectivities

using FixedSizeArrays
import GeometryTypes: Point, Simplex

abstract BoundaryCondition
"""
The ambient dimension of the problem. These types don't hold data and only serve
as a flag for function dispatch.
"""
abstract Dimension
immutable Dimension2 <: Dimension end
immutable Dimension3 <: Dimension end

type Dirichlet{T<:Function} <: BoundaryCondition
    val ::T
end

type Neumann{T<:Function} <: BoundaryCondition
    val ::T
end

type MeshBoundaryCondition{T<:BoundaryCondition}
    bcs     ::Array{T}
    markers ::Array{Int64}
end

"""
A MeshEntity is a unique tuple (d,i), where d is the topological dimension of
the entitiy and i is the index for the mesh entity.
"""
type MeshEntity{N,T} <: FixedVector{N,T}
    _::NTuple{N,T}
end
function call{T<:MeshEntity, F<:FixedVector}(::Type{T}, f::F...)
    MeshEntitiy{length(f),F}(f)
end

"""
The mesh connectivity describes how different mesh entities are connected to
each other. For example, it describes which vertices belong to a mesh face or
which edges belong to a mesh vertex. From [1], all the connectivity combinations
can be derived from knowing only the connectivity between the mesh cells and the
mesh vertices.

The notation for the connectivity of dimension d to dimension d' is d -> d'.

The indices array contains the indices for the dimension d' entities that are
connected to the dimension d entities. The kth element in the offsets array
gives the index into the indices array for the entities connected to the kth
entity of dimension d. See [1] for details.
"""
type MeshConnectivity
    indices ::Array{Int64}
    offsets ::Array{Int64}
end

type MeshTopology{D<:Dimension}
    dim          ::D # Ambient dimension of the model
    connectivity ::Dict{AbstractString, MeshConnectivity}
end

type MeshGeometry
    vertices     ::Array{Point}
    nodes        ::Array{Point}
    normals      ::Array{Point} # One normal for every mesh cell
end

type Mesh{D<:Dimension}
    topology ::MeshTopology{D}
    geometry ::MeshGeometry
end
