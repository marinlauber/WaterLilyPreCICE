using PreCICE
using StaticArrays

# can be cleaner
include("MeshBodies.jl")
include("KDTree.jl")

"""
    An structure to hold the coupling data for WaterLily-CalculiX coupling
"""
struct CalculiXInterface <: AbstractInterface
    ControlPointsID::AbstractArray
    ControlPoints::AbstractArray
    # quadPointID::AbstractArray
    # quadPoint::AbstractArray
    forces::AbstractArray
    deformation::AbstractArray
    map_id :: AbstractVector
    dt::Vector{Float64}
end

function CalculiXInterface(T=Float64;fname="geom.inp",map=(x,t)->x,scale=1.f0,kwargs...)  
    # construct the body
    #@TODO this is just to get the connectivity in the end...
    body = MeshBody(fname ; map, scale, T)
    # get nodes and elements IDS from precice
    numberOfVertices, dimensions = length(body.mesh.position), 3
    vertices = Array{Float64,2}(undef, numberOfVertices, dimensions)
    for i in 1:numberOfVertices
        vertices[i,:] .= body.mesh.position[i].data
    end
    ControlPointsID = PreCICE.setMeshVertices("Fluid-Mesh", vertices)
    # quadPoint = Array{Float64,2}(undef, length(body.mesh), dimensions)
    # quadPoint .= hcat(center.(body.mesh)...)'
    # quadPointID = copy(quadPoint)
    # quadPointID = PreCICE.setMeshVertices("Fluid-Mesh", quadPoint)

    # storage arrays
    forces = zeros(typeof(scale), size(vertices))
    deformation = zeros(typeof(scale), size(vertices))
    # @show  map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(body.mesh)))
    # mapping from center to nodes, needed for the forces
    # map_id = map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(body.mesh)))
    map_id = zeros(10)

    # initilise PreCICE
    PreCICE.initialize()
    dt = PreCICE.getMaxTimeStepSize()
    interface = CalculiXInterface(ControlPointsID, vertices, forces, deformation, map_id, [dt])
    # return coupling interface
    return interface, body
end

function readData!(interface::CalculiXInterface)
    # Read control point displacements
    interface.deformation .= PreCICE.readData("Fluid-Mesh", "Displacements", 
                                              interface.ControlPointsID, interface.dt[end])
end

function update!(interface::CalculiXInterface, sim; kwargs...)
    # update mesh position, measure is done elsewhere
    points = Point3f[]
    for (i,pnt) in enumerate(interface.mesh0.position)
        push!(points, Point3f(SA[pnt.data...] .+ interface.deformation[i,:]))
    end
    sim.body.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(interface.mesh0))
end

function get_forces!(interface::CalculiXInterface,flow::Flow,body;Kwargs...)
    t = sum(@views(interface.dt[1:end])) # the time
    interface.forces .= 0 # reset the forces
    # compute nodal forces
    # for (i,id) in interface.srfID
        # f = dS(@views(interface.mesh[id])).*interp()
        # interface.forces[interface.map_id[id],:] .+= transpose(f)./3 # add all the contribution from the faces to the nodes
    # end
end

function writeData!(interface::CalculiXInterface)
    # write the force at the integration points
    #@TODO: permutedims might not be necessary
    PreCICE.writeData("Fluid-Mesh", "Forces", interface.ControlPointsID, interface.forces)
end
