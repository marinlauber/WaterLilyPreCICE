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
    forces::AbstractArray
    deformation::AbstractArray
    map_id :: AbstractVector
    dt::Vector{Float64}
end

function CalculiXInterface(T=Float64; fname="geom.inp", map=(x,t)->x, scale=1.f0, kwargs...)  
    # construct the body
    #@TODO this is just to get the connectivity in the end...
    body = MeshBody(fname ; map, scale, T)

    # get nodes and elements IDS from precice
    numberOfVertices, dimensions = length(body.mesh.position), 3
    vertices = Array{Float64,2}(undef, numberOfVertices, dimensions)
    for i in 1:numberOfVertices
        vertices[i,:] .= body.mesh.position[i].data./scale # we need to have the same scale as in CalculiX
    end
    ControlPointsID = PreCICE.setMeshVertices("Fluid-Mesh", vertices)
   
    # storage arrays
    forces = zeros(Float64, size(vertices))
    deformation = zeros(Float64, size(vertices))

    # mapping from center to nodes, needed for the forces
    map_id = Base.map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(body.mesh)))
    
    # initilise PreCICE
    PreCICE.initialize()
    dt = PreCICE.getMaxTimeStepSize()
    
    # return coupling interface
    interface = CalculiXInterface(ControlPointsID, forces, deformation, map_id, [dt])
    return interface, body
end

function readData!(interface::CalculiXInterface)
    # Read control point displacements
    interface.deformation .= PreCICE.readData("Fluid-Mesh", "Displacements", 
                                              interface.ControlPointsID, interface.dt[end])
end

function update!(interface::CalculiXInterface, sim::CoupledSimulation; kwargs...)
    # update mesh position, measure is done elsewhere
    points = Point3f[]
    for (i,pnt) in enumerate(sim.body.mesh0.position)
        push!(points, Point3f(SA[pnt.data...] .+ sim.body.scale.*interface.deformation[i,:]))
    end
    # update
    sim.body.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(sim.body.mesh0))
end

import WaterLily: interp
function get_forces!(interface::CalculiXInterface, flow::Flow, body::MeshBody; δ=1.f0, kwargs...)
    t = sum(@views(interface.dt[1:end])) # the time
    interface.forces .= 0 # reset the forces
    # compute nodal forces
    for (i,id) in body.srfID
        tri = body.mesh[id]
        # map into correct part of the mesh, time does nothing
        f = dS(tri).*interp(body.map(center(tri).+1.5 .+ δ.*normal(tri),t),flow.p)
        interface.forces[interface.map_id[id],:] .+= transpose(f)./3 # add all the contribution from the faces to the nodes
    end
end

function writeData!(interface::CalculiXInterface)
    # write the force at the integration points
    #@TODO: permutedims might not be necessary
    PreCICE.writeData("Fluid-Mesh", "Forces", interface.ControlPointsID, interface.forces)
end
