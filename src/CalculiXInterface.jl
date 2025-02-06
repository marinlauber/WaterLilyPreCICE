using PreCICE
using StaticArrays

# can be cleaner
include("MeshBodies.jl")
include("KDTree.jl")

"""
    An structure to hold the coupling data for WaterLily-CalculiX coupling
"""
struct CalculiXInterface <: AbstractInterface
    U::Float64
    ControlPointsID::AbstractArray
    ControlPoints::AbstractArray
    quadPointID::AbstractArray
    quadPoint::AbstractArray
    forces::AbstractArray
    deformation::AbstractArray
    dt::Vector{Float64}
end

function CalculiXInterface(U,L;fname="geom.inp",map=(x,t)->x,scale=1.00,T=Float64)  
    
    # construct the body
    println("REading mesh file ", fname)
    println("Using map(0.) -> ", map(SA[0.,0.],0.0))
    println("Using scale ", scale)
    println("Using T ", T)
    body = MeshBody(fname;map,scale,T)

    # get nodes and elements IDS from precice
    #@TODO: position and quadPoitns might need correct reshaping
    numberOfVertices, dimensions = length(body.mesh.position), 3
    vertices = Array{Float64,2}(undef, numberOfVertices, dimensions)
    vertices .= hcat(body.mesh.position...)'
    ControlPointsID = PreCICE.setMeshVertices("Fluid-Mesh-Nodes", vertices)
    quadPoint = Array{Float64,2}(undef, length(body.mesh), dimensions)
    quadPoint .= hcat(center.(body.mesh)...)'
    quadPointID = PreCICE.setMeshVertices("Fluid-Mesh-Faces", quadPoint)

    # storage arrays
    forces = zeros(typeof(scale), size(quadPoint))
    deformation = zeros(typeof(scale), size(vertices))

    # initilise PreCICE
    PreCICE.initialize()
    dt = PreCICE.getMaxTimeStepSize()
    interface = CalculiXInterface(U, ControlPointsID, vertices, quadPointID, quadPoint, forces, deformation, [dt])
    
    # return coupling interface
    return interface, body
end

function readData!(interface::CalculiXInterface)
    # Read control point displacements
    interface.deformation .= PreCICE.readData("Fluid-Mesh-Nodes", "Displacements", 
                                              interface.ControlPointsID, interface.dt[end])
end


function update!(interface::CalculiXInterface,sim::Simulation;kwargs...)
    # update the geom as this has not been done yet
    # velocity = (position .- sim.body.position)./sim.flow.Δt[end]
    # points = interface.ControlPoints .+ interface.deformation
    # update the body
    # sim.body = GeometryBasics.Mesh(points,GeometryBasics.faces(sim.body.mesh))
end

function getInterfaceForces!(interface::CalculiXInterface,flow::Flow,body)
    # get the forces at the integration points
    interface.forces .= 0.0
end

function writeData!(interface::CalculiXInterface,sim::Simulation,store::Store)
    # write the force at the integration points
    #@TODO: permutedims might not be necessary
    PreCICE.writeData("Fluid-Mesh-Faces", "Forces", interface.quadPointID, interface.forces)
end
