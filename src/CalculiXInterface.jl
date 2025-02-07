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
    map_id :: AbstractVector
    dt::Vector{Float64}
end

function CalculiXInterface(U,L;fname="geom.inp",map=(x,t)->x,scale=1.00,T=Float64)  
    
    # construct the body
    println("Reading mesh file ", fname)
    println("Using map(0.) -> ", map(SA[0.,0.],0.0))
    println("Using scale ", scale)
    println("Using T ", T)
    body = MeshBody(fname;map,scale,T)

    # get nodes and elements IDS from precice
    #@TODO: position and quadPoitns might need correct reshaping
    numberOfVertices, dimensions = length(body.mesh.position), 3
    vertices = Array{Float64,2}(undef, numberOfVertices, dimensions)
    vertices .= hcat(body.mesh.position...)'
    ControlPointsID = PreCICE.setMeshVertices("Fluid-Mesh", vertices)
    quadPoint = Array{Float64,2}(undef, length(body.mesh), dimensions)
    quadPoint .= hcat(center.(body.mesh)...)'
    # quadPointID = PreCICE.setMeshVertices("Fluid-Mesh", quadPoint)

    # storage arrays
    forces = zeros(typeof(scale), size(vertices))
    deformation = zeros(typeof(scale), size(vertices))

    # mapping from center to nodes, needed for the forces
    map_id = map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(body.mesh)))


    # initilise PreCICE
    PreCICE.initialize()
    dt = PreCICE.getMaxTimeStepSize()
    interface = CalculiXInterface(U, ControlPointsID, vertices,
                                  quadPointID, quadPoint, forces, deformation, map_id, [dt])
    
    # return coupling interface
    return interface, body
end

function readData!(interface::CalculiXInterface)
    # Read control point displacements
    interface.deformation .= PreCICE.readData("Fluid-Mesh", "Displacements", 
                                              interface.ControlPointsID, interface.dt[end])
end


function update!(interface::CalculiXInterface, sim)
    t = sum(@views(interface.dt[1:end])) # the time
    interface.forces .= 0 # reset the forces
    # compute nodal forces
    points = Point3f[]
    for (i,pnt) in enumerate(interface.mesh0.position)
        push!(points, Point3f(SA[pnt.data...] .+ interface.deformation[i,:]))
    end
    sim.body.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(interface.mesh0))
end

function getInterfaceForces!(interface::CalculiXInterface,flow::Flow,body)
    # get the forces at the integration points
    interface.forces .= 0.0
end

function writeData!(interface::CalculiXInterface,sim::Simulation,store::Store)
    # write the force at the integration points
    #@TODO: permutedims might not be necessary
    PreCICE.writeData("Fluid-Mesh", "Forces", interface.ControlPointsID, interface.forces)
end
