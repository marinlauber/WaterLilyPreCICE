using PreCICE
using Interpolations

# this are the loading curves inside CalculiX
A1 = linear_interpolation([0.,10.,112.], [0.,1.,1.])
A2 = linear_interpolation([0.,10.,12.,112.], [0.,1.,1.,10.])
A3 = linear_interpolation([0.,10.,12.,112.], [0.,1.,1.,44.2])

function static_inflation(i,t)
    i==1 && return  0.38*A1(t)
    i==6 && return -0.38*A1(t)
    i in [2,3] && return -0.1*A3(t)
    i in [4,5] && return  0.48*A2(t)
    i in [7,8] && return  0.1*A3(t)
    i in [9,10] && return -0.48*A2(t)
end
function static_pressure!(forces,body,func;t=0)
    for (i,id) in body.srfID
        forces[id,:] .= WaterLilyPreCICE.normal(@views(body.mesh[id])).*func(i,t)
    end
end

struct LumpedInterface <: AbstractInterface
    mesh :: GeometryBasics.Mesh
    srfID :: AbstractVector
    forces :: AbstractArray
    ControlPointsID::AbstractArray
    quadPointID::AbstractArray
    dt :: Vector{Float64}
end
export LumpedInterface

function LumpedInterface(args...;fname="../Solid/geom.inp",T=Float64)  
    
    # load the file
    mesh,srf_id = load_inp(fname)
    
    # get nodes and elements IDS from precice
    numberOfVertices, dimensions = length(mesh.position), 3
    vertices = Array{T,2}(undef, numberOfVertices, dimensions)
    vertices .= hcat(mesh.position...)'
    ControlPointsID = PreCICE.setMeshVertices("Fluid-Mesh-Nodes", vertices)
    quadPoint = Array{T,2}(undef, length(mesh), dimensions)
    quadPoint .= hcat(center.(mesh)...)'
    quadPointID = PreCICE.setMeshVertices("Fluid-Mesh-Faces", quadPoint)

    # storage arrays
    forces = zeros(T, size(quadPoint))
    
    # initilise PreCICE
    PreCICE.initialize()
    dt = PreCICE.getMaxTimeStepSize()
    srf_id = mapreduce(((i,ids),)->map(T->(i,T),ids),vcat,enumerate(body.srfID))
    LumpedInterface(mesh,srf_id,forces,quadPointID,ControlPointsID,[dt])
end

function readData!(interface::LumpedInterface)
    # Read control point displacements
    PreCICE.readData("Fluid-Mesh-Nodes", "Displacements", interface.ControlPointsID, interface.dt[end])
end

function update!(interface::LumpedInterface)
    static_pressure!(interface.forces,interface.body,static_inflation,t=sum(@views(interface.dt)))
end

function writeData!(interface::LumpedInterface)
    # write the force at the quad points
    PreCICE.writeData("Fluid-Mesh-Faces", "Forces", interface.quadPointID, interface.forces)
end