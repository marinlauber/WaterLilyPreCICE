using PreCICE

mutable struct LumpedInterface <: AbstractInterface
    mesh0:: GeometryBasics.Mesh # initial geometry, never changed
    mesh :: GeometryBasics.Mesh
    srfID :: AbstractVector
    vertices :: AbstractArray # might not be needed
    ControlPointsID::AbstractArray
    forces :: AbstractArray
    quadPointID::AbstractArray
    func :: Function
    dt :: Vector{Float64}
end
export LumpedInterface

function LumpedInterface(;fname="../Solid/geom.inp",T=Float64,func=(i,t)->0)

    # keyword aguments might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # create coupling
    PreCICE.createParticipant("LPM", configFileName, 0, 1)
    
    # load the file
    mesh,srf_id = load_inp(fname)
    
    # get nodes and elements IDS from preCICE
    numberOfVertices, dimensions = length(mesh.position), 3
    vertices = Array{T,2}(undef, numberOfVertices, dimensions)
    for i in 1:numberOfVertices
        vertices[i,:] .= mesh.position[i].data
    end
    ControlPointsID = PreCICE.setMeshVertices("LPM-Nodes", vertices)
    quadPoint = Array{T,2}(undef, length(mesh), dimensions)
    for i in 1:numberOfVertices
        quadPoint[i,:] .= center(mesh[i])
    end
    quadPointID = PreCICE.setMeshVertices("LPM-Faces", quadPoint)
    @show ControlPointsID
    @show quadPointID
    # storage arrays
    forces = zeros(T, size(quadPoint))
    
    # initialise PreCICE
    PreCICE.initialize()
    # dt = PreCICE.getMaxTimeStepSize()
    srf_id = mapreduce(((i,ids),)->map(T->(i,T),ids),vcat,enumerate(srf_id))
    LumpedInterface(mesh,deepcopy(mesh),srf_id,vertices,ControlPointsID,forces,quadPointID,func,Int[0.0])
end
function readData!(interface::LumpedInterface)
    # Read control point displacements
    interface.vertices = PreCICE.readData("LPM-Nodes", "Displacements", interface.ControlPointsID, interface.dt[end])
    # update the mesh
    points = Point3f[]
    for (i,pnt) in enumerate(interface.mesh0.position)
        push!(points, Point3f(SA[pnt.data...] .+ interface.vertices[i,:]))
    end
    interface.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(interface.mesh0))
end

function update!(interface::LumpedInterface)
    t = sum(@views(interface.dt[1:end])) # the time
    for (i,id) in interface.srfID
        interface.forces[id,:] .= dS(@views(interface.mesh0[id])).*interface.func(i,t)
    end
end

function writeData!(interface::LumpedInterface)
    # write the force at the quad points
    PreCICE.writeData("LPM-Faces", "Forces", interface.quadPointID, interface.forces)
end
import WaterLily
function WaterLily.write!(w,a::LumpedInterface)
    k = w.count[1]
    points = hcat([[p.data...] for p ∈ a.mesh0.position]...)
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, Base.to_index.(face)) for face in faces(a.mesh0)]
    vtk = vtk_grid(w.dir_name*@sprintf("/%s_%06i", w.fname, k), points, cells) 
    for (name,func) in w.output_attrib
        vtk[name] = func(a)
    end
    vtk_save(vtk); w.count[1]=k+1
    w.collection[round(sum(a.dt),digits=4)]=vtk
end
