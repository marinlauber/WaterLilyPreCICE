using PreCICE

mutable struct LumpedInterface <: AbstractInterface
    mesh0:: GeometryBasics.Mesh # initial geometry, never changed
    mesh :: GeometryBasics.Mesh
    srfID :: AbstractVector
    vertices :: AbstractArray # might not be needed
    ControlPointsID::AbstractArray
    forces :: AbstractArray
    σ :: AbstractArray # face storage vector
    func :: Function
    map_id :: AbstractVector
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
    mesh,srf_id = load_inp(fname) # can we get rid of this?
        
    # initialise PreCICE
    PreCICE.initialize()
    
    # we need to initialize before we can get the mesh points and coordinates
    (ControlPointsID, ControlPoints) = PreCICE.getMeshVertexIDsAndCoordinates("Solid-Nodes")
    ControlPointsID = Array{Int32}(ControlPointsID)
    vertices = Array{T,2}(transpose(reshape(ControlPoints,reverse(size(ControlPoints)))))
    verts = GeometryBasics.Point3f[]
    for i in 1:size(vertices,1)
        push!(verts, GeometryBasics.Point3f(vertices[i,:]))
    end
    mesh = GeometryBasics.Mesh(verts,GeometryBasics.faces(mesh))

    # # get nodes and elements IDS from preCICE
    numberOfVertices, dimensions = length(mesh.position), 3

    # integration points
    quadPoints = Array{T,2}(undef, length(mesh), dimensions)
    for i in 1:numberOfVertices
        quadPoints[i,:] .= center(mesh[i])
    end

    # maping from center to nodes
    tmp = mapreduce(((i,F),)->vcat(map(T->(i,T),Base.to_index.(F).data)...),vcat,enumerate(faces(mesh)))
    map_id = map(i->getindex.(tmp,1)[findall(==(i),getindex.(tmp,2))],1:numberOfVertices)

    # forces are located on the nodes this time
    forces = zeros(T, size(ControlPoints))
    σ = Array{T,2}(undef, length(mesh), dimensions)
    srf_id = mapreduce(((i,ids),)->map(T->(i,T),ids),vcat,enumerate(srf_id))
    # return interface
    LumpedInterface(mesh,deepcopy(mesh),srf_id,vertices,ControlPointsID,forces,σ,func,map_id,Int[0.0])
end
function readData!(interface::LumpedInterface)
    # Read control point displacements
    interface.vertices = PreCICE.readData("Solid-Nodes", "Displacements", interface.ControlPointsID, interface.dt[end])
    # update the mesh
    points = Point3f[]
    for (i,pnt) in enumerate(interface.mesh0.position)
        push!(points, Point3f(SA[pnt.data...] .+ interface.vertices[i,:]))
    end
    interface.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(interface.mesh0))
end

function update!(interface::LumpedInterface)
    t = sum(@views(interface.dt[1:end])) # the time
    map_id = interface.map_id # pointer
    for (i,id) in interface.srfID
        interface.σ[id,:] .= dS(@views(interface.mesh[id])).*interface.func(i,t)
    end
    for i in 1:length(interface.mesh.position)
        # 1/3 of the force of each triangle is added to each node
        interface.forces[i,:] .= transpose(sum(interface.σ[map_id[i],:]./3,dims=1))
    end
end

function writeData!(interface::LumpedInterface)
    # write the force at the quad points
    PreCICE.writeData("Solid-Nodes", "Forces", interface.ControlPointsID, interface.forces)
end

import WaterLily
function WaterLily.write!(w,a::LumpedInterface)
    k = w.count[1]
    points = hcat([[p.data...] for p ∈ a.mesh0.position]...)
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, Base.to_index.(face)) for face in faces(a.mesh0)]
    vtk = vtk_grid(w.dir_name*@sprintf("/%s_%06i", w.fname, k), points, cells) 
    for (name,func) in w.output_attrib
        # point/vector data must be oriented in the same way as the mesh
        vtk[name] = ndims(func(a))==1 ? func(a) : permutedims(func(a))
    end
    vtk_save(vtk); w.count[1]=k+1
    w.collection[round(sum(a.dt),digits=4)]=vtk
end
