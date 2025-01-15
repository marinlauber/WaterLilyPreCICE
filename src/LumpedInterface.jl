using PreCICE
using CirculatorySystemModels
using ModelingToolkit
using OrdinaryDiffEq

mutable struct LumpedInterface <: AbstractInterface
    mesh0:: GeometryBasics.Mesh # initial geometry, never changed
    mesh :: GeometryBasics.Mesh
    srfID :: AbstractVector
    vertices :: AbstractArray # might not be needed
    ControlPointsID::AbstractArray
    forces :: AbstractArray
    func :: Function
    map_id :: AbstractVector
    dt :: Vector{Float64}
    Q :: Float64
    P :: Float64
    integrator :: Union{Nothing,OrdinaryDiffEq.ODEIntegrator}
end
export LumpedInterface

function LumpedInterface(;fname="../Solid/geom.inp",T=Float64,func=(i,t)->0,prob=nothing)

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

    # mapping from center to nodes, needed for the forces
    forces = zeros(T, size(ControlPoints))
    map_id = map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(mesh)))

    # link elements to surface IDs
    srf_id = mapreduce(((i,ids),)->map(T->(i,T),ids),vcat,enumerate(srf_id))

    # generate lumped model, if the 'prob' is not provided, we return a nothing
    integrator = init(prob, Vern7(), reltol=1e-6, abstol=1e-9)

    # return interface
    LumpedInterface(mesh,deepcopy(mesh),srf_id,vertices,ControlPointsID,
                    forces,func,map_id,Int[0.0],0.,0.,integrator)
end

# binding 
OrdinaryDiffEq.init(a::Nothing, args...; kwargs...) = nothing # no init if now problem is provided
OrdinaryDiffEq.step!(a::Nothing, args...; kwargs...) = nothing # no step if now integrator is provided

"""
    readData!(::LumpedInterface)

Read the coupling data (displacements) and update the mesh position.
"""
function readData!(interface::LumpedInterface)
    # Read control point displacements
    interface.Q = WaterLilyPreCICE.volume(interface.mesh)[1] # store the volume for flow rate compuation
    interface.vertices = PreCICE.readData("Solid-Nodes", "Displacements", interface.ControlPointsID, interface.dt[end])
    # update the mesh
    points = Point3f[]
    for (i,pnt) in enumerate(interface.mesh0.position)
        push!(points, Point3f(SA[pnt.data...] .+ interface.vertices[i,:]))
    end
    interface.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(interface.mesh0))
    # update flow rate
    interface.Q = (WaterLilyPreCICE.volume(interface.mesh)[1] .- interface.Q) / interface.dt[end]
end

"""
    update!(::LumpedInterface)

Updates the interface conditions (the forces) from the interface function.
"""
function update!(interface::LumpedInterface)
    t = sum(@views(interface.dt[1:end])) # the time
    interface.forces .= 0 # reset the forces
    # update 0D model
    OrdinaryDiffEq.step!(interface.integrator, interface.dt[end], false)
    # compute nodal forces
    for (i,id) in interface.srfID
        f = dS(@views(interface.mesh[id])).*interface.func(i,t)
        interface.forces[interface.map_id[id],:] .+= transpose(f)./3 # add all the contribution from the faces to the nodes
    end
end

"""
    writeData!(::LumpedInterface)

Write the coupling data.
"""
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
