using PreCICE
using CirculatorySystemModels
using ModelingToolkit
using OrdinaryDiffEq

mutable struct LumpedInterface{T} <: AbstractInterface
    mesh0           :: GeometryBasics.Mesh # initial geometry, never changed
    mesh            :: GeometryBasics.Mesh
    srf_el          :: NTuple
    deformation     :: AbstractArray # might not be needed
    ControlPointsID :: AbstractArray
    forces          :: AbstractArray
    func            :: Function
    map_id          :: AbstractVector
    dt              :: Vector{T}  # time step vector
    Q               :: Vector{T}  # flow rate vector
    P               :: Vector{T}  # pressure vector
    integrator      :: Union{Nothing,OrdinaryDiffEq.ODEIntegrator}
    mesh_store      :: GeometryBasics.Mesh # storage
end

function LumpedInterface(T=Float64; surface_mesh="../Solid/geom.inp", func=(i,t)->0, prob=nothing, kwargs...)

    # keyword aguments might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # create coupling
    PreCICE.createParticipant("LPM", configFileName, 0, 1)
    
    # load the file
    mesh,srf_id = load_inp(surface_mesh) # can we get rid of this?

    # check if we need to initialize the data
    if PreCICE.requiresInitialData()
        @assert true "this is not true"
    end
        
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

    # mapping from center to nodes, needed for the forces
    forces = zeros(T, size(ControlPoints))
    map_id = map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(mesh)))
    
    # generate lumped model, if the 'prob' is not provided, we return a nothing
    integrator = init(prob, Vern7(), reltol=1e-6, abstol=1e-9)
    u₀ = [0.,1.,2.] #deepcopy(integrator.u0)

    # return interface
    LumpedInterface(mesh,deepcopy(mesh),srf_id,vertices,ControlPointsID,
                    forces,func,map_id,T[0],T[0],T[0],integrator,deepcopy(mesh))
end

# binding 
OrdinaryDiffEq.init(a::Nothing, args...; kwargs...) = nothing # no init if now problem is provided
OrdinaryDiffEq.step!(a::Nothing, args...; kwargs...) = nothing # no step if now integrator is provided

"""
    readData!(::LumpedInterface)

Read the coupling data (displacements) and update the mesh position.
"""
function readData!(interface::LumpedInterface)
    # set time step
    dt_precice = PreCICE.getMaxTimeStepSize()
    #@TODO get max timestep from Lumped model
    push!(interface.dt, min(10, dt_precice)) # min physical time step

    if PreCICE.requiresWritingCheckpoint()
        # save the mesh at this step
        interface.mesh_store = deepcopy(interface.mesh)
        # save initial condition of ODE solver
    end
    # Read control point displacements
    interface.deformation = PreCICE.readData("Solid-Nodes", "Displacements",
                                             interface.ControlPointsID, interface.dt[end])
end

"""
    update!(::LumpedInterface)

Updates the interface conditions (the forces) from the interface function.
"""
function update!(interface::LumpedInterface)
    # store the volume for flow rate computation
    Vᵢ = WaterLilyPreCICE.volume(interface.mesh)[1]
    # update the mesh
    points = Point3f[]
    for (i,pnt) in enumerate(interface.mesh0.position)
        push!(points, Point3f(SA[pnt.data...] .+ interface.deformation[i,:]))
    end
    interface.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(interface.mesh0))
    # update flow rate @TODO, can we do that on the displacements directly?
    push!(interface.Q, (Vᵢ-WaterLilyPreCICE.volume(interface.mesh)[1]) / interface.dt[end])
    # update 0D model
    # interface.integrator... = interface.Q # modify flow rate
    # OrdinaryDiffEq.step!(interface.integrator, interface.dt[end], false)
    # compute forces
    get_forces!(interface)
end

function get_forces!(interface::LumpedInterface, t=sum(@views(interface.dt)); kwargs...)
    # compute nodal forces
    interface.forces .= 0 # reset the forces
    for (i,id) in interface.srf_el
        f = dS(@views(interface.mesh[id])).*interface.func(i,t,interface)
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
    
    # do the coupling
    PreCICE.advance(interface.dt[end])
    
    # read checkpoint if required or move on
    if PreCICE.requiresReadingCheckpoint()
        # revert the mesh
        interface.mesh = deepcopy(interface.mesh_store)
        # pop the flux and pressures
        pop!(interface.dt); pop!(interface.Q); pop!(interface.P)
        # revert state of ODE solver
    end
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
