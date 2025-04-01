using PreCICE
using CirculatorySystemModels
using ModelingToolkit
using OrdinaryDiffEq

mutable struct LumpedInterface{T} <: AbstractInterface
    mesh0           :: GeometryBasics.Mesh # initial geometry, never changed
    mesh            :: GeometryBasics.Mesh
    srf_id          :: NTuple
    deformation     :: AbstractArray{T} # might not be needed
    ControlPointsID :: AbstractArray{Int32}
    forces          :: AbstractArray{T}
    func            :: Function
    map_id          :: AbstractVector
    dt              :: AbstractVector{T}  # time step vector
    V               :: AbstractVector{T}  # flow rate vector
    P               :: AbstractVector{T}  # pressure vector
    integrator      :: Union{Nothing,OrdinaryDiffEq.ODEIntegrator}
    u₀              :: Union{Nothing,AbstractVector{T}}
    sol             :: AbstractVector
    mesh_store      :: GeometryBasics.Mesh # storage
end

function LumpedInterface(T=Float64; surface_mesh="../Solid/geom.inp", func=(i,t)->0, integrator=nothing, kwargs...)

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
    vertices .= 0 # reset since tey become the displacements
    map_id = map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(mesh)))
    
    # time step
    dt_precice = PreCICE.getMaxTimeStepSize()

    # generate lumped model, if the 'integrator' is not provided, we return a nothing
    u₀ = isnothing(integrator) ? nothing : deepcopy(integrator.u)
    V₀ = T[sum(volume(mesh))./3.]

    # return interfaceFQ
    LumpedInterface(mesh,deepcopy(mesh),srf_id,vertices,ControlPointsID,
                    forces,func,map_id,T[0],V₀,T[0],integrator,u₀,[],deepcopy(mesh))
end

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
        !isnothing(interface.integrator) && (interface.u₀ = deepcopy(interface.integrator.u))
    end
    # Read control point displacements
    interface.deformation .= PreCICE.readData("Solid-Nodes", "Displacements",
                                              interface.ControlPointsID, interface.dt[end])
end

"""
    update!(::LumpedInterface)

Updates the interface conditions (the forces) from the interface function.
"""
function update!(interface::LumpedInterface)
    # update the mesh
    points = Point3f[]
    for (i,pnt) in enumerate(interface.mesh0.position)
        push!(points, Point3f(SA[pnt.data...] .+ interface.deformation[i,:]))
    end
    interface.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(interface.mesh0))
    # stores the volume
    push!(interface.V, volume(interface))
    # update 0D model
    if !isnothing(interface.integrator)
        SciMLBase.set_u!(interface.integrator, [get_Q(interface),interface.u₀[2:end]...])
        integrate!(interface)
        push!(interface.P, interface.integrator.u[1])
    end
    # compute forces
    get_forces!(interface)
end

"""
    integrate!(a::LumpedInterface)

Integrate the 0D model, this function will step the ODE solver from integrator.t to 
integrator.t + a.dt[end] and stops exaclty there. It will also store the solution at
each coupling step.
"""
integrate!(a::LumpedInterface) = (OrdinaryDiffEq.step!(a.integrator, a.dt[end], true);
                                  push!(a.sol, [a.integrator.t, a.integrator.u]))

"""
    get_Q(::LumpedInterface)

Get the flow rate for the 0D model, either using a 1st or 2nd order scheme.
"""
get_Q(a::LumpedInterface) = length(a.V)<3 ? sum(@view(a.V[end-1:end]).*SA[-1.,1])/a.dt[end] : sum(@view(a.V[end-2:end]).*SA[1.,-4.,3.])/2a.dt[end]

"""
    get_forces!(::LumpedInterface)

Compute the forces on the interface.
"""
function get_forces!(interface::LumpedInterface, t=sum(@views(interface.dt)); kwargs...)
    # compute nodal forces
    interface.forces .= 0 # reset the forces
    for (i,id) in interface.srf_id
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
        pop!(interface.dt); pop!(interface.V); pop!(interface.P)
        # revert state of ODE solver
        if !isnothing(interface.integrator)
            pop!(interface.sol)
            SciMLBase.set_ut!(interface.integrator, interface.u₀[2], interface.u₀[1])
        end
    end
end

# volume computation for lumped interface, needs to switch signs for correct contributions
function volume(a::LumpedInterface)
    vol = zeros(SVector{3,Float64})
    max_id = getindex(maximum(a.srf_id),1)
    mn,mx = getindex.(extrema(a.srf_id),1)
    Np = (mx-2)/2
    srf = vcat(collect(1:Np+1),collect(2Np+2:3Np+2)) # inner surfaces only
    for (i,T) in enumerate(a.mesh)
        id = getindex(a.srf_id[i],1)
        sgn = ifelse(id ≤ max_id÷2, 1, -1)
        id ∈ srf && (vol += sgn .* WaterLilyPreCICE.center(T) .* WaterLilyPreCICE.dS(T))
    end
    return sum(vol)/3
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
