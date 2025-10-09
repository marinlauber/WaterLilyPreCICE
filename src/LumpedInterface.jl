using OrdinaryDiffEq
import OrdinaryDiffEqCore: ODEIntegrator

mutable struct LumpedInterface{T} <: AbstractInterface
    mesh0           :: GeometryBasics.Mesh # initial geometry, never changed
    mesh            :: GeometryBasics.Mesh
    srf_id          :: NTuple
    deformation     :: AbstractArray{T} # might not be needed
    ControlPointsID :: AbstractArray{Int32}
    forces          :: AbstractArray{T}
    func            :: Function
    map_id          :: AbstractVector
    Δt              :: AbstractVector{T}  # time step vector
    V               :: AbstractVector{T}  # flow rate vector
    P               :: AbstractVector{T}  # pressure vector
    integrator      :: Union{Nothing,ODEIntegrator}
    u₀              :: Union{Nothing,AbstractVector{T}}
    sol             :: AbstractVector
    rw_mesh         :: String
    read_data       :: String
    write_data      :: String
    iteration       :: Int
    step            :: Int
end
"""
    LumpedInterface()

Arguments:
- `T`: data type, default is `Float64`
- `surface_mesh`: path to the surface mesh in CalculiX format, default is `"../Solid/geom.inp"`
- `func`: function to compute the pressure force on the surface with id="i", the default is just to apply the pressure
   at the point `p` at time `t`: `(i,t,p)->p` this can be used to switch off and add pressures on different surfaces.
- `integrator`: an `OrdinaryDiffEq.ODEIntegrator` object to solve the 0D model, if `nothing` is provided,
   the interface will not solve any 0D model and just apply the pressure at the interface.
- `participant`: name of the participant in the PreCICE configuration file, default is `"LPM"`
- `rw_mesh`: name of the mesh to read and write data from/to, default is `"Solid-Mesh"`
- `read_data`: name of the data to read from the other participant, default is `"Displacements"`
- `write_data`: name of the data to write to the other participant, default is `"Forces"`
- `kwargs`: additional keyword arguments to pass to PreCICE
"""
function LumpedInterface(T=Float64; surface_mesh="../Solid/geom.inp", func=(i,t,p)->p, integrator=nothing,
                         participant="LPM", rw_mesh="Solid-Mesh", read_data="Displacements", write_data="Forces",
                         kwargs...)

    # keyword aguments might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # create coupling
    PreCICE.createParticipant(participant, configFileName, 0, 1)

    # load the file
    mesh0,srf_id = load_inp(surface_mesh) # can we get rid of this?

    # initialise PreCICE
    PreCICE.initialize()

    # we need to initialize before we can get the mesh points and coordinates
    (ControlPointsID, ControlPoints) = PreCICE.getMeshVertexIDsAndCoordinates(rw_mesh)
    ControlPointsID = Array{Int32}(ControlPointsID)
    vertices = Array{T,2}(transpose(reshape(ControlPoints,reverse(size(ControlPoints)))))
    verts = GeometryBasics.Point3f[]
    for i in 1:size(vertices,1)
        push!(verts, GeometryBasics.Point3f(vertices[i,:]))
    end
    mesh = GeometryBasics.Mesh(verts,GeometryBasics.faces(mesh0))

    # mapping from center to nodes, needed for the forces
    forces = zeros(T, size(ControlPoints))
    vertices .= 0 # reset since they become the displacements
    map_id = map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(mesh)))

    # time step
    dt_precice = PreCICE.getMaxTimeStepSize()

    # generate lumped model, if the 'integrator' is not provided, we return a nothing
    u₀ = isnothing(integrator) ? nothing : deepcopy(integrator.u)
    V₀ = T[sum(volume(mesh))./3.]
    sol = isnothing(integrator) ? [[]] : [[integrator.t, u₀...]]
    Δt = T[0] # current time is t=sum(Δt[1:end-1]), t+Δt = sum(Δt)
    P = T[0]

    # return interfaceFQ
    LumpedInterface(mesh0,deepcopy(mesh),srf_id,vertices,ControlPointsID,
                    forces,func,map_id,Δt,V₀,P,integrator,u₀,sol,rw_mesh,read_data,write_data,1,1)
end

"""
    readData!(::LumpedInterface)

Read the coupling data (displacements) and update the mesh position.
"""
function readData!(interface::LumpedInterface, dt_solver=1)
    # write checkpoint
    if PreCICE.requiresWritingCheckpoint()
        # save initial condition of ODE solver
        !isnothing(interface.integrator) && (interface.u₀ = [interface.integrator.t, interface.integrator.u...])
    end

    # set time step
    dt_precice = PreCICE.getMaxTimeStepSize()
    #@TODO get max timestep from Lumped model
    push!(interface.Δt, min(dt_solver, dt_precice)) # min physical time step

    # Read control point displacements
    interface.deformation .= PreCICE.readData(interface.rw_mesh, interface.read_data,
                                              interface.ControlPointsID, interface.Δt[end])
    # update the mesh so that any measure on it is correct
    points = Point3f[]
    for (i,pnt) in enumerate(interface.mesh0.position)
        push!(points, Point3f(SA[pnt.data...] .+ interface.deformation[i,:]))
    end
    interface.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(interface.mesh0))
end

"""
    update!(::LumpedInterface)

Updates the interface conditions (the forces) from the interface function.
"""
function update!(interface::LumpedInterface, pressure; integrate=true)
    # stores the volume
    push!(interface.V, volume(interface))
    # update 0D model
    (!isnothing(interface.integrator) && integrate) && integrate!(interface)
    # compute forces
    compute_forces!(interface, pressure)
end

"""
    integrate!(a::LumpedInterface)

Integrate the 0D model, this function will step the ODE solver from integrator.t to
integrator.t + a.Δt[end] and stops exaclty there. It will also store the solution at
each coupling step.
"""
function integrate!(a::LumpedInterface, u_new)
    # set initial conditions
    SciMLBase.set_ut!(a.integrator, u_new...)
    # solve that step up to t+Δt
    OrdinaryDiffEq.step!(a.integrator, a.Δt[end], true)
    # save the results
    push!(a.sol, [a.integrator.t, a.integrator.u...])
    push!(a.P, a.integrator.u[1])
end

"""
    get_Q(::LumpedInterface)

Get the flow rate for the 0D model, either using a 1st or 2nd order scheme.
"""
function get_Q(a::LumpedInterface)
    N = length(a.V)
    # dVdt|₀ = 0.0
    N==1 && return 0.0
    # return sum(a.V[end-1:end].*SA[-1.,1])/a.Δt[end]
    N<4 && return sum(@view(a.V[end-1:end]).*SA[-1.,1])/a.Δt[end]
    # assumes that dV/dt|₀ = 0 and used 2nd order scheme
    return sum(@view(a.V[end-2:end]).*SA[1.,-4.,3.])/2a.Δt[end]
end

"""
    compute_forces!(::LumpedInterface)

Compute the forces on the interface at t+Δt.
"""
function compute_forces!(interface::LumpedInterface, pressure; kwargs...)
    # how many nodes per face
    N = length(interface.map_id[1])
    # compute nodal forces
    interface.forces .= 0 # reset the forces
    t = sum(interface.Δt)
    for (i,id) in interface.srf_id
        f = dS(@views(interface.mesh[id])) .* interface.func(i,t,pressure)
        interface.forces[interface.map_id[id],:] .+= transpose(f)./N # add all the contribution from the faces to the nodes
    end
end

"""
    writeData!(::LumpedInterface)

Write the coupling data.
"""
function writeData!(interface::LumpedInterface)
    # write the force at the quad points
    PreCICE.writeData(interface.rw_mesh, interface.write_data,
                      interface.ControlPointsID, interface.forces)

    # do the coupling
    PreCICE.advance(interface.Δt[end]) # advance to t+Δt, check convergence and accelerate data

    #  revert to sum(Δt) = t
    if PreCICE.requiresReadingCheckpoint()
        # remove last time step since we are going back
        pop!(interface.Δt)
        interface.iteration += 1
    else
        # we can move on, reset counter and increment step
        interface.iteration = 1
        interface.step += 1
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
# # https://examples.vtk.org/site/Cxx/GeometricObjects/IsoparametricCellsDemo/
# import WaterLily
# function WaterLily.save!(w, a::LumpedInterface; vtkcell_type=celltype(a.mesh))
#     k = w.count[1]
#     points = hcat([[p.data...] for p ∈ a.mesh0.position]...)
#     cells = [MeshCell(vtkcell_type, Base.to_index.(face)) for face in faces(a.mesh0)]
#     vtk = vtk_grid(w.dir_name*@sprintf("/%s_%06i", w.fname, k), points, cells)
#     for (name,func) in w.output_attrib
#         # point/vector data must be oriented in the same way as the mesh
#         vtk[name] = ndims(func(a))==1 ? func(a) : permutedims(func(a))
#     end
#     vtk_save(vtk); w.count[1]=k+1
#     w.collection[round(sum(a.Δt),digits=4)]=vtk
# end
