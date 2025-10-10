# Interface

struct Interface <: AbstractInterface
    ControlPointsID :: AbstractArray
    forces          :: AbstractArray
    deformation     :: AbstractArray
    map_id          :: AbstractVector
    dt              :: Vector{Float64}
    rw_mesh         :: String
    read_data       :: String
    write_data      :: String
end


"""
    Interface(T=Float64;
              surface_mesh="geom.inp",
              center=0,
              scale=1.f0,
              boundary=true,
              thk=0,
              passive_bodies=nothing,
              rw_mesh="Fluid-Mesh",
              read_data="Displacements",
              write_data="Forces",
              kwargs...)

Constructor for a coupling Interface that uses PreCICE for coupling WaterLily:
    - `T`             : Type of the interface, default is Float64.
    - `surface_mesh`  : Path to the surface mesh file, default is "geom.inp".
    - `center`        : Center of the solid, default is 0. Can be used to move the structure in the domain.
    - `scale`         : Scaling factor for the mesh, default is 1.0.
    - `boundary`      : Is the mesh provided the boundary of the solid, default is true.
    - `thk`           : Thickness of the solid, default is 0. If `boundary` is false, this is used to set the thickness.
    - `passive_bodies`: Passive bodies to add to the interface, they are immersed, but no FSI occurs for those.
    - `rw_mesh`       : Name of the mesh in PreCICE, default is "Fluid-Mesh".
    - `read_data`     : Name of the data to read from PreCICE, default is "Displacements".
    - `write_data`    : Name of the data to write to PreCICE, default is "Forces".

"""
function Interface(T=Float64; surface_mesh="geom.inp", center=0, scale=1.f0, boundary=true, thk=0, passive_bodies=nothing,
                    rw_mesh="Fluid-Mesh", read_data="Displacements", write_data="Forces", kwargs...)

    # load the file
    mesh, srf_id = load_inp(surface_mesh) # can we get rid of this?

    # pass nodes and elements IDs to preCICE, here we use the unscaled and un-mapped mesh
    numberOfVertices, dimensions = length(mesh.position), 3
    vertices = Array{Float64,2}(undef, dimensions, numberOfVertices)
    for i in 1:numberOfVertices
        vertices[:,i] .= mesh.position[i].data # we need to have the same scale as in CalculiX
    end
    ControlPointsID = PreCICE.setMeshVertices(rw_mesh, reshape(vertices, (:,3)))

    # construct the actual mesh now that the mapping has been computed
    verts = Point{3,T}[]
    for i in 1:numberOfVertices
        # prepare the mesh, here we move it to the center of the domain
        push!(verts, Point{3,T}(vertices[:,i].*scale .+ center))
    end
    mesh = GeometryBasics.Mesh(verts,GeometryBasics.faces(mesh))
    bbox = Rect(mesh.position)
    bbox = Rect(bbox.origin.-max(4,thk),bbox.widths.+max(8,2thk))
    velocity = GeometryBasics.Mesh(zero(verts),GeometryBasics.faces(mesh))
    body = MeshBody(mesh,deepcopy(mesh),velocity,srf_id,(x,t)->x,bbox,T(scale),T(thk/2),boundary)

    # storage arrays, @TODO these needs to be float64 for precice
    forces = zeros(Float64, numberOfVertices, dimensions)
    deformation = zeros(Float64, numberOfVertices, dimensions)

    # mapping from center to nodes, needed for the forces
    map_id = Base.map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(body.mesh)))

    # initilise PreCICE
    PreCICE.initialize()
    dt = PreCICE.getMaxTimeStepSize()

    # add some passive_bodies if we need
    !isnothing(passive_bodies) && for b in passive_bodies
        body += b # use Setbody
    end

    # return coupling interface
    interface = Interface(ControlPointsID, forces, deformation, map_id, [dt], rw_mesh, read_data, write_data)
    return interface, body
end

function readData!(interface::S) where S<:Interface
    # Read control point displacements
    interface.deformation .= PreCICE.readData(interface.rw_mesh, interface.read_data,
                                              interface.ControlPointsID, interface.dt[end])
end

"""
    update!(::MeshBody{T}, ::Interface; kwargs...)

Updates a `MeshBody` with the current state of the `Interface`. The mesh position is updated based on the deformation
field provided by the `Interface`. The mesh is assumed to be in the reference configuration, and the deformation is applied to it.
"""
function update!(body::MeshBody{T}, interface::S, dt; kwargs...) where {S<:Interface,T}
    # update mesh position, measure is done elsewhere
    points = Point{3,T}[]
    for (i,pnt) in enumerate(body.mesh0.position) # initial mesh is in the ref config.
        push!(points, Point{3,T}(SA[pnt.data...] .+ body.scale.*interface.deformation[i,:]))
    end
    # update the MeshBody with the new points, time here must be converted into fluid time
    update!(body, points, dt)
end
function update!(body::WaterLily.SetBody, interface::S, dt; kwargs...) where S<:Interface
    update!(body.a, interface, dt; kwargs...)
    update!(body.b, interface, dt; kwargs...)
end
function update!(::AutoBody, ::Interface, dt; kwargs...) end # do nothing for other bodies

"""
    get_forces!(interface::S, flow::Flow, body::MeshBody; δ=1.f0, kwargs...)

Computes the forces on the mesh body based on the flow field and the interface. The forces are computed at the integration points
and mapped to the nodes of the mesh body. The forces are stored in the `interface.forces` array, which is reset before computation.
"""
function compute_forces!(interface::S, flow::Flow, body::MeshBody; δ=1, kwargs...) where S<:Interface
    # how many nodes per face
    N = length(interface.map_id[1])
    interface.forces .= 0 # reset the forces
    # compute nodal forces
    for id in 1:length(body.mesh)
        tri = body.mesh[id]
        # map into correct part of the mesh
        f = get_p(tri, flow.p, δ, Val{body.boundary}())
        interface.forces[interface.map_id[id],:] .-= transpose(f)./N # add all the contribution from the faces to the nodes
    end
end
function compute_forces!(interface::S, flow::Flow, body::WaterLily.SetBody; δ=1.f0, kwargs...) where S<:Interface
    compute_forces!(interface, flow, body.a; δ, kwargs...)
    compute_forces!(interface, flow, body.b; δ, kwargs...)
end
function compute_forces!(::Interface, ::Flow, ::AutoBody; kwargs...) end # skip if not a MeshBody

"""
    writeData!(interface::S)

Writes the forces at the integration points to the PreCICE interface. The forces are written to the mesh defined by `interface.rw_mesh`
"""
function writeData!(interface::S) where S<:Interface
    # write the force at the integration points
    PreCICE.writeData(interface.rw_mesh, interface.write_data, interface.ControlPointsID, interface.forces)
end
