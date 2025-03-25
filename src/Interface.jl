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

function Interface(T=Float64; surface_mesh="geom.inp", center=0, scale=1.f0, boundary=true, thk=0, 
                    rw_mesh="Fluid-Mesh", read_data="Displacements", write_data="Forces", kwargs...)  
    
    # load the file
    mesh, srf_id = load_inp(surface_mesh) # can we get rid of this?
    
    # pass nodes and elements IDS to preCICE, here we use the unscaled and un-maped mesh
    numberOfVertices, dimensions = length(mesh.position), 3
    vertices = Array{Float64,2}(undef, dimensions, numberOfVertices)
    for i in 1:numberOfVertices
        vertices[:,i] .= mesh.position[i].data # we need to have the same scale as in CalculiX
    end
    ControlPointsID = PreCICE.setMeshVertices(rw_mesh, reshape(vertices, (:,3)))

    # construct the actual mesh now that the mapping has been compuated
    verts = GeometryBasics.Point3f[]
    for i in 1:numberOfVertices
        # prepare the mesh, here we move it to the center of the domain
        push!(verts, GeometryBasics.Point3f(vertices[:,i].*scale .+ center))
    end
    mesh = GeometryBasics.Mesh(verts,GeometryBasics.faces(mesh))
    bbox = Rect(mesh.position)
    bbox = Rect(bbox.origin.-max(4,thk),bbox.widths.+max(8,2thk))
    velocity = GeometryBasics.Mesh(zero(verts),GeometryBasics.faces(mesh))
    body = MeshBody(mesh,deepcopy(mesh),velocity,srf_id,(x,t)->x,bbox,T(scale),T(thk/2),boundary)

    # storage arrays
    forces = zeros(Float64, numberOfVertices, dimensions)
    deformation = zeros(Float64, numberOfVertices, dimensions)

    # mapping from center to nodes, needed for the forces
    map_id = Base.map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(body.mesh)))
    
    # initilise PreCICE
    PreCICE.initialize()
    dt = PreCICE.getMaxTimeStepSize()
    
    # return coupling interface
    interface = Interface(ControlPointsID, forces, deformation, map_id, [dt], rw_mesh, read_data, write_data)
    return interface, body
end

function readData!(interface::S) where S<:Interface
    # Read control point displacements
    interface.deformation .= PreCICE.readData(interface.rw_mesh, interface.read_data, 
                                              interface.ControlPointsID, interface.dt[end])
end

function update!(interface::S, sim::CoupledSimulation; kwargs...) where S<:Interface
    # update mesh position, measure is done elsewhere
    points = Point3f[]
    for (i,pnt) in enumerate(sim.body.mesh0.position) # initial mesh is in the ref config.
        push!(points, Point3f(SA[pnt.data...] .+ sim.body.scale.*interface.deformation[i,:]))
    end
    # update
    sim.body.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(sim.body.mesh))
    bbox = Rect(points)
    sim.body.bbox = Rect(bbox.origin.-max(4,2sim.body.half_thk),bbox.widths.+max(8,4sim.body.half_thk))
end

import WaterLily: interp
function get_forces!(interface::S, flow::Flow, body::MeshBody; δ=1.f0, kwargs...) where S<:Interface
    t = sum(@views(interface.dt[1:end])) # the time
    interface.forces .= 0 # reset the forces
    # compute nodal forces
    for id in 1:length(body.mesh)
        tri = body.mesh[id]
        # map into correct part of the mesh, time does nothing
        f = get_p(tri, flow.p, δ, Val{body.boundary}())
        interface.forces[interface.map_id[id],:] .-= transpose(f)./3 # add all the contribution from the faces to the nodes
    end
end

function writeData!(interface::S) where S<:Interface
    # write the force at the integration points
    PreCICE.writeData(interface.rw_mesh, interface.write_data, interface.ControlPointsID, interface.forces)
end
