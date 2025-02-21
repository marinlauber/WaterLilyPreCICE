# CalculiXInterface


# struct CalculiXInterface <: AbstractInterface
#     ControlPointsID::AbstractArray
#     forces::AbstractArray
#     deformation::AbstractArray
#     map_id :: AbstractVector
#     dt::Vector{Float64}
# end
# this is only slighlty different than the classical interface
# struct CalculiXInterface <: AbstractInterface
#    a :: AbstractInterface
# end
# # overlead properties
# Base.getproperty(C::CalculiXInterface, s::Symbol) = s in propertynames(C) ? getfield(C, s) : getfield(C.a, s)

function CalculiXInterface(T=Float64; surface_mesh="geom.inp", center=0, scale=1.f0, boundary=true, thk=0, 
                           rw_mesh="Solid-Mesh", read_data="Displacements", write_data="Forces", kwargs...)  

    # load the file
    mesh,srf_id = load_inp(surface_mesh) # can we get rid of this?
        
    # initialise PreCICE
    PreCICE.initialize()
    
    # we need to initialize before we can get the mesh points and coordinates
    (ControlPointsID, ControlPoints) = PreCICE.getMeshVertexIDsAndCoordinates(rw_mesh)
    ControlPointsID = Array{Int32}(ControlPointsID)
    vertices = Array{T,2}(transpose(reshape(ControlPoints,reverse(size(ControlPoints)))))
    verts = GeometryBasics.Point3f[]
    for i in 1:size(vertices,1)
        # prepare the mesh, here we move it to the center of the domain
        push!(verts, GeometryBasics.Point3f(vertices[i,:].*scale .+ center))
    end
    mesh = GeometryBasics.Mesh(verts,GeometryBasics.faces(mesh))
    bbox = Rect(mesh.position)
    bbox = Rect(bbox.origin.-max(4,thk),bbox.widths.+max(8,2thk))
    body = MeshBody(mesh,deepcopy(mesh),srf_id,(x,t)->x,bbox,T(scale),T(thk/2),boundary)

    # storage arrays
    forces = zeros(Float64, numberOfVertices, dimensions)
    deformation = zeros(Float64, numberOfVertices, dimensions)

    # mapping from center to nodes, needed for the forces
    map_id = Base.map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(body.mesh)))
    
    # initilise PreCICE
    dt = PreCICE.getMaxTimeStepSize()
    
    # return coupling interface
    interface = Interface(ControlPointsID, forces, deformation, map_id, [dt], rw_mesh, read_data, write_data)
    return interface, body
end

# function readData!(interface::CalculiXInterface)
#     # Read control point displacements
#     interface.deformation .= PreCICE.readData("Solid-Mesh", "Displacements", 
#                                               interface.ControlPointsID, interface.dt[end])
# end

# function update!(interface::S, sim::CoupledSimulation; kwargs...) where S<:Union{Interface,CalculiXInterface}
#     # update mesh position, measure is done elsewhere
#     points = Point3f[]
#     for (i,pnt) in enumerate(sim.store.b.mesh.position)
#         push!(points, Point3f(SA[pnt.data...] .+ sim.body.scale.*interface.deformation[i,:]))
#     end
#     # update
#     sim.body.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(sim.body.mesh))
#     bbox = Rect(points)
#     sim.body.bbox = Rect(bbox.origin.-max(4,2sim.body.half_thk),bbox.widths.+max(8,4sim.body.half_thk))
# end

# import WaterLily: interp
# function get_forces!(interface::S, flow::Flow, body::MeshBody; δ=1.f0, kwargs...) where S<:Union{Interface,CalculiXInterface}
#     t = sum(@views(interface.dt[1:end])) # the time
#     interface.forces .= 0 # reset the forces
#     # compute nodal forces
#     for id in 1:length(body.mesh)
#         tri = body.mesh[id]
#         # map into correct part of the mesh, time does nothing
#         f = get_p(tri, flow.p, δ)
#         interface.forces[interface.map_id[id],:] .+= transpose(f)./3 # add all the contribution from the faces to the nodes
#     end
# end

# function writeData!(interface::CalculiXInterface)
#     # write the force at the integration points
#     PreCICE.writeData("Solid-Mesh", "Forces", interface.ControlPointsID, interface.forces)
# end
