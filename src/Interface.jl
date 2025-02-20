using PreCICE
using StaticArrays

# can be cleaner
include("MeshBodies.jl")

struct Interface <: AbstractInterface
    ControlPointsID :: AbstractArray
    forces          :: AbstractArray
    deformation     :: AbstractArray
    map_id          :: AbstractVector
    dt              :: Vector{Float64}
end

function Interface(T=Float64; surface_mesh="geom.inp", center=0, scale=1.f0, boundary=true, thk=0, kwargs...)  
    
    # load the file
    mesh,srf_id = load_inp(surface_mesh) # can we get rid of this?
    
    # pass nodes and elements IDS to preCICE, here we use the unscaled and un-maped mesh
    numberOfVertices, dimensions = length(mesh.position), 3
    vertices = Array{Float64,2}(undef, numberOfVertices, dimensions)
    for i in 1:numberOfVertices
        vertices[i,:] .= mesh.position[i].data # we need to have the same scale as in CalculiX
    end
    ControlPointsID = PreCICE.setMeshVertices("Fluid-Mesh", vertices)

    # construct the actual mesh now that the mapping has been compuated
    verts = GeometryBasics.Point3f[]
    for i in 1:numberOfVertices
        # prepare the mesh, here we move it to the center of the domain
        push!(verts, GeometryBasics.Point3f(vertices[i,:].*scale .+ center))
    end
    mesh = GeometryBasics.Mesh(verts,GeometryBasics.faces(mesh))
    bbox = Rect(mesh.position)
    bbox = Rect(bbox.origin.-max(4,thk),bbox.widths.+max(8,2thk))
    body = MeshBody(mesh,srf_id,(x,t)->x,bbox,T(scale),T(thk/2),boundary)

    # storage arrays
    forces = zeros(Float64, size(vertices))
    deformation = zeros(Float64, size(vertices))

    # mapping from center to nodes, needed for the forces
    map_id = Base.map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(body.mesh)))
    
    # initilise PreCICE
    PreCICE.initialize()
    dt = PreCICE.getMaxTimeStepSize()
    
    # return coupling interface
    interface = Interface(ControlPointsID, forces, deformation, map_id, [dt])
    return interface, body
end

function readData!(interface::Interface)
    # Read control point displacements
    interface.deformation .= PreCICE.readData("Fluid-Mesh", "Displacements", 
                                              interface.ControlPointsID, interface.dt[end])
end

function update!(interface::Interface, sim::CoupledSimulation; kwargs...)
    # update mesh position, measure is done elsewhere
    points = Point3f[]
    for (i,pnt) in enumerate(sim.store.b.mesh.position)
        push!(points, Point3f(SA[pnt.data...] .+ sim.body.scale.*interface.deformation[i,:]))
    end
    # update
    sim.body.mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(sim.body.mesh))
    bbox = Rect(points)
    sim.body.bbox = Rect(bbox.origin.-max(4,2sim.body.half_thk),bbox.widths.+max(8,4sim.body.half_thk))
end

import WaterLily: interp
function get_forces!(interface::Interface, flow::Flow, body::MeshBody; δ=1.f0, kwargs...)
    t = sum(@views(interface.dt[1:end])) # the time
    interface.forces .= 0 # reset the forces
    # compute nodal forces
    for (i,id) in body.srfID
        tri = body.mesh[id]
        # map into correct part of the mesh, time does nothing
        f = get_p(tri,flow.p,δ)
        interface.forces[interface.map_id[id],:] .+= transpose(f)./3 # add all the contribution from the faces to the nodes
    end
end

function writeData!(interface::Interface)
    # write the force at the integration points
    PreCICE.writeData("Fluid-Mesh", "Forces", interface.ControlPointsID, interface.forces)
end
