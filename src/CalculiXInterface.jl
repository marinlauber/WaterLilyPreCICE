# CalculiXInterface

# this is a constructor for a the calculix interface, which is the same as the classical interface
# except that the mesh is read from a CalculiX input file directly, and not from the surface file
function CalculiXInterface(T=Float64; surface_mesh="geom.inp", center=0, scale=1.f0, boundary=true, thk=0, 
                           rw_mesh="Solid-Mesh", read_data="Displacements", write_data="Forces", kwargs...)  

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
        # prepare the mesh, here we move it to the center of the domain
        push!(verts, GeometryBasics.Point3f(vertices[i,:].*scale .+ center))
    end
    mesh = GeometryBasics.Mesh(verts,GeometryBasics.faces(mesh0))
    bbox = Rect(mesh.position)
    bbox = Rect(bbox.origin.-max(4,thk),bbox.widths.+max(8,2thk))
    velocity = GeometryBasics.Mesh(zero(verts),GeometryBasics.faces(mesh))
    body = MeshBody(mesh,mesh0,velocity,srf_id,(x,t)->x,bbox,T(scale),T(thk/2),boundary)

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
