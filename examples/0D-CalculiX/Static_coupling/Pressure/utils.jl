# make the mesh and helper files
function generate(configFileName)
    # load the mesh
    mesh,srf_id = load_inp("../Solid/geom.inp")

    # create coupling
    PreCICE.createParticipant("LPM", configFileName, 0, 1)

    # initialise PreCICE
    PreCICE.initialize()

    # we need to initialize before we can get the mesh points and coordinates
    (ControlPointsID, ControlPoints) = PreCICE.getMeshVertexIDsAndCoordinates("Solid-Mesh")
    ControlPointsID = Array{Int32}(ControlPointsID)

    # mapping from center to nodes, needed for the forces
    forces = zeros(Float64, size(ControlPoints))

    # map triangle to nodes, needed for the forces
    map_id = map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(mesh)))

    # return the good stuff
    return mesh, map_id, srf_id, forces, ControlPointsID
end

# compute forces
function compute_forces!(forces, pressure, mesh, srf_id, map_id)
    forces .= 0 # reset the forces
    for (i,id) in srf_id
        f = dS(@views(mesh[id])) .* pressure # pressure jump
        forces[map_id[id],:] .+= transpose(f)./3 # add all the contribution from the faces to the nodes
    end
end

# oriented area of a triangle
@inbounds @inline dS(tri::GeometryBasics.Ngon{3}) = 0.5f0SVector(cross(tri.points[2]-tri.points[1],tri.points[3]-tri.points[1]))

# update the mesh with the new displacements
function update_mesh(mesh0, deformation)
    # update the mesh so that any measure on it is correct
    points = Point3f[]
    for (i,pnt) in enumerate(mesh0.position)
        push!(points, Point3f(SA[pnt.data...] .+ deformation[i,:]))
    end
    return GeometryBasics.Mesh(points,GeometryBasics.faces(mesh0))
end

# compute the volume of the mesh, assumes the sphere is centered at 0 and is uniformly deformed
@inline get_volume(mesh) = 4/3*π*(√sum(abs2, mesh.position[1] .- 0))^3