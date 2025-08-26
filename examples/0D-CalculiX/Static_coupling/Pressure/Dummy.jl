using GeometryBasics,WriteVTK,StaticArrays,PreCICE
using LinearAlgebra: cross
using Printf

function load_inp(fname; facetype=GLTriangleFace, pointtype=Point3f)
    #INP file format
    @assert endswith(fname,".inp") "file type not supported"
    fs = open(fname)

    points = pointtype[]
    faces = facetype[]
    node_idx = Int[]
    srf_id = Tuple[]
    cnt = 0

    # read the first 3 lines if there is the "*heading" keyword
    line = readline(fs)
    contains(line,"*heading") && (line = readline(fs))
    BlockType = contains(line,"*NODE") ? Val{:NodeBlock}() : Val{:DataBlock}()

    # read the file
    while !eof(fs)
        line = readline(fs)
        contains(line,"*ELSET, ELSET=") && (cnt+=1)
        BlockType, line = parse_blocktype!(BlockType, fs, line)
        if BlockType == Val{:NodeBlock}()
            push!(node_idx, parse(Int,split(line,",")[1])) # keep track of the node index of the inp file
            push!(points, pointtype(parse.(eltype(pointtype),split(line,",")[2:4])))
        elseif BlockType == Val{:ElementBlock}()
            nodes = parse.(Int,split(line,",")[2:end])
            push!(faces, TriangleFace{Int}(facetype([findfirst(==(node),node_idx) for node in nodes])...)) # parse the face
        elseif BlockType == Val{:ElSetBlock}()
            push!(srf_id, (cnt, parse.(Int,split(line,",")[1])));
        else
            continue
        end
    end
    close(fs) # close file stream
    # case where there is no surface ID and we pass it a single tuple of all the faces
    srf_id = length(srf_id)==0 ? ntuple(i->(1,i),length(faces)) : ntuple(i->(srf_id[i].+(0,1-srf_id[1][2])),length(srf_id))
    return Mesh(points, faces), srf_id
end
function parse_blocktype!(block, io, line)
    contains(line,"*NODE") && return block=Val{:NodeBlock}(),readline(io)
    contains(line,"*ELEMENT") && return block=Val{:ElementBlock}(),readline(io)
    contains(line,"*ELSET, ELSET=") && return block=Val{:ElSetBlock}(),readline(io)
    return block, line
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

let

    # keyword argument might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # load the mesh
    mesh0,srf_id = load_inp("../Solid/geom.inp")

    # create coupling
    PreCICE.createParticipant("LPM", configFileName, 0, 1)

    # initialise PreCICE
    PreCICE.initialize()

    # we need to initialize before we can get the mesh points and coordinates
    (ControlPointsID, ControlPoints) = PreCICE.getMeshVertexIDsAndCoordinates("Solid-Mesh")
    ControlPointsID = Array{Int32}(ControlPointsID)
    vertices = Array{Float64,2}(transpose(reshape(ControlPoints,reverse(size(ControlPoints)))))
    verts = GeometryBasics.Point3f[]
    for i in 1:size(vertices,1)
        push!(verts, GeometryBasics.Point3f(vertices[i,:]))
    end
    mesh = GeometryBasics.Mesh(verts,GeometryBasics.faces(mesh0))

    # mapping from center to nodes, needed for the forces
    forces = zeros(Float64, size(ControlPoints))

    # map triangle to nodes, needed for the forces
    map_id = map(((i,F),)->vcat(Base.to_index.(F).data...),enumerate(faces(mesh)))

    # solver setting
    solver_dt = 1.0
    time = [0.0]

    # initialise storage
    mesh_storage = deepcopy(mesh)
    r0 = √sum(abs2, mesh.position[1] .- 0)
    volume = [4/3*π*r0^3]
    pressure = []
    Q = []

    # iteration storage
    pressure_iter = []
    volume_iter = []
    p₁ = p₀ = 0.
    iteration = 1.0
    Q_target = 0.5
    Cₕ = 0.025 # relaxation factor for the pressure

    # main time loop
    while PreCICE.isCouplingOngoing()

        if PreCICE.requiresWritingCheckpoint()
            # save state at sum(time) = t
            mesh_storage = deepcopy(mesh)
        end

        # set time step
        precice_dt = PreCICE.getMaxTimeStepSize()
        dt = min(precice_dt, solver_dt)
        push!(time, dt) # sum(time) = t+Δt (end of this time step)

        # read data from other participant
        displacements = PreCICE.readData("Solid-Mesh", "Displacements", ControlPointsID, dt)
        # update the mesh
        mesh = update_mesh(mesh0, displacements)
        r = √sum(abs2, mesh.position[1] .- 0)
        println(" New radius at t=$(sum(time)): ", r)

        # linear pressure ramp from 0 to 1
        # p₁ = 0.01*sum(time)

        # volume of sphere
        push!(volume_iter, 4/3*π*r^3)

        # target volume
        target_volume = volume[1] + Q_target*sum(time) # target flow rate is 0.5
        println(" Current volume: ",volume_iter[end]," target volume: ",target_volume)

        if iteration < 3
            p₁ = p₀ = Cₕ*(target_volume - volume_iter[end])
        else
            p₁ = p₀ + Cₕ*(target_volume - volume_iter[end])
            ∂p = pressure_iter[end] - pressure_iter[end-1]
            ∂v = volume_iter[end] - volume_iter[end-1]
            println(" ∂p/∂v: ", round(∂p/∂v, digits=6), "   Cₕ: ", round(Cₕ, digits=6))
            # p₁ = p₀ + (∂p/∂v)*(target_volume - volume[end])
            p₀ = p₁
        end
        push!(pressure_iter, p₁)

        println(" P₁: ", p₁)
        println(" relative volume residual: ", abs(target_volume - volume[end])/target_volume)

        # we then need to recompute the forces with the correct volume and pressure
        forces .= 0 # reset the forces
        for (i,id) in srf_id
            f = dS(@views(mesh[id])) .* p₁ # pressure is 0.01, regardless of the time
            forces[map_id[id],:] .+= transpose(f)./3 # add all the contribution from the faces to the nodes
        end

        # write the force at the nodes
        PreCICE.writeData("Solid-Mesh", "Forces", ControlPointsID, forces)

        # do the coupling
        PreCICE.advance(dt) # advance to t+Δt

        # read checkpoint if required or move on
        if PreCICE.requiresReadingCheckpoint()
            iteration += 1
            # revert to sum(time) = t
            mesh = deepcopy(mesh_storage)
            pop!(time) # remove last time step since we are going back
        else
            # only update the mesh once per time window
            # mesh = update_mesh(mesh0, displacements)
            pressure_iter = [] # reset
            volume_iter = []
            iteration = 1
            push!(volume, 4/3*π*r^3)
            push!(pressure, p₁)
            push!(Q, (volume[end]-volume[end-1])/dt)
        end

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            vol = 4/3*π*(√sum(abs2, mesh.position[1] .- 0))^3
            println(" final relative volume residual: ", abs(target_volume - volume[end])/target_volume)
            println()
        end
    end
    @show volume
    @show pressure
    @show Q
    # save the nodal force at the end of the simulation
    # open("../Solid/cload.nam", "w") do io
    #     for i in 1:size(forces,1), j in 1:3
    #         println(io, @sprintf("%d, %d, %1.8f", i, j, forces[i,j]))
    #     end
    # end
    PreCICE.finalize()
end