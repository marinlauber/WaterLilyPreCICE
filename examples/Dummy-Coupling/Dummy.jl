using WaterLilyPreCICE,StaticArrays

# just a mesh
L,Re,U = 2^5,100,1
map(x,t) = x .- SA[L,L,L/2]
body = MeshBody(joinpath(@__DIR__,"../../meshes/sphere.inp");map,scale=L/2)

let
    # keyword aguments might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # create coupling
    PreCICE.createParticipant("Dummy", configFileName, 0, 1)

    # get nodes and elements IDS from precice
    numberOfVertices, dimensions = length(body.mesh.position), 3
    vertices = Array{Float64,2}(undef, numberOfVertices, dimensions)
    for i in 1:numberOfVertices
        vertices[i,:] .= body.mesh.position[i].data
    end
    ControlPointsID = PreCICE.setMeshVertices("Solid-Mesh", vertices)

    # storage arrays
    deformation = zero(vertices)

    # initilise PreCICE
    PreCICE.initialize()
    dt = PreCICE.getMaxTimeStepSize()

    # run precice
    while PreCICE.isCouplingOngoing()

        # set time step
        dt_precice = PreCICE.getMaxTimeStepSize()
        dt = dt_precice

        # read the data from the other participant
        if PreCICE.requiresWritingCheckpoint()
            @show "store!(store,sim)"
        end
        readData = PreCICE.readData("Solid-Mesh", "Forces", ControlPointsID, dt)

        # update the this participant
        deformation .+= 1.0

        # write data to the other participant
        PreCICE.writeData("Solid-Mesh", "Displacements", ControlPointsID, deformation)
        
        # advance coupling
        dt = PreCICE.advance(dt)
        # read checkpoint if required or move on
        if PreCICE.requiresReadingCheckpoint()
            @show "revert!(store,sim)"
        end
    end
    PreCICE.finalize()
end
println("Dummy: Closing Julia solver...")
