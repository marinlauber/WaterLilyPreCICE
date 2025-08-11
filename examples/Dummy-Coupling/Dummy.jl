using WaterLilyPreCICE,StaticArrays

# just a mesh
L,Re,U = 2^5,100,1
body = MeshBody(joinpath(@__DIR__,"../../meshes/sphere.inp");scale=L/2)

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
    vertices = Array{Float64,2}(undef, dimensions, numberOfVertices)
    # @TODO use mesh0 since it is unscaled and unmapped
    for i in 1:numberOfVertices
        vertices[:,i] .= body.mesh0.position[i].data
    end
    ControlPointsID = PreCICE.setMeshVertices("Solid-Mesh", reshape(vertices, (:,3)))

    # storage arrays
    deformation = zeros(Float64, numberOfVertices, dimensions)

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
        deformation[:,1] .+= dt/L # small increment

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
