using WaterLilyPreCICE,StaticArrays,WriteVTK

# function make_sphere(;L=32,Re=250,U=1)
#     # move the geometry to the center of the domain
#     map(x,t) = x .- SA[L,L,L/2]
#     body = MeshBody(joinpath(@__DIR__,"sphere.inp");map,scale=L/2)
#     # generate sim
#     Simulation((4L,2L,L), (U,0,0), L; ν=U*L/Re, body)
# end

# # make a writer with some attributes to output to the file
# velocity(a::Simulation) = a.flow.u |> Array;
# pressure(a::Simulation) = a.flow.p |> Array;
# _body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array;)
# _vbody(a::Simulation) = a.flow.V |> Array;
# mu0(a::Simulation) = a.flow.μ₀ |> Array;
# custom_attrib = Dict("u" => velocity, "p" => pressure, "d" => _body, "v" => _vbody, "μ₀" => mu0,)

# # make the sim
# sim = make_sphere(L=32)
# # make the paraview writer
# wr = vtkWriter("WaterLily-Dummy";attrib=custom_attrib)

# just a mesh
L,Re,U = 2^5,100,1
body = MeshBody("sphere.inp")

flow = Flow((4L,2L,L), (U,0,0); f=Array, Δt=0.25, ν=0.35)
pois = MultiLevelPoisson(flow.p,flow.μ₀,flow.σ)

let
    # keyword aguments might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end
    # create coupling
    PreCICE.createParticipant("WaterLily", configFileName, 0, 1)

    # # get nodes and elements IDS from precice
    # body = sim.body;
    numberOfVertices, dimensions = length(body.mesh.position), 3
    vertices = Array{Float64,2}(undef, numberOfVertices, dimensions)
    for i in 1:numberOfVertices
        vertices[i,:] .= body.mesh.position[i].data
    end
    ControlPointsID = PreCICE.setMeshVertices("Fluid-Mesh", vertices)

    # # storage arrays
    forces = zero(vertices)

    # initilise PreCICE
    PreCICE.initialize()
    dt = PreCICE.getMaxTimeStepSize()

    # run precice
    iter,every = 0,1
    while PreCICE.isCouplingOngoing()

        # set time step
        @show "getMaxTimeStepSize!"
        dt_precice = PreCICE.getMaxTimeStepSize()
        dt = dt_precice

        # read the data from the other participant
        @show "readData!"
        if PreCICE.requiresWritingCheckpoint()
            @show "store!(store,sim)"
        end
        readData = PreCICE.readData("Fluid-Mesh", "Displacements", ControlPointsID, dt)
        @show readData[1:4,:]

        # update the this participant
        @show "sim_step!"
        @inside flow.p[I] = 0.0
        mom_step!(flow,pois)
        # WaterLily.sim_step!(sim;remeasure=false)
        forces .-= 1.0

        # write data to the other participant
        @show "writeData!" 
        PreCICE.writeData("Fluid-Mesh", "Forces", ControlPointsID, forces)
        
        # advance coupling
        @show "advance!"
        dt = PreCICE.advance(dt)
        # read checkpoint if required or move on
        if PreCICE.requiresReadingCheckpoint()
            @show "revert!(store,sim)"
        end

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            # save the data
            # mod(iter,every)==0 && write!(wr, sim)
            iter += 1
        end
    end
    # end
    # close(wr)
    PreCICE.finalize()
end
println("WaterLily: Closing Julia solver...")