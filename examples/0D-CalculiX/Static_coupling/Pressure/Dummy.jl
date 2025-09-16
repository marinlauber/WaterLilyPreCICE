using GeometryBasics,WriteVTK,StaticArrays,PreCICE
using Printf,CSV,DataFrames
using LinearAlgebra: cross

# reading inp files
include("IO.jl")
include("utils.jl")

let

    # keyword argument might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # initialized precice and return mesh, mapping and forces
    mesh0, map_id, srf_id, forces, ControlPointsID = generate(configFileName)

    # solver setting
    solver_dt = 1.0
    time = [0.0]

    # initialise storage
    scale_vol = 80
    vol0 = scale_vol*get_volume(mesh0)

    # iteration storage
    storage_step = []
    PLV₁ = PLV₀ = 1.0 # must be nonzero
    Pact = 0.0
    iteration = step = 1
    EDV = 120.0
    Cₕ = 0.1 # relaxation factor for the pressure
    simple = true # use a simple fixed-point or a secant method
    mmHg2Pa = 133.322387415
    ACTUATE = false # do we provide an external pressure?

    # main time loop
    while PreCICE.isCouplingOngoing()

        if PreCICE.requiresWritingCheckpoint()
            # save state at sum(time) = t
        end

        # set time step
        precice_dt = PreCICE.getMaxTimeStepSize()
        dt = min(precice_dt, solver_dt)
        push!(time, dt) # sum(time) = t+Δt (end of this time step)

        # read data from other participant
        displacements = PreCICE.readData("Solid-Mesh", "Displacements", ControlPointsID, dt)
        # update the mesh
        mesh = update_mesh(mesh0, displacements)
        vol = scale_vol*get_volume(mesh)

        # target volume
        VLV_0D = vol0 + sum(time)*(EDV-vol0) # target EDV is 120ml

        # the actuation pressure at the aortic pressure must match the EDV of 120ml
        # ACTUATE && (Pact = sum(time)*(80.0-8.0097))

        # fixed-point for the pressure
        # (Pact > PLV₀) && (PLV₀ += (Pact-PLV₀)) # ensure PLV₀ > Pact
        PLV₁ = PLV₀ + Cₕ*(VLV_0D - vol)
        PLV₀ = PLV₁

        # we then need to recompute the forces with the correct volume and pressure
        compute_forces!(forces, mmHg2Pa*(PLV₁-Pact), mesh, srf_id, map_id)

        # write the force at the nodes
        PreCICE.writeData("Solid-Mesh", "Forces", ControlPointsID, forces)

        # do the coupling
        PreCICE.advance(dt) # advance to t+Δt

        # save the data
        push!(storage_step, [step, iteration, sum(time), PLV₁, Pact, vol,
                             VLV_0D, 0., 0., 0., 0.])

        # read checkpoint if required or move on
        if PreCICE.requiresReadingCheckpoint()
            # revert to sum(time) = t
            pop!(time) # remove last time step since we are going back
            iteration += 1
        else
            # we can move on, reset counter and increment step
            iteration = 1
            step += 1
        end

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            out_data = reduce(vcat,storage_step')
            CSV.write("sphere_output.csv", DataFrame(out_data,:auto),
                      header=["timestep","iter","time","PLV_3D", "PACT_3D", "VLV_3D",
                              "VLV_0D", "PAO_0D","QAO_0D", "QMV_0D", "PLV_0D"])
        end
    end
    # save the nodal force at the end of the simulation
    # open("../Solid/cload.nam", "w") do io
    #     for i in 1:size(forces,1), j in 1:3
    #         println(io, @sprintf("%d, %d, %1.8f", i, j, forces[i,j]))
    #     end
    # end
    PreCICE.finalize()
end