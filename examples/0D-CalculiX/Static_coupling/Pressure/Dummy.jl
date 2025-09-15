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
    vol0 = 100*get_volume(mesh0)

    # iteration storage
    storage_step,pressure_iter,volume_iter = [],[],[]
    PLV₁ = PLV₀ = 0.
    iteration = step = 1
    Q_target = 60.0
    Cₕ = 0.64 # relaxation factor for the pressure
    simple = true # use a simple fixed-point or a secant method

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
        vol = 100*get_volume(mesh)

        # volume of sphere
        push!(volume_iter, vol)

        # target volume
        target_volume = vol0 + Q_target*sum(time) # target flow rate is 0.5

        # fixed-point for the pressure
        if simple==true
            PLV₁ = PLV₀ + Cₕ*(target_volume - volume_iter[end])
        else
            ∂p = pressure_iter[end] - pressure_iter[end-1]
            ∂v = volume_iter[end] - volume_iter[end-1]
            PLV₁ = PLV₀ + (∂p/∂v)*(target_volume - volume[end])
        end
        PLV₀ = PLV₁
        push!(pressure_iter, PLV₁)

        # we then need to recompute the forces with the correct volume and pressure
        compute_forces!(forces, PLV₁, mesh, srf_id, map_id)

        # write the force at the nodes
        PreCICE.writeData("Solid-Mesh", "Forces", ControlPointsID, forces)

        # do the coupling
        PreCICE.advance(dt) # advance to t+Δt

        # save the data
        push!(storage_step, [step, iteration, sum(time), PLV₁, 0.0, vol,
                             0., 0., 0., 0., 0.])

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