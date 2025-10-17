using WaterLilyPreCICE,OrdinaryDiffEq
using CSV,DataFrames,Printf

# compute the volume of the mesh, assumes the sphere is centered at 0 and is uniformly deformed
@inline get_volume(mesh) = 4/3*π*(√sum(abs2, mesh.position[1] .- 0))^3

let
    # initialized precice and return mesh, mapping and forces
    interface = LumpedInterface(surface_mesh="../CalculiX/geom.inp")

    # iteration storage
    storage_step = []
    Cₕ = 0.180                # relaxation factor for the pressure
    mmHg2Pa = 133.322387415

    # setup
    PLV₁ = PLV₀ = 1.0         # must be nonzero
    Pact = 0.0
    EDV = 120.0

    # solver parameters
    simple = true             # use a simple fixed-point or a secant method
    ACTUATE = false           # do we provide an external pressure?

    # initialise
    scale_vol = 90
    vol0 = scale_vol*get_volume(interface.mesh)
    # print zero pressure volume
    @printf("Initial volume: %.2f ml\n", vol0)

    # main time loop
    while PreCICE.isCouplingOngoing()

        # the actuation pressure at the aortic pressure must match the EDV of 120ml at sum(time) = t
        ACTUATE && (Pact = sum(interface.Δt)*68.0)

        # read the data from the other participant, set sum(time) = t+Δt (end of this time step)
        readData!(interface)

        # volume of the sphere
        vol = scale_vol*get_volume(mesh)

        # target volume
        VLV_0D = vol0 + sum(interface.Δt)*(EDV-vol0) # target EDV is 120ml

        # fixed-point for the pressure
        PLV₁ = max(PLV₀ + Cₕ*(VLV_0D - vol), 0.001)
        interface.step==1 && (PLV₁ = P₀) # first time step, used EDP to get EDV
        PLV₀ = PLV₁ # for next iteration or next time step

        # we then need to recompute the forces with the correct volume and pressure
        compute_forces!(forces; p=mmHg2Pa*(PLV₁-Pact))

        # write data to the other participant, advance coupling
        writeData!(interface)

        # save the data
        push!(storage_step, [interface.step, interface.iteration, sum(interface.Δt), PLV₁, Pact, vol,
                             VLV_0D, 0., 0., 0., 0.])

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            out_data = reduce(vcat,storage_step')
            CSV.write("sphere_output.csv", DataFrame(out_data,:auto),
                      header=["timestep","iter","time","PLV_3D", "PACT_3D", "VLV_3D",
                              "VLV_0D", "PAO_0D","QAO_0D", "QMV_0D", "PLV_0D"])
        end
    end
    # save the nodal force at the end of the simulation
    # open("../CalculiX/cload.nam", "w") do io
    #     for i in 1:size(forces,1), j in 1:3
    #         println(io, @sprintf("%d, %d, %1.8f", i, j, forces[i,j]))
    #     end
    # end
    PreCICE.finalize()
end