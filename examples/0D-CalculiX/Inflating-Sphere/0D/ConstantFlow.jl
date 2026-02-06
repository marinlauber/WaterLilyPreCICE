using WaterLilyPreCICE
using CSV,DataFrames

# compute the volume of the mesh, assumes the sphere is centered at 0 and is uniformly deformed
@inline get_volume(mesh) = 4/3*π*(√sum(abs2, mesh.position[1] .- 0))^3

let
    # initialized precice and return mesh, mapping and forces
    interface = LumpedInterface(surface_mesh="../CalculiX/geom.inp")

    # iteration storage
    storage_step = []
    mmHg2Pa = 133.322387415
    
    # setup
    Plv_k = Plv_0 = 1.0         # must be nonzero
    Pact = 0.0
    
    # solver parameters
    ω⁰ = 0.180                # relaxation factor for the pressure
    ACTUATE = false           # do we provide an external pressure?
    
    # initialise
    scale_vol = 90
    EDV = 120.0                                 # ml
    ESV = scale_vol*get_volume(interface.mesh)  # ml
    @printf("Initial volume: %.2f ml\n", ESV)

    # main time loop
    while PreCICE.isCouplingOngoing()

        # the actuation pressure and the ventricular pressure must match the EDV of 120ml at sum(time) = t
        ACTUATE && (Pact = 68sum(interface.Δt))

        # read the data from the other participant, set sum(time) = t+Δt (end of this time step)
        readData!(interface)

        # volume of the sphere
        Vlv_3D = scale_vol*get_volume(interface.mesh)

        # target volume
        Vlv_0D = ESV + sum(interface.Δt)*(EDV-ESV) # target EDV is 120ml

        # fixed-point for the pressure
        Plv_k = max(Plv_0 + ω⁰*(Vlv_0D - Vlv_3D), 0.001)
        interface.step==1 && (Plv_k = 1) # first time step, used EDP to get EDV
        Plv_0 = Plv_k # for next iteration or next time step

        # we then need to recompute the forces with the correct volume and pressure
        compute_forces!(interface; p=mmHg2Pa*(Plv_k-Pact))

        # write data to the other participant, advance coupling
        writeData!(interface)

        # save the data
        push!(storage_step, [interface.step, interface.iteration, sum(interface.Δt), 
                             Plv_k, Pact, Vlv_3D, Vlv_0D, 0., 0., 0., 0.])

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            out_data = reduce(vcat,storage_step')
            CSV.write("ConstantFlow.csv", DataFrame(out_data,:auto),
                      header=["timestep","iter","time","Plv_3D", "Pact_3D", "Vlv_3D",
                              "Vlv_0D", "Pao_0D","Qao_0D", "Qmv_0D", "Plv_0D"])
        end
    end
    PreCICE.finalize()
end