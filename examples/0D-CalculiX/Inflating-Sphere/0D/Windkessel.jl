using WaterLilyPreCICE,OrdinaryDiffEq
using CSV,DataFrames,Printf

# Double Hill function inspired by Stergiopulos et al. (DOI:10.1152/ajpheart.1996.270.6.H2050)
function Elastance(t;Emin=0.0,Emax=1.0,a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    (Emax-Emin) * α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)  + Emin
end

# open-loop windkessel
function Windkessel!(du,u,p,t)
    # unpack
    (VLV, Pao, PLV) = u
    (Pfill,Rmv_fwd,Rmv_bwd,Rao_fwd,Rao_bwd,R,C)  = p

    # calculate Qmv, Pfilling>PLV; forward transmitral flow, PLV>Pfilling - backward transmitral flow
    Qmv = Pfill ≥ PLV ? (Pfill-PLV)/Rmv_fwd : (PLV-Pfill)/Rmv_bwd

    # calculate Qao, PLV>Pao; forward aortic flow, PLV>Pao; backward aortic flow
    Qao = PLV ≥ Pao ? (PLV-Pao)/Rao_fwd : (Pao-PLV)/Rao_bwd

    # rates
    du[1] = Qmv - Qao          # dVLV/dt=Qmv-Qao
    du[2] = Qao/C - Pao/(R*C)  # dPao/dt=Qao/C-Pao/RC
    du[3] = 0.0                # un-used u[3] hold the ventricular pressure
end

# compute the volume of the mesh, assumes the sphere is centered at 0 and is uniformly deformed
@inline get_volume(mesh) = 4/3*π*(√sum(abs2, mesh.position[1] .- 0))^3

# main loop
let
    # iteration storage
    storage_step = []
    Cₕ = 0.180                 # relaxation factor for the pressure
    mmHg2Pa = 133.322387415

    # Cardiac parameters for Lumped model
    Pfilling = 6.0              #mmHg; to reach EDV 120ml with sphere
    EDV = 120                   #ml; end-diastolic volume. We will use EDV with Pvenous to calculate Emin

    # Valve resistances
    Rmv_fwd = 0.01             #mmHg/ml/s; resistance in forward flow direction
    Rmv_bwd = 1e10             #mmHg/ml/s; leak resistance
    Rao_fwd = 0.01             #mmHg/ml/s; resistance in forward flow direction
    Rao_bwd = 1e10             #mmHg/ml/s; leak resistance

    # Arterial model parameters
    R_WK2 = 1                   #mmHg/ml/s
    C_WK2 = 2                   #ml/mmHg

    # setup
    PLV₁ = PLV₀ = 1.0
    P₀ = 6.01
    u₀ = [EDV, 60, P₀]           # initial conditions
    tspan = (0.0, 10.0)
    params = [Pfilling, Rmv_fwd, Rmv_bwd, Rao_fwd, Rao_bwd, R_WK2, C_WK2]

    # generate a problem to solve
    prob = ODEProblem(Windkessel!, u₀, tspan, params)

    # full control over iterations
    integrator = init(prob, Tsit5(), dtmax=0.001, reltol=1e-6, abstol=1e-9, save_everystep=false)

    # coupling interface
    interface = LumpedInterface(surface_mesh="../CalculiX/geom.inp",integrator=integrator)

    # initialise
    scale_vol = 90.0
    vol0 = scale_vol*get_volume(interface.mesh0) # convert to ml
    # print zero pressure volume
    @printf("Initial volume: %.2f ml\n", vol0)

    while PreCICE.isCouplingOngoing()

        # pressure at this step, meaning sum(interface.Δt) = t
        Pact = 68*Elastance(sum(interface.Δt))

        # read the data from the other participant, set sum(time) = t+Δt (end of this time step)
        readData!(interface)

        # solve the ODE to get VLV and Pao at t+Δt, fill the initial condition with current state
        integrate!(interface, [[interface.u₀[2], interface.u₀[3], PLV₁+Pact], interface.u₀[1]])

        # the target and current volume
        VLV_0D = interface.integrator.u[1]
        vol = scale_vol*get_volume(interface.mesh)

        # fixed-point for the pressure
        PLV₁ = max(PLV₀ + Cₕ*(VLV_0D - vol), 0.001)
        interface.step==1 && (PLV₁ = P₀) # first time step, used EDP to get EDV
        PLV₀ = PLV₁ # for next iteration or next time step

        # we then need to recompute the forces with the correct volume and pressure
        compute_forces!(interface, mmHg2Pa*PLV₁)

        # write data to the other participant, advance coupling
        writeData!(interface)

        # just to keep track of the values
        Pao,PLV = interface.integrator.u[2:3]
        Qmv = Pfilling ≥ PLV ? (Pfilling-PLV)/Rmv_fwd : (PLV-Pfilling)/Rmv_bwd
        Qao = PLV ≥ Pao ? (PLV-Pao)/Rao_fwd : (Pao-PLV)/Rao_bwd

        # save the data
        push!(storage_step, [interface.step, interface.iteration, sum(interface.Δt), PLV₁+Pact, Pact, vol,
                             VLV_0D, interface.integrator.u[2], Qao, Qmv, PLV, Pfilling])

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            out_data = reduce(vcat,storage_step')
            CSV.write("sphere_output.csv", DataFrame(out_data,:auto),
                      header=["timestep","iter","time","PLV_3D", "PACT_3D", "VLV_3D",
                              "VLV_0D", "PAO_0D","QAO_0D", "QMV_0D", "PLV_0D","Pfill_0D"])
        end
    end
    PreCICE.finalize()
end