using WaterLilyPreCICE,OrdinaryDiffEq
using CSV,DataFrames,Printf

# Double Hill function inspired by Stergiopulos et al. (DOI:10.1152/ajpheart.1996.270.6.H2050)
function Elastance(t;Emin=0.0,Emax=1.0,a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    (Emax-Emin) * α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)  + Emin
end

# open-loop windkessel
function Windkessel!(du,u,p,t)
    # unpack
    (Vlv,Pa,Plv) = u
    (Ra,Ca,Rv,Rp,Pv)  = p

    # flow at the two vales
    Qmv = Pv ≥ Plv ? (Pv - Plv)/Rv : (Plv - Pv)/1e10
    Qao = Plv ≥ Pa ? (Plv - Pa)/Ra : (Pa - Plv)/1e10

    # rates
    du[1] = Qmv - Qao             # dVlv/dt=Qmv-Qao
    du[2] = Qao/Ca - Pa/(Rp*Ca)   # dPa/dt=Qao/C-Pao/RC
    du[3] = 0.0                   # un-used u[3] hold the ventricular pressure
end

# compute the volume of the mesh, assumes the sphere is centered at 0 and is uniformly deformed
@inline get_volume(mesh) = 4/3*π*(√sum(abs2, mesh.position[1] .- 0))^3

# main loop
let
    # iteration storage
    storage_step = []
    ω⁰ = 0.180              # relaxation factor for the pressure
    mmHg2Pa = 133.322387415
    Pv = 6.0                #mmHg; to reach EDV 120ml with sphere
    EDV = 120               #ml; end-diastolic volume. We will use EDV with Pvenous to calculate Emin
    Rv = 0.01               #mmHg/ml/s; resistance in forward flow direction
    Ra = 0.01               #mmHg/ml/s; resistance in forward flow direction
    Rp = 1                  #mmHg/ml/s
    Ca = 2                  #ml/mmHg

    # setup
    Plv_k = Plv_0 = 1.0
    Pv = 6.01
    u₀ = [EDV, 60, Pv]           # initial conditions
    tspan = (0.0, 10.0)
    params = (Ra,Ca,Rv,Rp,Pv)

    # generate a problem to solve
    prob = ODEProblem(Windkessel!, u₀, tspan, params)

    # full control over iterations
    integrator = init(prob, Tsit5(), dtmax=0.001, reltol=1e-6, abstol=1e-9, save_everystep=false)

    # coupling interface
    interface = LumpedInterface(surface_mesh="../CalculiX/geom.inp", integrator=integrator)

    # initialise
    scale_vol = 90.0
    ESV = scale_vol*get_volume(interface.mesh0) # convert to ml
    @printf("Initial volume: %.2f ml\n", ESV)

    while PreCICE.isCouplingOngoing()

        # pressure at this step, meaning sum(interface.Δt) = t
        Pact = 68*Elastance(sum(interface.Δt))

        # read the data from the other participant, set sum(time) = t+Δt (end of this time step)
        readData!(interface)

        # solve the ODE to get VLV and Pao at t+Δt, fill the initial condition with current state
        integrate!(interface, [[interface.u₀[2:end-1]..., Plv_k+Pact], interface.u₀[1]])

        # the target and current volume
        Vlv_0D = interface.integrator.u[1]
        Vlv_3D = scale_vol*get_volume(interface.mesh)

        # fixed-point for the pressure
        Plv_k = max(Plv_0 + ω⁰*(Vlv_0D - Vlv_3D), 0.001)
        interface.step==1 && (Plv_k = Pv) # first time step, used EDP to get EDV
        Plv_0 = Plv_k # for next iteration or next time step

        # we then need to recompute the forces with the correct volume and pressure
        compute_forces!(interface; p=mmHg2Pa*Plv_k)

        # write data to the other participant, advance coupling
        writeData!(interface)

        # just to keep track of the values
        Pa,Plv = interface.integrator.u[2:3]
        Qmv = Pv ≥ Plv ? (Pv-Plv)/Rv : (Plv-Pv)/1e10
        Qao = Plv ≥ Pa ? (Plv-Pa)/Ra : (Pa-Plv)/1e10

        # save the data
        push!(storage_step, [interface.step, interface.iteration, sum(interface.Δt),
                             Plv_k+Pact, Pact, Vlv_3D, Vlv_0D, Pa, Qao, Qmv, Plv, Pv])

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            out_data = reduce(vcat,storage_step')
            CSV.write("Windkessel_open.csv", DataFrame(out_data,:auto),
                      header=["timestep","iter","time","Plv_3D", "Pact_3D", "Vlv_3D",
                              "Vlv_0D", "Pao_0D","Qao_0D", "Qmv_0D", "Plv_0D","Pfill_0D"])
        end
    end
    PreCICE.finalize()
end