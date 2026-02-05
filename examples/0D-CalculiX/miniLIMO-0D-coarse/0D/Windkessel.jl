using WaterLilyPreCICE,OrdinaryDiffEq
using CSV,DataFrames,Printf,Interpolations

 # vtk attributes
vtk_srf(a::LumpedInterface) = Float32[el[1] for el in a.srf_id]
vtk_center(a::LumpedInterface) = WaterLilyPreCICE.center.(a.mesh)
vtk_normal(a::LumpedInterface) = WaterLilyPreCICE.normal.(a.mesh)
vtk_dS(a::LumpedInterface) = WaterLilyPreCICE.dS.(a.mesh)
vtk_u(a::LumpedInterface) = a.deformation
vtk_f(a::LumpedInterface) = a.forces
custom = Dict("SRF" =>vtk_srf,"center"=>vtk_center,"normal"=>vtk_normal,
              "dS"=>vtk_dS,"u"=>vtk_u,"f"=>vtk_f)

# plot!(t,Elastance.(t),label="Double hill")
ϕᵢ(t;tC=0.10,tR=0.25,TC=0.15,TR=0.45) = 0.0<=(t-tC)%1<=TC ? 0.5*(1-cos(π*((t-tC)%1)/TC)) : (0.0<=(t-tR)%1<=TR ? 0.5*(1+cos(π*((t-tR)%1)/TR)) : 0)

# open-loop windkessel
function Windkessel!(du,u,p,t)
    # unpack
    (Vlv,Pa,Plv) = u
    (Ra,Ca,Rv,Cv,Rp,Pv)  = p

    # flow at the two vales
    Qmv = Pv ≥ Plv ? (Pv - Plv)/Rv : (Plv - Pv)/1e10
    Qao = Plv ≥ Pa ? (Plv - Pa)/Ra : (Pa - Plv)/1e10

    # rates
    du[1] = Qmv - Qao             # dVlv/dt=Qmv-Qao
    du[2] = Qao/Ca - Pa/(Rp*Ca)   # dPa/dt=Qao/C-Pao/RC
    du[3] = 0.0                   # un-used u[3] hold the ventricular pressure
end

# now the pressure
function dynamic_coupling(i,t;Plv,Pact)
    i == 1 && return  Plv
    i == 8 && return -Plv
    i ∈ 2:4 && return Pact
    i ∈ 5:7 && return Plv-Pact
    i ∈ 9:11 && return -Pact
    i ∈ 12:14 && return -(Plv-Pact)
end

# main loop
let
    # iteration storage
    storage_step = []
    ω⁰ = 5.20; ωᵏ = 0.0                  # relaxation factor for the pressure
    r⁰ = 0.f0; rᵏ = 0.f0
    mmHg2kPa = 0.133322387415
    EDV =  194.2                       #ml; end-diastolic volume.
    scale = 10.9

    # Kasra's parameters
    Pa2mmHg = 0.00750062 # Pa/mmHg
    m3_to_ml = 1.0e6          # m³ to ml
    Ra = 8.0e6*Pa2mmHg/m3_to_ml     # Pa.s/m³ -> mmHg.s/ml
    Rp = 1.0e8*Pa2mmHg/m3_to_ml     # Pa.s/m³
    Rv = 5.0e5*Pa2mmHg/m3_to_ml     # Pa.s/m³
    Ca = 8.0e-9*m3_to_ml/Pa2mmHg    # m³/Pa
    Cv = 5.0e-8*m3_to_ml/Pa2mmHg    # m³/Pa not used in openloop
    Pv = 6.01                       # venous pressure in mmHg

    # setup
    Plv_k = Plv_0 = 0.25/(mmHg2kPa/scale)
    u₀ = [EDV, 80, 6.01]           # initial conditions
    tspan = (0.0, 1000.0)
    params = (Ra,Ca,Rv,Cv,Rp,Pv)

    # generate a problem to solve
    prob = ODEProblem(Windkessel!, u₀, tspan, params)

    # full control over iterations
    integrator = init(prob, Tsit5(), dtmax=0.001, reltol=1e-6, abstol=1e-9, save_everystep=false)

    # coupling interfacep
    interface = LumpedInterface(surface_mesh="../CalculiX/geom_deformed.inp",
                                func=dynamic_coupling,
                                integrator=integrator)

    # initialise
    vol0 = scale^3*WaterLilyPreCICE.volume(interface) # convert to ml
    # print zero pressure volume
    @printf("Initial volume: %.6f ml\n", vol0)

    # make the writer
    wr = vtkWriter("pouch"; attrib=custom)
    save!(wr,interface)

    # main loop
    while PreCICE.isCouplingOngoing()

        # pressure at this step, meaning sum(interface.Δt) = t
        Pact = 184.0 * ϕᵢ(sum(interface.Δt)/10;tC=0.1,tR=0.4,TC=0.3,TR=0.3)

        # read the data from the other participant, set sum(time) = t+Δt (end of this time step)
        readData!(interface)

        # solve the ODE to get VLV and Pao at t+Δt, fill the initial condition with current state
        integrate!(interface, [[interface.u₀[2:end-1]..., Plv_k], interface.u₀[1]]; Δt=interface.Δt[end]/10)

        # the target and current volume
        vlv_0D = interface.integrator.u[1]
        vol_3D = scale^3*WaterLilyPreCICE.volume(interface)

        # fixed-point for the pressure, Aitken’s ∆² acceleration
        rᵏ = vlv_0D - vol_3D
        ωᵏ = ω⁰ #r⁰≈0 ? ω⁰ :
        ωₙ = -ω⁰*(r⁰ / (rᵏ-r⁰))
        Plv_k = Plv_0 + ωᵏ*rᵏ
        println(" Relaxation factor: $(round(ωᵏ,digits=3)), Residual: $(round(rᵏ,digits=6)) ml")
        Plv_0 = Plv_k # for next iteration or next time step
        r⁰ = rᵏ       # previous residual

        # we then need to recompute the forces with the correct volume and pressure
        compute_forces!(interface;Plv=Plv_k*mmHg2kPa/scale,Pact=Pact*mmHg2kPa/scale)

        # write data to the other participant, advance coupling
        writeData!(interface)

        # just to keep track of the values
        Pa,Plv = interface.integrator.u[2:3]
        Qmv = Pv ≥ Plv ? (Pv-Plv)/Rv : (Plv-Pv)/1e10
        Qao = Plv ≥ Pa ? (Plv-Pa)/Ra : (Pa-Plv)/1e10

        # save the data
        push!(storage_step, [interface.step, interface.iteration, sum(interface.Δt)/10,
                             Plv_k, Pact, vol_3D, vlv_0D, Pa, Qao, Qmv, Plv, Pv])

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            (length(interface.Δt)+1)%10==0 && save!(wr,interface)
            println(" Ventricular pressure: $(round(Plv_k,digits=3)) mmHg, ",
                      "actuation pressure: $(round(Pact,digits=3)) mmHg, ",
                      "Interface volume: $(round(vol_3D,digits=3)) ml",
                      "target volume: $(round(vlv_0D,digits=3)) ml")
            println(" Integrator time : $(round(interface.integrator.t,digits=3)) s,",
                      " Steps: $(interface.step), Iterations: $(interface.iteration)")
            flush(stdout)
            out_data = reduce(vcat,storage_step')
            CSV.write("sphere_output.csv", DataFrame(out_data,:auto),
                      header=["timestep","iter","time","PLV_3D", "PACT_3D", "VLV_3D",
                              "VLV_0D", "PAO_0D","QAO_0D", "QMV_0D", "PLV_0D","Pfill_0D"])
        end
    end
    close(wr)
    PreCICE.finalize()
end