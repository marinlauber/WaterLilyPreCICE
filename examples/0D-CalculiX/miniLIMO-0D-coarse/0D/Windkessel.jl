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

# Double Hill function inspired by Stergiopulos et al. (DOI:10.1152/ajpheart.1996.270.6.H2050)
function Elastance(t;Emin=0.0,Emax=1.0,a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    (Emax-Emin) * α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)  + Emin
end

# closed-loop windkessel
function Windkessel!(du,u,p,t)
    # unpack
    (Vlv,Pa,Pv,Plv) = u
    (Ra,Ca,Rv,Cv,Rp) = p

    # flow at the two vales
    Qmv = Plv < Pv ? (Pv - Plv)/Rv : (Plv - Pv)/1e10
    Qao = Plv > Pa ? (Plv - Pa)/Ra : (Pa - Plv)/1e10

    # rates
    du[1] = Qmv - Qao                 # dVlv/dt
    du[2] = Qao/Ca + (Pv-Pa)/(Rp*Ca)  # dPa/dt
    du[3] = (Pa-Pv)/(Rp*Cv) - Qmv/Cv  # dPv/dt
    du[4] = 0.0                       # un-used u[4] hold the ventricular pressure
end

# now the pressure
function dynamic_coupling(i,t;Plv,Pact)
    i==1 && return  Plv
    i==8 && return -Plv
    i in [2,3,4] && return Plv-Pact
    i in [5,6,7] && return Pact
    i in [9,10,11] && return -(Plv-Pact)
    i in [12,13,14] && return -Pact
end

Plv_func = linear_interpolation([0.,10.,100], [0.,1.,1.])
Pact_func = linear_interpolation([0.,10.,100.], [0.,0.,20.])
activate(t;tc=10) = t < tc ? 0 : 1

# main loop
let
     # iteration storage
    storage_step = []
    Cₕ = 0.180              # relaxation factor for the pressure
    mmHg2kPa = 0.133322387415
    EDV = 120               #ml; end-diastolic volume. We will use EDV with Pvenous to calculate Emin
    Rv = 0.01               #mmHg/ml/s; resistance in forward flow direction
    Ra = 0.01               #mmHg/ml/s; resistance in forward flow direction
    Rp = 1                  #mmHg/ml/s
    Ca = 2                  #ml/mmHg
    Cv = 6.0                #

    # setup
    PLV₁ = PLV₀ = 1.0
    P₀ = 6.01
    u₀ = [EDV, 60, 6.0, P₀]           # initial conditions
    tspan = (0.0, 100.0)
    params = (Ra,Ca,Rv,Cv,Rp)

    # generate a problem to solve
    prob = ODEProblem(Windkessel!, u₀, tspan, params)

    # full control over iterations
    # integrator = init(prob, Tsit5(), dtmax=0.01, reltol=1e-6, abstol=1e-9)
    integrator = init(prob, dtmax=0.001, save_everystep=false)

    # coupling interface
    interface = LumpedInterface(surface_mesh="../CalculiX/geom.inp",func=dynamic_coupling, integrator=integrator)

    # initialise
    scale_vol = 12.5^3
    vol0 = scale_vol*WaterLilyPreCICE.volume(interface) # convert to ml
    # print zero pressure volume
    @printf("Initial volume: %.2f ml\n", vol0)

    # make the writer
    wr = vtkWriter("pouch"; attrib=custom)
    save!(wr,interface)

    # main loop
    while PreCICE.isCouplingOngoing()

        # pressure at this step, meaning sum(interface.Δt) = t
        Pact = 27.5 * activate(sum(interface.Δt)) * Elastance(sum(interface.Δt) / 10)
        # Pact = 6.01 * mmHg2kPa * Pact_func(sum(interface.Δt))

        # read the data from the other participant, set sum(time) = t+Δt (end of this time step)
        readData!(interface)

        # solve the ODE to get VLV and Pao at t+Δt, fill the initial condition with current state
        integrate!(interface, [[interface.u₀[2:end-1]..., PLV₁+Pact], interface.u₀[1]])

        # the target and current volume
        VLV_0D = interface.integrator.u[1]
        vol = scale_vol*WaterLilyPreCICE.volume(interface)
        vol = vol - (vol0 - EDV) # we are actually interested in the change from the initial volume

        # fixed-point for the pressure
        # PLV₁ = max(PLV₀ + Cₕ*(VLV_0D - vol), 0.001)
        # interface.step==1 && (PLV₁ = P₀) # first time step, used EDP to get EDV
        # PLV₀ = PLV₁ # for next iteration or next time step

        # we then need to recompute the forces with the correct volume and pressure
        PLV₁ = 6.01 * Plv_func(sum(interface.Δt))
        compute_forces!(interface;Plv=PLV₁*mmHg2kPa,Pact=Pact*mmHg2kPa)

        # write data to the other participant, advance coupling
        writeData!(interface)

        # # just to keep track of the values
        Pa,Pv,Plv = interface.integrator.u[2:4]
        Qmv = Pv ≥ Plv ? (Pv-Plv)/Rv : (Plv-Pv)/1e10
        Qao = Plv ≥ Pa ? (Plv-Pa)/Ra : (Pa-Plv)/1e10

        # save the data
        push!(storage_step, [interface.step, interface.iteration, sum(interface.Δt), PLV₁, Pact, vol,
                             VLV_0D, Pa, Qao, Qmv, Plv, Pv])

        # if we have converged, save if required
        if PreCICE.isTimeWindowComplete()
            (length(interface.Δt)+1)%50==0 && save!(wr,interface)
            println(" Interface volume: $(round(vol,digits=3)) ml")
            out_data = reduce(vcat,storage_step')
            CSV.write("sphere_output.csv", DataFrame(out_data,:auto),
                      header=["timestep","iter","time","PLV_3D", "PACT_3D", "VLV_3D",
                              "VLV_0D", "PAO_0D","QAO_0D", "QMV_0D", "PLV_0D","Pfill_0D"])
        end
    end
    close(wr)
    PreCICE.finalize()
end