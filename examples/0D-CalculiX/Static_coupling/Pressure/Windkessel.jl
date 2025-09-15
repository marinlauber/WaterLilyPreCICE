using GeometryBasics,WriteVTK,StaticArrays,PreCICE
using Printf,CSV,DataFrames
using LinearAlgebra: cross
using OrdinaryDiffEq,Roots

# reading inp files
include("IO.jl")
include("utils.jl")

# Double Hill function inspired by Stergiopulos et al. (DOI:10.1152/ajpheart.1996.270.6.H2050)
function Elastance(t;Emin=0.05,Emax=2,a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    (Emax-Emin) * α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)  + Emin
end # amplitude is 2
# times = collect(0:0.05:1) .+ 0.65
# plot(times, Elastance.(times), xlabel="time (s)", ylabel="Elastance (mmHg/ml)", label="E(t)")

# pressure volume loop function
@inline computePLV(t,V;Emin=0.05,Emax=2,V0=20) = Elastance(t;Emin,Emax) * (V-V0)

# open-loop windkessel
function Windkessel!(du,u,p,t)
    # unpack
    (VLV, Pao, PLV) = u
    (Pfill,Rmv_fwd,Rmv_bwd,Rao_fwd,Rao_bwd,R,C)  = p

    # first calculate PLV from elastance and VLV
    # PLV = computePLV(t,VLV)

    # calculate Qmv, Pfilling>PLV; forward transmitral flow, PLV>Pfilling - backward transmitral flow
    Qmv = Pfill ≥ PLV ? (Pfill-PLV)/Rmv_fwd : (PLV-Pfill)/Rmv_bwd

    # calculate Qao, PLV>Pao; forward aortic flow, PLV>Pao; backward aortic flow
    Qao = PLV ≥ Pao ? (PLV-Pao)/Rao_fwd : (Pao-PLV)/Rao_bwd

    # rates
    du[1] = Qmv - Qao          #dVLV/dt=Qmv-Qao
    du[2] = Qao/C - Pao/(R*C)  #dPao/dt=Qao/C-Pao/RC
    du[3] = 0.0                # un-used u[3] hold the ventricular pressure
end

# solver the problem
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
    vol0 = 100*get_volume(mesh0) # convert to ml

    # iteration storage
    storage_step,pressure_iter,volume_iter = [],[],[]
    iteration = step = 1
    Q_target = 20.0  # ml/s between t=0 and t=1
    Cₕ = 0.64    # relaxation factor for the pressure
    simple = true    # use a simple fixed-point or a secant method

    # set model parameters
    # Cardiac parameters
    Emax = 2                    #mmHg/ml; slope of the ESPVR
    V0 = 20                     #ml; intercept with volume axis of the ESPVR
    Pfilling = 5                #mmHg; venous filling pressure
    EDV = 120                   #ml; end-diastolic volume. We will use EDV with Pvenous to calculate Emin
    Emin = Pfilling/(EDV-V0)    #mmHg/ml
    HR = 60                     #heart rate in beats/min

    # Valve resistances
    Rmv_fwd = 0.002             #mmHg/ml/s; resistance in forward flow direction
    Rmv_bwd = 1e10              #mmHg/ml/s; leak resistance
    Rao_fwd = 0.002             #mmHg/ml/s; resistance in forward flow direction
    Rao_bwd = 1e10              #mmHg/ml/s; leak resistance

    # Arterial model parameters
    R_WK2 = 1                   #mmHg/ml/s
    C_WK2 = 2                   #ml/mmHg

    #Setup
    PLV₁ = PLV₀ = Pfilling
    u₀ = [EDV, 40, PLV₁] # initial conditions
    tspan = (0.0, 10.0)
    params = [Pfilling, Rmv_fwd, Rmv_bwd, Rao_fwd, Rao_bwd, R_WK2, C_WK2]

    # generate a problem to solve
    prob = ODEProblem(Windkessel!, u₀, tspan, params)

    # full control over iterations
    integrator = init(prob, Tsit5(), dtmax=0.001, reltol=1e-6, abstol=1e-9,
                      save_everystep=false)
    t_sol = []
    # initial storage
    ut_store = [integrator.t, integrator.u...]

    # main time loop
    while PreCICE.isCouplingOngoing()

        if PreCICE.requiresWritingCheckpoint()
            # save state at sum(time) = t
            ut_store = [integrator.t, integrator.u...]
        end

        # set time step
        precice_dt = PreCICE.getMaxTimeStepSize()
        dt = min(precice_dt, solver_dt)
        push!(time, dt) # sum(time) = t+Δt (end of this time step)

        # read data from other participant
        displacements = PreCICE.readData("Solid-Mesh", "Displacements", ControlPointsID, dt)
        # update the mesh
        mesh = update_mesh(mesh0, displacements)
        vol = 100*get_volume(mesh) # scale to ml

        # volume of sphere
        push!(volume_iter, vol)

        # solve the ODE to get VLV and Pao at t+Δt
        OrdinaryDiffEq.step!(integrator, dt, true) # stop exactly at t+Δt
        push!(t_sol, [integrator.t, integrator.u[1:2]...])

        # target volume
        target_volume = vol0 + Q_target*sum(time) # target flow rate is Q_target
        # target_volume = integrator.u[1] # volume from 0D

        # just to keep track of the values
        # PLV = computePLV(sum(time),target_volume)
        Pao = integrator.u[2]
        PLV = integrator.u[3]
        Qmv = Pfilling ≥ PLV ? (Pfilling-PLV)/Rmv_fwd : (PLV-Pfilling)/Rmv_bwd
        Qao = PLV ≥ Pao ? (PLV-Pao)/Rao_fwd : (Pao-PLV)/Rao_bwd

        # fixed-point for the pressure
        if simple==true
            PLV₁ = PLV₀ + Cₕ*(target_volume - volume_iter[end])
        else
            ∂p = pressure_iter[end] - pressure_iter[end-1]
            ∂v = volume_iter[end] - volume_iter[end-1]
            PLV₁ = PLV₀ + (∂p/∂v)*(target_volume - volume[end])
        end
        PLV₀ = PLV₁

        # actuation pressure at this time
        #@TODO check the time scaling
        Pact = 0*Elastance(sum(time); Emin=Emin, Emax=Emax)
        push!(pressure_iter, PLV₁-Pact)
        if PLV₁ < Pact
            @warn "The pressure in the ventricle must be higher than the actuation pressure"
        end
        
        # we then need to recompute the forces with the correct volume and pressure
        forces .= 0 # reset the forces
        for (i,id) in srf_id
            f = dS(@views(mesh[id])) .* max(0,(PLV₁-Pact)) # pressure jump
            forces[map_id[id],:] .+= transpose(f)./3 # add all the contribution from the faces to the nodes
        end

        # write the force at the nodes
        PreCICE.writeData("Solid-Mesh", "Forces", ControlPointsID, forces)

        # do the coupling
        PreCICE.advance(dt) # advance to t+Δt

        # save the data
        push!(storage_step, [step, iteration, sum(time), PLV₁, Pact, vol,
                             integrator.u[1:2]..., Qao, Qmv, PLV])

        # read checkpoint if required or move on
        if PreCICE.requiresReadingCheckpoint()
            # revert to sum(time) = t
            pop!(time) # remove last time step since we are going back
            # revert integrator, put the correct pressure there
            SciMLBase.set_ut!(integrator, [ut_store[2], ut_store[3], PLV₁], ut_store[1])
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
    PreCICE.finalize()
end