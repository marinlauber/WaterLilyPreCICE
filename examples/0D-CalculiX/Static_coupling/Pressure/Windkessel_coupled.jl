using GeometryBasics,WriteVTK,StaticArrays,PreCICE
using Printf,CSV,DataFrames
using LinearAlgebra: cross
using OrdinaryDiffEq,Roots

# reading inp files
include("IO.jl")
include("utils.jl")

# Double Hill function inspired by Stergiopulos et al. (DOI:10.1152/ajpheart.1996.270.6.H2050)
function Elastance(t;Emin=0.0,Emax=1,a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    (Emax-Emin) * α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)  + Emin
end # amplitude is 2
# plot(collect(0:0.005:1), Elastance.(collect(0:0.005:1)))

# pressure volume loop function
@inline computePLV(t,V;Emin=0.05,Emax=2,V0=20) = Elastance(t;Emin,Emax) * (V-V0)

# open-loop windkessel
function Windkessel!(du,u,p,t)
    # unpack
    (VLV, Pao, PLV) = u
    (Pfill,Rmv_fwd,Rmv_bwd,Rao_fwd,Rao_bwd,R,C)  = p

    # first calculate PLV from elastance and VLV
    PLV = computePLV(t,VLV)

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
    scale_vol = 80.0
    vol0 = scale_vol*get_volume(mesh0) # convert to ml

    # iteration storage
    storage_step = []
    iteration = step = 1
    Cₕ = 0.1    # relaxation factor for the pressure
    mmHg2Pa = 133.322387415

    # Cardiac parameters
    Pfilling = 5.0              #mmHg; @TODO to reach EDV 120ml with sphere
    EDV = 120                   #ml; end-diastolic volume. We will use EDV with Pvenous to calculate Emin

    # Valve resistances
    Rmv_fwd = 0.002             #mmHg/ml/s; resistance in forward flow direction
    Rmv_bwd = 1e10              #mmHg/ml/s; leak resistance
    Rao_fwd = 0.002             #mmHg/ml/s; resistance in forward flow direction
    Rao_bwd = 1e10              #mmHg/ml/s; leak resistance

    # Arterial model parameters
    R_WK2 = 1                   #mmHg/ml/s
    C_WK2 = 2                   #ml/mmHg

    #Setup
    PLV₁ = PLV₀ = 0.0
    Pact = 0.0
    u₀ = [EDV, 60, Pfilling] # initial conditions
    tspan = (0.0, 10.0)
    params = [Pfilling, Rmv_fwd, Rmv_bwd, Rao_fwd, Rao_bwd, R_WK2, C_WK2]

    # generate a problem to solve
    prob = ODEProblem(Windkessel!, u₀, tspan, params)

    # full control over iterations
    integrator = init(prob, Tsit5(), dtmax=0.001, reltol=1e-6, abstol=1e-9, save_everystep=false)
    # initial storage
    ut_store = [integrator.t, integrator.u...]

    # helper to store the previous forces to bypass convergence
    force_previous = copy(forces)

    # main time loop
    while PreCICE.isCouplingOngoing()

        if PreCICE.requiresWritingCheckpoint()
            # save state at sum(time) = t
            ut_store = [integrator.t, integrator.u...]
        end

        # actuation pressure at this time
        # (sum(time)>-1) && (Pact = 70*Elastance(sum(time))) # this is a time =t, not t+Δt

        # set time step
        precice_dt = PreCICE.getMaxTimeStepSize()
        dt = min(precice_dt, solver_dt)
        push!(time, dt) # sum(time) = t+Δt (end of this time step)

        # read data from other participant
        displacements = PreCICE.readData("Solid-Mesh", "Displacements", ControlPointsID, dt)
        # update the mesh
        mesh = update_mesh(mesh0, displacements)
        vol = scale_vol*get_volume(mesh) # scale to ml

        # the first time-step is used to reach the EDPVR of 120ml@18mmHg
        # if sum(time)<-2
            # VLV_0D = vol0 + sum(time)*(EDV-vol0)
        # else
            # solve the ODE to get VLV and Pao at t+Δt
        OrdinaryDiffEq.step!(integrator, dt, true) # stop exactly at t+Δt
        VLV_0D = integrator.u[1] # volume target from 0D
        # end

        # fixed-point for the pressure to match volume
        PLV₁ = PLV₀ + Cₕ*(VLV_0D - vol)

        # we then need to recompute the forces with the correct volume and pressure
        compute_forces!(forces, mmHg2Pa*PLV₁, mesh, srf_id, map_id)
        # if we are converged in pressure, just exit the coupling loop
        # if abs(PLV₀-PLV₁)/PLV₁ < 1e-3
            # forces .= force_previous # use previous forces if not converged
        # end # TODO cannot use this as it break for iso-volumetric
        PLV₀ = PLV₁ # for next iteration or next time step

        # write the force at the nodes
        PreCICE.writeData("Solid-Mesh", "Forces", ControlPointsID, forces)

        # do the coupling
        PreCICE.advance(dt) # advance to t+Δt

        # just to keep track of the values
        Pao,PLV = integrator.u[2:3]
        Qmv = Pfilling ≥ PLV ? (Pfilling-PLV)/Rmv_fwd : (PLV-Pfilling)/Rmv_bwd
        Qao = PLV ≥ Pao ? (PLV-Pao)/Rao_fwd : (Pao-PLV)/Rao_bwd
        # save the data
        push!(storage_step, [step, iteration, sum(time), PLV₁, Pact, vol,
                             VLV_0D, integrator.u[2], Qao, Qmv, PLV])

        # read checkpoint if required or move on
        if PreCICE.requiresReadingCheckpoint() # revert to sum(time) = t
            pop!(time) # remove last time step since we are going back
            # revert integrator, put the correct pressure there
            # SciMLBase.set_ut!(integrator, [ut_store[2], ut_store[3], PLV₁], ut_store[1])
            SciMLBase.set_ut!(integrator, [ut_store[2], ut_store[3], ut_store[4]], ut_store[1])
            iteration += 1
            force_previous .= forces # previous forces
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