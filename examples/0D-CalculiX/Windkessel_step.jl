using StaticArrays,Plots,OrdinaryDiffEq

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

# Double Hill function inspired by Stergiopulos et al. (DOI:10.1152/ajpheart.1996.270.6.H2050)
function Elastance(t;Emin=0.05,Emax=2,a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    (Emax-Emin) * α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)  + Emin
end

# plot to check
t = collect(0:1/HR:4)
plot(t, Elastance.(t), label="Elastance",lw=2, xlabel="Time (t/T)", ylabel="Elastance (mmHg/ml)", legend=:topleft)

# pressure volume loop function
@inline computePLV(t,V;Emin=0.05,Emax=2,V0=20) = Elastance(t;Emin,Emax) * (V-V0)
plot!(t,computePLV.(t,EDV),label="PLV",lw=2)

function Windkessel!(du,u,p,t)
    # unpack
    (VLV, Pao) = u
    (Pfill,Rmv_fwd,Rmv_bwd,Rao_fwd,Rao_bwd,R,C)  = p

    #first calculate PLV from elastance and VLV 
    PLV = computePLV(t,VLV)

    # calculate Qmv, Pfilling>PLV; forward transmitral flow, PLV>Pfilling - backward transmitral flow
    Qmv = Pfill ≥ PLV ? (Pfill-PLV)/Rmv_fwd : (PLV-Pfill)/Rmv_bwd
    
    #calculate Qao, PLV>Pao; forward aortic flow, PLV>Pao; backward aortic flow
    Qao = PLV ≥ Pao ? (PLV-Pao)/Rao_fwd : (Pao-PLV)/Rao_bwd

    # rates
    du[1] = Qmv - Qao          #dVLV/dt=Qmv-Qao
    du[2] = Qao/C - Pao/(R*C)  #dPao/dt=Qao/C-Pao/RC
end

#Setup
u₀ = [EDV, 60] # initial conditions
tspan = (0.0, 4.0)
params = [Pfilling, Rmv_fwd, Rmv_bwd, Rao_fwd, Rao_bwd, R_WK2, C_WK2]

#Pass to solver
prob = ODEProblem(Windkessel!, u₀, tspan, params)
# https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts
@time sol = solve(prob, Tsit5(), dtmax=0.02)

# full control over iterations
integrator = init(prob, Tsit5(), dtmax=0.02, reltol=1e-6, abstol=1e-9,
                    save_everystep=false)
t_sol = []

# sol_new = DiffEqBase.build_solution(prob, integrator.alg, integrator.t,
                    #  integrator.u, retcode =:Success)

@time while integrator.t < tspan[end]
    OrdinaryDiffEq.step!(integrator, integrator.dt, false)
    push!(t_sol, [integrator.t, integrator.u...])
    # Be cautious: one should not directly mutate the t and u fields of the integrator
    SciMLBase.set_ut!(integrator, integrator.u, integrator.t)
end
t = getindex.(t_sol, 1)

p1 = plot(sol.t, computePLV.(sol.t, getindex.(sol.u, 1)),label="P_\\ LV",lw=2)
plot!(p1, t, computePLV.(t, getindex.(t_sol, 2)), label="PLV",lw=2, ls=:dash)
plot!(p1, sol, idxs=[2] ,linewidth=2, xaxis="Time (t/T)", yaxis="Pressure (mmHg)", label="P_\\ AO")
plot!(p1, t, getindex.(t_sol, 3), label="PA0",lw=2, ls=:dash)

p2 = plot(getindex.(sol.u, 1), computePLV.(sol.t, getindex.(sol.u, 1)),
     label=:none, lw=2, xlims=(0,150), ylims=(0,100), xlabel="Volume")
plot(p1,p2;layout=(1,2))


# formulations in P are easier for coupled 0D-3D models
function Windkessel_P!(du,u,p,t)
    # unpack
    (PLV, Qoa) = u
    (Pfill,Rmv_fwd,Rmv_bwd,Rao_fwd,Rao_bwd,R,C)  = p
    @show Pfill,PLV,Pao

    # what the flow rate around the mitral valve 
    Qmv = Pfill ≥ PLV ? (Pfill-PLV)/Rmv_fwd : (PLV-Pfill)/Rmv_bwd
    @show Qmv

    # what way is the flow going in the aortic valve?
    # Qao = PLV ≥ Pao ? (PLV-Pao)/Rao_fwd : (Pao-PLV)/Rao_bwd
    # @show Qao

    # given the flow rate into the aorta, what is the ventricular pressure?
    PLV = PLV ≥ Pao ? -Qao*Rao_fwd + Pao : Qmv*Rmv_fwd + Pfill
    # u[1] = PLV
    @show PLV

    # rates
    # du[1] = 0                  #dPLV/dt=PLV-Pao
    # du[2] = Qao/C - Pao/(R*C)  #dPao/dt=Qao/C-Pao/RC

    # rates
    du[1] = 0 # dPLV/dt
    du[2] = 0 # dQao/dt
end

u = [2Pfilling, Pfilling+1] # initial conditions
du = zeros(2)
p = [Pfilling, Rmv_fwd, Rmv_bwd, Rao_fwd, Rao_bwd, R_WK2, C_WK2]
for i in 1:8
    Windkessel_P!(du,u,p,0.0)
end