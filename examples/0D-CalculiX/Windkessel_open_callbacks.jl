using StaticArrays,Plots,OrdinaryDiffEq,Roots

# set model parameters
# Cardiac parameters
Emax = 2                    #mmHg/ml; slope of the ESPVR
V0 = 40                     #ml; intercept with volume axis of the ESPVR
Pfilling = 6.5                #mmHg; venous filling pressure
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

# pressure volume loop function
@inline computePLV(t,V;Emin=0.05,Emax=2,V0=20) = Elastance(t;Emin,Emax) * (V-V0)

# open-loop windkessel
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
tspan = (0.0, 15.0)
params = [Pfilling, Rmv_fwd, Rmv_bwd, Rao_fwd, Rao_bwd, R_WK2, C_WK2]

# callbackschange the filling pressure sudenly
condition(u, t, integrator) = t > 4
affect!(integrator) = integrator.p[1]= 2*Pfilling  # change filling pressure after 4 seconds
cb = DiscreteCallback(condition, affect!)

#Pass to solver
prob = ODEProblem(Windkessel!, u₀, tspan, params)
# https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts
# @time sol = solve(prob, Tsit5(), dtmax=0.001, callback = cb)

# full control over iterations
integrator = init(ODEProblem(Windkessel!, u₀, tspan, params), Tsit5(),
                  dtmax=0.001, reltol=1e-6, abstol=1e-9,
                  save_everystep=false, callback=cb)
t_sol = []

# integrate until the end of the time span
@time while integrator.t < tspan[end]
    OrdinaryDiffEq.step!(integrator, 0.01, true) # stop exactly there
    push!(t_sol, [integrator.t, integrator.u...])
    # Be cautious: one should not directly mutate the t and u fields of the integrator
    SciMLBase.set_ut!(integrator, integrator.u, integrator.t)
end
t = getindex.(t_sol, 1)

# plot
p1 = plot(sol.t, computePLV.(sol.t, sol[1,:]),label="P_\\ LV",lw=2)
plot!(p1, t, computePLV.(t, getindex.(t_sol, 2)), label="P_\\ LV",lw=2,ls=:dashdot)
plot!(p1, sol, idxs=[2] ,linewidth=2, xaxis="Time (t/T)", yaxis="Pressure (mmHg)",
      label="P_\\ AO", ylims=(0,200))
plot!(p1, sol.t, Pfilling.*ones(length(sol.t)), label="P_\\ Fill", lw=2)
plot!(p1, t, getindex.(t_sol, 3), label="P_\\ AO",lw=2,ls=:dashdot)
p2 = plot(getindex.(sol.u, 1), computePLV.(sol.t, getindex.(sol.u, 1)),
     label=:none, lw=2, xlims=(0,200), ylims=(0,200), xlabel="Volume")
p3 = plot(sol.t, sol[1,:],label="V_\\ LV",lw=2)
Pao = sol[2,:]
PLV = computePLV.(sol.t, sol[1,:])
Qmv = map(p->Pfilling ≥ p ? (Pfilling-p)/Rmv_fwd : (p-Pfilling)/Rmv_bwd, PLV)
Qao = map(T->T[1]≥T[2] ? (T[1]-T[2])/Rao_fwd : (T[2]-T[1])/Rao_bwd, zip(PLV,Pao))
plot!(p3, t, getindex.(t_sol, 2),label="V_\\ LV",lw=2,ls=:dashdot)
p3 = plot(sol.t, Qmv, label="Q_\\ mv", lw=2, xaxis="Time (t/T)", yaxis="Flow rate (ml/s)")
plot!(p3, sol.t, Qao, label="Q_\\ ao", lw=2, ylims=(0,1e3))
# bunch together
plot(p1,p2,p3;layout=(1,3),size=(1200,400))
# savefig("windkessel_0D_open_loop_callbacks.png")
