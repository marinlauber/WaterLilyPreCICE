using StaticArrays,Plots,OrdinaryDiffEq,Roots

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
savefig("Elastance.png")

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
tspan = (0.0, 10.10)
params = [Pfilling, Rmv_fwd, Rmv_bwd, Rao_fwd, Rao_bwd, R_WK2, C_WK2]

#Pass to solver
prob = ODEProblem(Windkessel!, u₀, tspan, params)
# https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts
@time sol = solve(prob, Tsit5(), dtmax=0.001)

# full control over iterations
integrator = init(prob, Tsit5(), dtmax=0.001, reltol=1e-6, abstol=1e-9,
                    save_everystep=false)
t_sol = []

# sol_new = DiffEqBase.build_solution(prob, integrator.alg, integrator.t,
                    #  integrator.u, retcode =:Success)

@time while integrator.t < tspan[end]
    OrdinaryDiffEq.step!(integrator, 0.1, true) # stop exactly there
    push!(t_sol, [integrator.t, integrator.u...])
    # Be cautious: one should not directly mutate the t and u fields of the integrator
    SciMLBase.set_ut!(integrator, integrator.u, integrator.t)
end
t = getindex.(t_sol, 1)

p1 = plot(sol.t, computePLV.(sol.t, sol[1,:]),label="P_\\ LV",lw=2)
# plot!(p1, sol.t, sol[1,:], label="V_\\ LV",lw=2)
# plot!(p1, t, computePLV.(t, getindex.(t_sol, 2)), label="PLV",lw=2, ls=:dash)
plot!(p1, sol, idxs=[2] ,linewidth=2, xaxis="Time (t/T)", yaxis="Pressure (mmHg)",
      label="P_\\ AO", ylims=(0,100))
# plot!(p1, t, getindex.(t_sol, 3), label="PA0",lw=2, ls=:dash)

p2 = plot(getindex.(sol.u, 1), computePLV.(sol.t, getindex.(sol.u, 1)),
     label=:none, lw=2, xlims=(0,150), ylims=(0,100), xlabel="Volume")
plot(p1,p2;layout=(1,2))
# savefig("Windkessel.png")

# # formulations in P are easier for coupled 0D-3D models
# function Windkessel_P!(du,u,p,t)
#     # unpack
#     (Plv, Pao) = u
#     (Pfill,Rmv_fwd,Rmv_bwd,Rao_fwd,Rao_bwd,R,C,dVdt_func)  = p

#     # dVdt is a prescribed function
#     dVdt = dVdt_func(t)

#     Pout = 0
#     # pressure in the ventricle
#     PLV = dVdt ≥ 0 ? Pfill - dVdt*Rmv_fwd  : -dVdt*Rao_fwd + Pao
#     u[1] = PLV

#     # what way is the flow going in the aortic valve?
#     Qao = dVdt ≤ 0 ? -dVdt : (Pao-PLV)/Rao_bwd # diastole, very small flow

#     # rates
#     du[1] = 0                        # useless
#     du[2] = Qao/C - (Pao-Pout)/(R*C) # dPao/dt
# end

# Pfilling  = 5
# tspan = (0,10)
# u = [5, 60] # initial conditions
# # we prescribe the volume from the first solution as dVdt 
# dVdt(t;ϵ=5e-2) = (getindex(sol(t+ϵ),1) - getindex(sol(t-ϵ),1)) / (2ϵ)
# p = [Pfilling, Rmv_fwd, Rmv_bwd, Rao_fwd, Rao_bwd, R_WK2, C_WK2, dVdt]

# prob = ODEProblem(Windkessel_P!, u, tspan, p)
# integrator = init(prob, Tsit5(), dtmax=0.02, reltol=1e-6, abstol=1e-9,
#                     save_everystep=false)
# t_sol = []
# @time while integrator.t < tspan[end]
#     OrdinaryDiffEq.step!(integrator, integrator.dt, false)
#     push!(t_sol, [integrator.t, integrator.u...])
#     # Be cautious: one should not directly mutate the t and u fields of the integrator
#     # SciMLBase.set_ut!(integrator, integrator.u, integrator.t)
#     SciMLBase.set_u!(integrator, integrator.u)
#     SciMLBase.set_t!(integrator, integrator.t)
# end
# t = getindex.(t_sol, 1)

# # https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts
# @time sol_t = solve(prob, Tsit5(), dtmax=0.001)
# p1=plot(sol_t.t, sol_t[1,:], label="P_\\ LV", linewidth=2, xaxis="Time (t/T)", yaxis="Pressure (mmHg)")
# plot!(p1,sol_t.t, sol_t[2,:], label="P_\\ AO", linewidth=2, xlims=(0,10), ylims=(0,100))
# plot!(p1,sol_t.t, 20.0.+dVdt.(sol_t.t)/10, label="dVLV/dt",lw=2)
# xlims!(p1,(8,10))
# p2=plot(getindex.(sol.(sol_t.t),1), sol_t[1,:],label=:none, xlims=(0,150), ylim=(0,100), xaxis="Volume (ml)", yaxis="Pressure (mmHg)")
# plot(p1,p2;layout=(1,2))
# # savefig("Windkessel_in_P.png")



# function Windkessel_3!(du,u,p,t)
#     # unpack
#     (Vlv,Pa,Pv) = u
#     (Ra,Ca,Rv,Cv,Rp) = p

#     # dVdt is a prescribed function
#     Plv = computePLV(t,Vlv;Emin=0.05,Emax=2,V0=20)

#     # flow at the two vales
#     Qa = Plv > Pa ? (Plv - Pa)/Ra : (Pa - Plv)/1e10
#     Qv = Plv < Pv ? (Pv - Plv)/Rv : (Plv - Pv)/1e10

#     # rates
#     du[1] = Qv - Qa                  # dVlv/dt    
#     du[2] = Qa/Ca + (Pv-Pa)/(Rp*Ca)  # dPa/dt
#     du[3] = (Pa-Pv)/(Rp*Cv) - Qv/Cv  # dVv/dt
# end

# tspan = (0,20)
# u = [60, 70, 8] # initial conditions for Plv, Pa, Pv
# Ra = 8e6 * 1.333e-8
# Rp = 3e8 * 1.333e-8
# Rv = 1e6 * 1.333e-8
# Ca = 8e-9 * 1.333e8
# Cv = 5e-8 * 1.333e8
# p = (Ra,Ca,Rv,Cv,Rp)

# prob = ODEProblem(Windkessel_3!, u, tspan, p)
# # https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts
# @time sol_3 = solve(prob, Tsit5(), dtmax=1e-4)
# p1=plot(sol_3.t, computePLV.(sol_3.t, sol_3[1,:]), label="P_\\ LV", lw=2, xaxis="Time (t/T)", yaxis="Pressure (mmHg)")
# plot!(p1,sol_3.t, sol_3[1,:],label="V_\\ LV", lw=2)
# plot!(p1,sol_3.t, sol_3[2,:], label="P_\\ Aor", lw=2, ls=:dash) #, xlims=(0,10), ylims=(0,100))
# plot!(p1,sol_3.t, sol_3[3,:], label="P_\\ Ven", lw=2, ls=:dot) #, xlims=(0,10), ylims=(0,100))
# xlims!(p1,(18,20.5)); ylims!(p1,(0,100))
# p2=plot(sol_3[1,:], computePLV.(sol_3.t, sol_3[1,:]),label=:none, xlims=(0,150), ylim=(0,150),
#         xaxis="Volume (ml)", yaxis="Pressure (mmHg)")
# plot(p1,p2;layout=(1,2),size=(800,400))
# # savefig("Windkessel_in_3.png")


# function Windkessel_3P!(du,u,p,t)
#     #unpack
#     (Plv,Pa,Pv) = u
#     (Ra,Ca,Rv,Cv,Rp,dVdt_func,dPlvdt_func) = p

#     # dVdt is a prescribed function
#     dVdt = dVdt_func(t)
    
#     # what way is the flow going in the aortic valve?
#     Qa = dVdt < -eps() ? -dVdt : (Pa - Plv)/1e5 # diastole, very small flow
#     Qv = dVdt >  eps() ?  dVdt : (Plv - Pv)/1e5 # systole, very small flow

#     # pressure in the ventricle
#     Plv = dVdt > eps() ? Pv - dVdt*Rv : (dVdt < -eps() ? -dVdt*Ra + Pa : Plv)
#     abs(dVdt) > 1e-5 && (u[1] = Plv) # store to access after
    
#     # dPlv/dt is non-zero only when we have dVdt~0
#     Cp = 2.4e11
#     dPdt = abs(dVdt) < 1e-5 ? (find_zero(dpdt->abs2(Cp*dVdt-dpdt), 0),-dPlvdt_func(t)) : (0,0)
#     @show dPdt

#     # rates
#     du[1] = -dPdt[2]                     # dPlv/dt  
#     du[2] = Qa/Ca - (Pa-Pv)/(Rp*Ca)  # dPa/dt
#     du[3] = (Pa-Pv)/(Rp*Cv) - Qv/Cv  # dVv/dt  
# end

# using ForwardDiff
# tspan = (0,18)
# # we prescribe the volume from the first solution as dVdt 
# dVdt(t;ϵ=1e-5) = (getindex(sol_3(t+ϵ),1) - getindex(sol_3(t-ϵ),1)) / (2ϵ)
# dPlvdt(t;ϵ=1e-5) = (computePLV(t+ϵ,getindex(sol_3(t+ϵ),1))-computePLV(t-ϵ,getindex(sol_3(t-ϵ),1)))/(2ϵ)

# # dVdt(t) = ForwardDiff.derivative(τ->sol_3(τ), t)
# p = (Ra,Ca,Rv,Cv,Rp,dVdt,dPlvdt)

# prob = ODEProblem(Windkessel_3P!, u, tspan, p)
# @time sol_3P = solve(prob, Tsit5(), dtmax=0.01)

# p1=plot(sol_3P.t, sol_3P[1,:], label="P_\\ LV", linewidth=2, xaxis="Time (t/T)", yaxis="Pressure (mmHg)",
#         m=:circle, ms=2)
# plot!(p1,sol_3.t, sol_3[1,:], label="V_\\ LV",lw=2)
# plot!(p1,sol_3P.t, sol_3P[2,:], label="P_\\ Aor", linewidth=2, ls=:dash)
# plot!(p1,sol_3P.t, sol_3P[3,:], label="P_\\ Ven", linewidth=2, ls=:dot)
# plot!(p1,sol_3P.t, 20.0.+dVdt.(sol_3P.t)/3, label="dVLV/dt",lw=2)
# plot!(p1,sol_3.t[2:end-1], 40.0.+dPlvdt.(sol_3.t[2:end-1])/10, label="dPLV/dt",lw=2)
# xlims!(p1,(10,18.5)); ylims!(p1,(0,100))
# # p2=plot(getindex.(sol_3.(sol_3P.t),1), sol_3P[1,:],label=:none, xlims=(0,150), ylim=(0,150),
#         # xaxis="Volume (ml)", yaxis="Pressure (mmHg)")
# # plot(p1,p2;layout=(1,2),size=(800,400))
# # savefig("Windkessel_in_3P.png")
