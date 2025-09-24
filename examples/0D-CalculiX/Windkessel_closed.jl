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
Rao_fwd = 0.002             #mmHg/ml/s; resistance in forward flow direction

# Arterial model parameters
R_WK2 = 1                   #mmHg/ml/s
C_WK2 = 2                   #ml/mmHg

# Double Hill function inspired by Stergiopulos et al. (DOI:10.1152/ajpheart.1996.270.6.H2050)
function Elastance(t;Emin=0.05,Emax=2,a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    (Emax-Emin) * α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)  + Emin
end

# pressure volume loop function
@inline computePLV(t,V;Emin=0.05,Emax=2,V0=20) = Elastance(t;Emin,Emax) * (V-V0)

# closed-loop windkessel
function Windkessel_3!(du,u,p,t)
    # unpack
    (Vlv,Pa,Pv) = u
    (Ra,Ca,Rv,Cv,Rp) = p

    # dVdt is a prescribed function
    Plv = computePLV(t,Vlv;Emin=0.05,Emax=2,V0=20)

    # flow at the two vales
    Qmv = Plv < Pv ? (Pv - Plv)/Rv : (Plv - Pv)/1e10
    Qao = Plv > Pa ? (Plv - Pa)/Ra : (Pa - Plv)/1e10

    # rates
    du[1] = Qmv - Qao                  # dVlv/dt
    du[2] = Qao/Ca + (Pv-Pa)/(Rp*Ca)  # dPa/dt
    du[3] = (Pa-Pv)/(Rp*Cv) - Qmv/Cv  # dPv/dt
end

tspan = (0,20)
u = [60, 70, 8] # initial conditions for Plv, Pa, Pv
Ra = 8e6 * 1.333e-8
Rp = 3e8 * 1.333e-8
Rv = 1e6 * 1.333e-8
Ca = 8e-9 * 1.333e8
Cv = 5e-8 * 1.333e8
p = (Ra,Ca,Rv,Cv,Rp)

prob = ODEProblem(Windkessel_3!, u, tspan, p)
# https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts
@time sol_3 = solve(prob, Tsit5(), dtmax=1e-4)
p1=plot(sol_3.t, computePLV.(sol_3.t, sol_3[1,:]), label="P_\\ LV", lw=2, xaxis="Time (t/T)", yaxis="Pressure (mmHg)")
plot!(p1,sol_3.t, sol_3[1,:],label="V_\\ LV", lw=2)
plot!(p1,sol_3.t, sol_3[2,:], label="P_\\ Aor", lw=2, ls=:dash) #, xlims=(0,10), ylims=(0,100))
plot!(p1,sol_3.t, sol_3[3,:], label="P_\\ Ven", lw=2, ls=:dot) #, xlims=(0,10), ylims=(0,100))
xlims!(p1,(0,20.5)); ylims!(p1,(0,100))
p2=plot(sol_3[1,:], computePLV.(sol_3.t, sol_3[1,:]),label=:none, xlims=(0,150), ylim=(0,150),
        xaxis="Volume (ml)", yaxis="Pressure (mmHg)")
plot(p1,p2;layout=(1,2),size=(800,400))
# savefig("Windkessel_in_3.png")