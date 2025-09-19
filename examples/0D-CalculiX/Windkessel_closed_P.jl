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

# pressure volume loop function
@inline computePLV(t,V;Emin=0.05,Emax=2,V0=20) = Elastance(t;Emin,Emax) * (V-V0)

# closed-loop windkessel in pressure
function Windkessel_3P!(du,u,p,t)
    #unpack
    (Plv,Pa,Pv) = u
    (Ra,Ca,Rv,Cv,Rp,dVdt_func,dPlvdt_func) = p

    # dVdt is a prescribed function
    dVdt = dVdt_func(t)

    # what way is the flow going in the aortic valve?
    Qa = dVdt < -eps() ? -dVdt : (Pa - Plv)/1e5 # diastole, very small flow
    Qv = dVdt >  eps() ?  dVdt : (Plv - Pv)/1e5 # systole, very small flow

    # pressure in the ventricle
    Plv = dVdt > eps() ? Pv - dVdt*Rv : (dVdt < -eps() ? -dVdt*Ra + Pa : Plv)
    abs(dVdt) > 1e-5 && (u[1] = Plv) # store to access after

    # dPlv/dt is non-zero only when we have dVdt~0
    Cp = 2.4e11
    dPdt = abs(dVdt) < 1e-5 ? (find_zero(dpdt->abs2(Cp*dVdt-dpdt), 0),-dPlvdt_func(t)) : (0,0)
    @show dPdt

    # rates
    du[1] = -dPdt[2]                     # dPlv/dt
    du[2] = Qa/Ca - (Pa-Pv)/(Rp*Ca)  # dPa/dt
    du[3] = (Pa-Pv)/(Rp*Cv) - Qv/Cv  # dVv/dt
end

# solve in volume first to get dVdt
tspan = (0,20)
u = [60, 70, 8] # initial conditions for Plv, Pa, Pv
Ra = 8e6 * 1.333e-8
Rp = 3e8 * 1.333e-8
Rv = 1e6 * 1.333e-8
Ca = 8e-9 * 1.333e8
Cv = 5e-8 * 1.333e8

# we prescribe the volume from the first solution as dVdt
dVdt(t;ϵ=1e-5) = (getindex(sol_3(t+ϵ),1) - getindex(sol_3(t-ϵ),1)) / (2ϵ)
dPlvdt(t;ϵ=1e-5) = (computePLV(t+ϵ,getindex(sol_3(t+ϵ),1))-computePLV(t-ϵ,getindex(sol_3(t-ϵ),1)))/(2ϵ)
p = (Ra,Ca,Rv,Cv,Rp,dVdt,dPlvdt)

prob = ODEProblem(Windkessel_3P!, u, tspan, p)
@time sol_3P = solve(prob, Tsit5(), dtmax=0.01)

p1=plot(sol_3P.t, sol_3P[1,:], label="P_\\ LV", linewidth=2, xaxis="Time (t/T)", yaxis="Pressure (mmHg)",
        m=:circle, ms=2)
plot!(p1,sol_3.t, sol_3[1,:], label="V_\\ LV",lw=2)
plot!(p1,sol_3P.t, sol_3P[2,:], label="P_\\ Aor", linewidth=2, ls=:dash)
plot!(p1,sol_3P.t, sol_3P[3,:], label="P_\\ Ven", linewidth=2, ls=:dot)
plot!(p1,sol_3P.t, 20.0.+dVdt.(sol_3P.t)/3, label="dVLV/dt",lw=2)
plot!(p1,sol_3.t[2:end-1], 40.0.+dPlvdt.(sol_3.t[2:end-1])/10, label="dPLV/dt",lw=2)
xlims!(p1,(10,18.5)); ylims!(p1,(0,100))
# p2=plot(getindex.(sol_3.(sol_3P.t),1), sol_3P[1,:],label=:none, xlims=(0,150), ylim=(0,150),
        # xaxis="Volume (ml)", yaxis="Pressure (mmHg)")
# plot(p1,p2;layout=(1,2),size=(800,400))
# savefig("Windkessel_in_3P.png")
