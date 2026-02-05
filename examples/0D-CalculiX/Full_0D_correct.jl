using OrdinaryDiffEq, Plots

# ϕ actuation function --> For the active part of elastance
# Single-line function definition
# Julia allows unicode in identifiers
ϕᵢ(t;tC=0.10,tR=0.25,TC=0.15,TR=0.45) = 0.0<=(t-tC)%1<=TC ? 0.5*(1-cos(π*((t-tC)%1)/TC)) : (0.0<=(t-tR)%1<=TR ? 0.5*(1+cos(π*((t-tR)%1)/TR)) : 0)  # THB here?


# time-varying elastance function
# E_pass : passive elastance
# E_act_max : max active elastance
# t_c : contraction time
# T_C : contraction duration
# T_R : relaxation duration

function Eᵢ(t, E_pass, E_act_max, t_c, T_C, T_R)
    return E_pass + E_act_max * ϕᵢ(t,tC=t_c,tR=t_c+T_C,TC=T_C,TR=T_R)
end

# the minus [-] sign is p1 and p2 are alog the flow direction
@inline Rᵢ(p1, p2, r_min, r_max) = p1 ≥ p2 ? r_min : -r_max # ? if : else  # This implementation here?




# du : preallocated vector that the function must fill with derivatives
# u : state vector
# p : tuple of parameters
# ! julia convention: means that the function will modify "du" instead of returning a new vector: Important for performance during integration
function cardiovascular_odes!(du, u, p, t)
    VLA, VLV, VRA, VRV,
    pSYSAR, pSYSVEN, pPULAR, pPULVEN,
    QSYSAR, QSYSVEN, QPULAR, QPULVEN = u  # define the state vector

    (
        R_SYS_AR, R_PUL_AR, R_SYS_VEN, R_PUL_VEN,
        C_SYS_AR, C_PUL_AR, C_SYS_VEN, C_PUL_VEN,
        L_SYS_AR, L_PUL_AR, L_SYS_VEN, L_PUL_VEN,
        E_pass_LA, E_pass_RA, E_pass_RV,
        E_act_LA, E_act_RA, E_act_RV,
        V0_LA, V0_RA, V0_RV,
        tC_LA, TC_LA, TR_LA,
        tC_RA, TC_RA, TR_RA,
        tC_RV, TC_RV, TR_RV,
        R_min, R_max,
        T_HB, p_ex
    ) = p  # define the parameter vector

    # Derived LV parameters (not directly in table) --> Kostas: Taken from "Modeling the cardiac response to hemodynamic changes associated with COVID--19: a computational study"
    E_pass_LV = 0.08   # ok
    E_act_LV  = 2.75   # ok
    V0_LV     = 5.0  # ok

    # Elastances
    E_LA = Eᵢ(t/T_HB, E_pass_LA, E_act_LA, tC_LA, TC_LA, TR_LA)
    E_LV = Eᵢ(t/T_HB, E_pass_LV, E_act_LV, tC_RV, TC_RV, TR_RV)
    E_RA = Eᵢ(t/T_HB, E_pass_RA, E_act_RA, tC_RA, TC_RA, TR_RA)
    E_RV = Eᵢ(t/T_HB, E_pass_RV, E_act_RV, tC_RV, TC_RV, TR_RV)

    # Chamber pressures
    pLA = p_ex + E_LA * (VLA - V0_LA) # p_ex = 0 ?
    pLV = p_ex + E_LV * (VLV - V0_LV)
    pRA = p_ex + E_RA * (VRA - V0_RA)
    pRV = p_ex + E_RV * (VRV - V0_RV)

    # Valves (forward flow only)
    QMV = (pLA - pLV) / Rᵢ(pLA, pLV, R_min, R_max)
    QAV = (pLV - pSYSAR) / Rᵢ(pLV, pSYSAR, R_min, R_max)
    QTV = (pRA - pRV) / Rᵢ(pRA, pRV, R_min, R_max)
    QPV = (pRV - pPULAR) / Rᵢ(pRV, pPULAR, R_min, R_max)

    # Volume changes
    du[1] = QPULVEN - QMV  # dV_LA / dt
    du[2] = QMV - QAV # dV_LV / dt
    du[3] = QSYSVEN - QTV # dV_RA / dt
    du[4] = QTV - QPV # dV_RV / dt

    # Vascular pressure dynamics
    du[5] = (QAV - QSYSAR) / C_SYS_AR  # dp_sys_ar / dt
    du[6] = (QSYSAR - QSYSVEN) / C_SYS_VEN  # dp_sys_ven / dt
    du[7] = (QPV - QPULAR) / C_PUL_AR # dp_pul_ar / dt
    du[8] = (QPULAR - QPULVEN) / C_PUL_VEN # dp_pul_ven / dt

    # Flow (inertial) dynamics
    du[9]  = (-QSYSAR * R_SYS_AR - (pSYSVEN - pSYSAR)) / L_SYS_AR # dQ_sys_ar / dt
    du[10] = (-QSYSVEN * R_SYS_VEN - (pRA - pSYSVEN)) / L_SYS_VEN # dQ_sys_ven / dt
    du[11] = (-QPULAR * R_PUL_AR - (pPULVEN - pPULAR)) / L_PUL_AR # dQ_pul_ar / dt
    du[12] = (-QPULVEN * R_PUL_VEN - (pLA - pPULVEN)) / L_PUL_VEN # dQ_pul_ven / dt
end

# ---------------------------------------------------------------
# Parameter set from Table 3
# ---------------------------------------------------------------

p = [
    # Resistances [mmHg·s·mL⁻¹]
    0.8, 0.1625, 0.26, 0.1625,        # R_SYS_AR, R_PUL_AR, R_SYS_VEN, R_PUL_VEN  OK
    # Capacitances [mL·mmHg⁻¹]
    1.2, 10.0, 130.0, 16.0,            # C_SYS_AR, C_PUL_AR, C_SYS_VEN, C_PUL_VEN  OK
    # Inertances [mmHg·s²·mL⁻¹]
    5e-3, 5e-4, 5e-4, 5e-4,           # L_SYS_AR, L_PUL_AR, L_SYS_VEN, L_PUL_VEN  OK
    # Passive Elastances [mmHg·mL⁻¹]
    0.09, 0.07, 0.05,                 # E_pass_LA, E_pass_RA, E_pass_RV    OK
    # Active Elastances [mmHg·mL⁻¹]
    0.07, 0.06, 0.55,                 # E_act_LA, E_act_RA, E_act_RV   OK
    # Reference volumes [mL]
    4.0, 4.0, 10.0,                   # V0_LA, V0_RA, V0_RV   OK
    # Timing [s]
    0.8, 0.17, 0.17,    # LA         # tC_LA, TC_LA, TR_LA  OK
    0.8, 0.17, 0.17,   # RA         # tC_RA, TC_RA, TR_RA   OK
    0.0, 0.34, 0.17,    # RV         # tC_RV, TC_RV, TR_RV  OK
    # Valves
    0.0075, 75006.2,                  # R_min, R_max   OK
    # Heartbeat period and external pressure
    1.0, 0.0                          # T_HB, p_ex   OK
]

# ---------------------------------------------------------------
# Initial conditions and simulation
# ---------------------------------------------------------------
#      VLA,  VLV,  VRA,  VRV, pSYSAR, pSYSVEN, pPULAR, pPULVEN, QSYSAR, QSYSVEN, QPULAR, QPULVEN
u0 = [65.0, 120.0, 65.0, 145.0, 80.0,     30.0,     35.0,    24.0,    0.0,    0.0,     0.0,    0.0]  # OK
tspan = (0.0, 50.0)

# callbacks change HR suddenly
condition1(u, t, integrator) = t > 25.0
affect1!(integrator) = integrator.p[1]=  0.6
cb1 = DiscreteCallback(condition1, affect1!)



prob = ODEProblem(cardiovascular_odes!, u0, tspan, p)
sol = solve(prob, Tsit5(), saveat=0.01, callback=cb1, save_everystep=false) # Tsit5 -> 5-th order RK

# plot
tvals = sol.t
println(tvals[2]-tvals[1])
println(length(tvals))

plt1 = plot(tvals, sol[1:2, :]', label=["V_LA" "V_LV"],xlims = (46, 50), lw=2, xlabel="Time [s]", ylabel="Volume [mL]")
display(plt1)

# Visualize the elastances
E_LV = Eᵢ.(tvals/1.0, 0.08, 2.75, 0.0, 0.34, 0.17)
E_RV =  Eᵢ.(tvals/1.0, 0.05, 0.55, 0.0, 0.34, 0.17)
E_LA = Eᵢ.(tvals/1.0, 0.09, 0.07, 0.8, 0.17, 0.17)
E_RA = Eᵢ.(tvals/1.0, 0.07, 0.06, 0.8, 0.17, 0.17)
# TODO: Visualize R_PUL_VEN
pLV = E_LV .* (sol[2, :] .- 5.0) # call E_i for every t in tvals and .* element wise multiplication
pRV = E_RV .* (sol[4, :] .- 10.0)
pLA = E_LA .* (sol[1, :] .- 4.0)
pRA = E_LA .* (sol[3, :] .- 4.0)


plt2 = plot(tvals, [pLV, pRV, pRA, pLA], label=["P_LV" "P_RV" "P_RA" "P_LA"],
            xlabel="Time [s]", ylabel="Pressure [mmHg]", lw=2)
display(plt2)

plt3 = plot(tvals, [E_LV, E_RV], label=["E_LV" "E_RV"],
            xlabel="Time [s]", ylabel="Elastances", lw=2)
display(plt3)

plt4 = plot(tvals, sol[5,:], label="P_SYS_AR", xlims=(45, 50),
            xlabel="Time [s]", ylabel="Pressure [mmHg]", lw=2)
plot!(plt4, tvals, pLV, ylabel="p_LV", xlims=(45, 50), xlabel="Time [s]", lw=2, label="pLV")
display(plt4)

# plt5 = plot(tvals, sol[9:12,:]', label=["Q_SYS_AR" "Q_SYS_VEN" "Q_PUL_AR" "Q_PUL_VEN"],
            # xlabel="Time [s]", ylabel="Flow Rate [mL/s]", lw=2)
# display(plt5)

plt6 = plot(sol[2, :], pLV, label="LV", xlabel="Volume [mL]", ylabel="Pressure [mmHg]", lw=2)
display(plt6)




