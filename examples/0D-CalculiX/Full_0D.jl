using OrdinaryDiffEq, Plots

# ϕ actuation function
ϕᵢ(t;tC=0.10,tR=0.25,TC=0.15,TR=0.45) = 0.0<=(t-tC)%1<=TC ? 0.5*(1-cos(π*((t-tC)%1)/TC)) : (0.0<=(t-tR)%1<=TR ? 0.5*(1+cos(π*((t-tR)%1)/TR)) : 0)

# time-varying elastance function
function Eᵢ(t, E_pass, E_act_max, t_c, T_C, T_R)
    return E_pass + E_act_max * ϕᵢ(t,tC=t_c,tR=t_c+T_C,TC=T_C,TR=T_R)
end

# the minus [-] sign is p1 and p2 are alog the flow direction
@inline Rᵢ(p1, p2, r_min, r_max) = p1 ≥ p2 ? r_min : -r_max

function cardiovascular_odes!(du, u, p, t)
    VLA, VLV, VRA, VRV,
    pSYSAR, pSYSVEN, pPULAR, pPULVEN,
    QSYSAR, QSYSVEN, QPULAR, QPULVEN = u

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
    ) = p

    # Derived LV parameters (not directly in table)
    E_pass_LV = 0.05
    E_act_LV  = 2.0
    V0_LV     = 10.0

    # Elastances
    E_LA = Eᵢ(t/T_HB, E_pass_LA, E_act_LA, tC_LA, TC_LA, TR_LA)
    E_LV = Eᵢ(t/T_HB, E_pass_LV, E_act_LV, tC_RV, TC_RV, TR_RV)
    E_RA = Eᵢ(t/T_HB, E_pass_RA, E_act_RA, tC_RA, TC_RA, TR_RA)
    E_RV = Eᵢ(t/T_HB, E_pass_RV, E_act_RV, tC_RV, TC_RV, TR_RV)

    # Chamber pressures
    pLA = p_ex + E_LA * (VLA - V0_LA)
    pLV = p_ex + E_LV * (VLV - V0_LV)
    pRA = p_ex + E_RA * (VRA - V0_RA)
    pRV = p_ex + E_RV * (VRV - V0_RV)

    # Valves (forward flow only)
    QMV = (pLA - pLV) / Rᵢ(pLA, pLV, R_min, R_max)
    QAV = (pLV - pSYSAR) / Rᵢ(pLV, pSYSAR, R_min, R_max)
    QTV = (pRA - pRV) / Rᵢ(pRA, pRV, R_min, R_max)
    QPV = (pRV - pPULAR) / Rᵢ(pRV, pPULAR, R_min, R_max)

    # Volume changes
    du[1] = QPULVEN - QMV
    du[2] = QMV - QAV
    du[3] = QSYSVEN - QTV
    du[4] = QTV - QPV

    # Vascular pressure dynamics
    du[5] = (QAV - QSYSAR) / C_SYS_AR
    du[6] = (QSYSAR - QSYSVEN) / C_SYS_VEN
    du[7] = (QPV - QPULAR) / C_PUL_AR
    du[8] = (QPULAR - QPULVEN) / C_PUL_VEN

    # Flow (inertial) dynamics
    du[9]  = (-QSYSAR * R_SYS_AR - (pSYSVEN - pSYSAR)) / L_SYS_AR
    du[10] = (-QSYSVEN * R_SYS_VEN - (pRA - pSYSVEN)) / L_SYS_VEN
    du[11] = (-QPULAR * R_PUL_AR - (pPULVEN - pPULAR)) / L_PUL_AR
    du[12] = (-QPULVEN * R_PUL_VEN - (pLA - pPULVEN)) / L_PUL_VEN
end

# ---------------------------------------------------------------
# Parameter set from Table 3
# ---------------------------------------------------------------
# p = (
#     # Resistances [mmHg·s·mL⁻¹]
#     0.8, 0.1625, 0.26, 0.1625,        # R_SYS_AR, R_PUL_AR, R_SYS_VEN, R_PUL_VEN
#     # Capacitances [mL·mmHg⁻¹]
#     1.2, 10.0, 60.0, 16.0,            # C_SYS_AR, C_PUL_AR, C_SYS_VEN, C_PUL_VEN
#     # Inertances [mmHg·s²·mL⁻¹]
#     5e-3, 5e-4, 5e-4, 5e-4,           # L_SYS_AR, L_PUL_AR, L_SYS_VEN, L_PUL_VEN
#     # Passive Elastances [mmHg·mL⁻¹]
#     0.09, 0.07, 0.05,                 # E_pass_LA, E_pass_RA, E_pass_RV
#     # Active Elastances [mmHg·mL⁻¹]
#     0.07, 0.06, 0.55,                 # E_act_LA, E_act_RA, E_act_RV
#     # Reference volumes [mL]
#     4.0, 4.0, 10.0,                   # V0_LA, V0_RA, V0_RV
#     # Timing [s]
#     0.6, 0.104, 0.68,    # LA         # tC_LA, TC_LA, TR_LA
#     0.56, 0.064, 0.64,   # RA         # tC_RA, TC_RA, TR_RA
#     0.0, 0.272, 0.12,    # RV         # tC_RV, TC_RV, TR_RV
#     # Valves
#     0.0075, 75000.0,                  # R_min, R_max
#     # Heartbeat period and external pressure
#     0.8, 0.0                          # T_HB, p_ex
# )

p = (
    # Resistances [mmHg·s·mL⁻¹]
    0.416, 0.048, 0.26, 0.036,        # R_SYS_AR, R_PUL_AR, R_SYS_VEN, R_PUL_VEN
    # Capacitances [mL·mmHg⁻¹]
    1.62, 5.0, 60.0, 16.0,            # C_SYS_AR, C_PUL_AR, C_SYS_VEN, C_PUL_VEN
    # Inertances [mmHg·s²·mL⁻¹]
    5e-3, 5e-4, 5e-4, 5e-4,           # L_SYS_AR, L_PUL_AR, L_SYS_VEN, L_PUL_VEN
    # Passive Elastances [mmHg·mL⁻¹]
    0.09, 0.07, 0.05,                 # E_pass_LA, E_pass_RA, E_pass_RV
    # Active Elastances [mmHg·mL⁻¹]
    0.07, 0.06, 0.55,                 # E_act_LA, E_act_RA, E_act_RV
    # Reference volumes [mL]
    4.0, 4.0, 10.0,                   # V0_LA, V0_RA, V0_RV
    # Timing [s]
    0.6, 0.104, 0.68,    # LA         # tC_LA, TC_LA, TR_LA
    0.56, 0.064, 0.64,   # RA         # tC_RA, TC_RA, TR_RA
    0.0, 0.272, 0.12,    # RV         # tC_RV, TC_RV, TR_RV
    # Valves
    0.0075, 75000.0,                  # R_min, R_max
    # Heartbeat period and external pressure
    0.8, 0.0                          # T_HB, p_ex
)

# ---------------------------------------------------------------
# Initial conditions and simulation
# ---------------------------------------------------------------
#      VLA,  VLV,  VRA,  VRV, pSYSAR, pSYSVEN, pPULAR, pPULVEN, QSYSAR, QSYSVEN, QPULAR, QPULVEN
u0 = [80.0, 120.0, 80.0, 120.0, 90.0,     5.0,     15.0,    8.0,    0.0,    0.0,     0.0,    0.0]
tspan = (0.0, 50.0)

prob = ODEProblem(cardiovascular_odes!, u0, tspan, p)
sol = solve(prob, Tsit5(), saveat=0.01)

# plot
tvals = sol.t
plt1 = plot(tvals, sol[1:4, :]', label=["V_LA" "V_LV" "V_RA" "V_RV"], xlims=(45,50), lw=2, xlabel="Time [s]", ylabel="Volume [mL]")

pLV = Eᵢ.(tvals/0.8, 0.05, 2.0, 0.0, 0.272, 0.12) .* (sol[2, :] .- 10.0)
pRV = Eᵢ.(tvals/0.8, 0.07, 0.55, 0.0, 0.272, 0.12) .* (sol[4, :] .- 10.0)
pLA = Eᵢ.(tvals/0.8, 0.07, 0.09, 0.0, 0.272, 0.12) .* (sol[1, :] .- 4.0)
pRA = Eᵢ.(tvals/0.8, 0.07, 0.06, 0.0, 0.272, 0.12) .* (sol[3, :] .- 4.0)

plt2 = plot(tvals, [pLV, pRV, pRA, pLA], label=["P_LV" "P_RV" "P_RA" "P_LA"],
            xlims=(45,50), xlabel="Time [s]", ylabel="Pressure [mmHg]", lw=2)
plt3 = plot(tvals, sol[5:8,:]', label=["P_SYS_AR" "P_SYS_VEN" "P_PUL_AR" "P_PUL_VEN"],
            xlabel="Time [s]", ylabel="Pressure [mmHg]", lw=2)
plt4 = plot(tvals, sol[9:12,:]', label=["Q_SYS_AR" "Q_SYS_VEN" "Q_PUL_AR" "Q_PUL_VEN"],
            xlabel="Time [s]", ylabel="Flow Rate [mL/s]", lw=2)

plt5 = plot( sol[2, :], pLV, label="LV", xlabel="Volume [mL]", ylabel="Pressure [mmHg]", lw=2)
plot!(plt5,  sol[4, :], pRV, label="RV", lw=2)
plt6 = plot(sol[3, :], pRA, label="RA", xlabel="Volume [mL]", ylabel="Pressure [mmHg]", lw=2)
plot!(plt6,  sol[1, :], pLA, label="LA", lw=2)

plot(plt1, plt2, plt3, plt4, plt5, plt6, layout=(2,3), size=(900,400))