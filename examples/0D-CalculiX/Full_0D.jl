using OrdinaryDiffEq, Plots

# ---------------------------------------------------------------
# Activation transient (Eq. 8)
# ---------------------------------------------------------------
function phi_activation(t, t_c, t_R, T_C, T_R, T_HB)
    r = mod(t - t_c, T_HB)
    if 0.0 <= r < T_C
        return 0.5 * (1 - cos(pi * r / T_C))
    end
    s = mod(t - t_R, T_HB)
    if 0.0 <= s < T_R
        return 0.5 * (1 + cos(pi * s / T_R))
    end
    return 0.0
end

# ---------------------------------------------------------------
# Time-varying elastance
# ---------------------------------------------------------------
function elastance(t, E_pass, E_act_max, t_c, t_R, T_C, T_R, T_HB)
    return E_pass + E_act_max * phi_activation(t, t_c, t_R, T_C, T_R, T_HB)
end

# ---------------------------------------------------------------
# Cardiovascular ODE system
# ---------------------------------------------------------------
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
        tC_LA, TR_LA, tR_LA,
        tC_RA, TR_RA, tR_RA,
        tC_RV, TR_RV, tR_RV,
        R_min, R_max,
        T_HB, p_ex
    ) = p

    # Derived LV parameters (not directly in table)
    E_pass_LV = 0.05
    E_act_LV  = 2.0
    V0_LV     = 10.0

    # Elastances
    E_LA = elastance(t, E_pass_LA, E_act_LA, tC_LA, tR_LA, TR_LA, TR_LA, T_HB)
    E_LV = elastance(t, E_pass_LV, E_act_LV, tC_RV, tR_RV, TR_RV, TR_RV, T_HB)
    E_RA = elastance(t, E_pass_RA, E_act_RA, tC_RA, tR_RA, TR_RA, TR_RA, T_HB)
    E_RV = elastance(t, E_pass_RV, E_act_RV, tC_RV, tR_RV, TR_RV, TR_RV, T_HB)

    # Chamber pressures
    pLA = p_ex + E_LA * (VLA - V0_LA)
    pLV = p_ex + E_LV * (VLV - V0_LV)
    pRA = p_ex + E_RA * (VRA - V0_RA)
    pRV = p_ex + E_RV * (VRV - V0_RV)

    # Valves (forward flow only)
    QMV = max((pLA - pLV) / R_min, 0.0)
    QAV = max((pLV - pSYSAR) / R_min, 0.0)
    QTV = max((pRA - pRV) / R_min, 0.0)
    QPV = max((pRV - pPULAR) / R_min, 0.0)

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

    return nothing
end

# ---------------------------------------------------------------
# Parameter set from Table 3
# ---------------------------------------------------------------
p = (
    # Resistances [mmHg·s·mL⁻¹]
    0.8, 0.1625, 0.26, 0.1625,
    # Capacitances [mL·mmHg⁻¹]
    1.2, 10.0, 60.0, 16.0,
    # Inertances [mmHg·s²·mL⁻¹]
    5e-3, 5e-4, 5e-4, 5e-4,
    # Passive Elastances [mmHg·mL⁻¹]
    0.09, 0.07, 0.05,
    # Active Elastances [mmHg·mL⁻¹]
    0.07, 0.06, 0.55,
    # Reference volumes [mL]
    4.0, 4.0, 10.0,
    # Timing [s]
    0.6, 0.104, 0.68,    # LA
    0.56, 0.064, 0.64,   # RA
    0.0, 0.272, 0.12,    # RV
    # Valves
    0.0075, 75000.0,
    # Heartbeat period and external pressure
    0.8, 0.0
)

# ---------------------------------------------------------------
# Initial conditions and simulation
# ---------------------------------------------------------------
# VLA, VLV, VRA, VRV, pSYSAR, pSYSVEN, pPULAR, pPULVEN, QSYSAR, QSYSVEN, QPULAR, QPULVEN
u0 = [20.0, 60.0, 20.0, 40.0, 90.0, 5.0, 15.0, 8.0, 0.0, 0.0, 0.0, 0.0]
tspan = (0.0, 5.0)

prob = ODEProblem(cardiovascular_odes!, u0, tspan, p)
sol = solve(prob, Tsit5(); saveat=0.001)

p1=plot(sol)
# ---------------------------------------------------------------
# Plot example elastance & ventricular loop
# ---------------------------------------------------------------
tvals = sol.t
VLV = sol[2, :]
VRV = sol[4, :]

# recompute elastances for LV & RV
E_LV = [elastance(t, 0.05, 2.0, 0.0, 0.12, 0.272, 0.12, 0.8) for t in tvals]
E_RV = [elastance(t, 0.05, 0.55, 0.0, 0.12, 0.272, 0.12, 0.8) for t in tvals]

pLV = E_LV .* (VLV .- 10.0)
pRV = E_RV .* (VRV .- 10.0)

plt1 = plot(tvals, E_LV, label="E_LV", lw=2, xlabel="Time [s]", ylabel="Elastance [mmHg/mL]")
plot!(plt1, E_RV, label="E_RV", lw=2, title="Ventricular Elastance")

plt2 = plot(VLV, pLV, label="LV loop", xlabel="Volume [mL]", ylabel="Pressure [mmHg]", lw=2)
plot!(plt2, VRV, pRV, label="RV loop", lw=2, title="Pressure–Volume Loops")

plot(plt1, plt2, layout=(1,2), size=(900,400))