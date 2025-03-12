using StaticArrays,Plots,OrdinaryDiffEq,PropDicts,Roots

# function residual(v, Vlv, Vrv, p, t)
#     t₁ = e(t)*p.Eesspt*(v-p.Vdspt)     + (1-e(t))*p.P0spt*(exp(p.λspt*(v-p.V0spt))-1.0)
#     t₂ = e(t)*p.Eeslvf*(Vlv-v-p.Vdlvf) + (1-e(t))*p.P0lvf*(exp(p.λlvf*(Vlv-v-p.V0lvf))-1.0)
#     t₃ = e(t)*p.Eesrvf*(Vrv-v-p.Vdrvf) + (1-e(t))*p.P0rvf*(exp(p.λrvf*(Vrv-v-p.V0rvf))-1.0)
#     return -(t₁ - t₂ + t₃)
# end

# from https://doi.org/10.1016/j.medengphy.2003.10.
@inline e(t;A=1,B=80,C=0.375,N=1) = A*exp(-(B*t-C)^N)

# #  Minimal haemodynamic system model including ventricular interaction and valve dynamics https://doi.org/10.1016/j.medengphy.2003.10.
# function System(du,u,t,p)
#     # unpack
#     (Vpa,Vpu,Vlv,Vao,Vvc,Vrv) = u

#     # pericardium P-V relationship
#     Vpcd = Vlv + Vrv
#     Ppcd = p.P0pcd * (exp(p.λpcd*(Vpcd-p.V0pcd)) - 1.0)
#     Pperi = Ppcd + p.Pth

#     # find the septum volume via root finding
#     p.Vspt = find_zero(v->residual(v, Vlv, Vrv, p, t), p.Vspt)

#     # ventricular  pressure-volume relationship
#     Vlvf = Vlv + p.Vspt
#     Vrvf = Vrv - p.Vspt
#     Plvf = e(t)*p.Eeslvf*(Vlvf-p.Vdlvf)   + (1-e(t))*p.P0lvf*(exp(p.λlvf*(Vlvf-p.V0lvf))-1.0)
#     Prvf = e(t)*p.Eesrvf*(Vrvf-p.Vdrvf)   + (1-e(t))*p.P0rvf*(exp(p.λrvf*(Vrvf-p.V0rvf))-1.0)
#     Pspt = e(t)*p.Eesspt*(p.Vspt-p.Vdspt) + (1-e(t))*p.P0spt*(exp(p.λspt*(p.Vspt-p.V0spt))-1.0)
#     Pspt = Plvf - Prvf
#     Plv = Plvf + Pperi
#     Prv = Prvf + Pperi

#     # peripheral chamber pressure volume relationship
#     Ppa = p.Eespa * (Vpa - p.Vdpa)
#     Pao = p.Eesao * (Vao - p.Vdao)
#     Ppu = p.Eespu * (Vpu - p.Vdpu)
#     Pvc = p.Eesvc * (Vvc - p.Vdvc)

#     # populate the dX/dt vector
#     du[1] = (Prv-Ppa)/p.Rpv  - (Ppa-Ppu)/p.Rpul  # dVpa/dt
#     du[2] = (Ppa-Ppu)/p.Rpul - (Ppu-Plv)/p.Rmt   # dVpu/dt
#     du[3] = (Ppu-Plv)/p.Rmt  - (Plv-Pao)/p.Rav   # dVlv/dt
#     du[4] = (Plv-Pao)/p.Rav  - (Pao-Pvc)/p.Rsys  # dVao/dt
#     du[5] = (Pao-Pvc)/p.Rsys - (Pvc-Prv)/p.Rtc   # dVvc/dt
#     du[6] = (Pvc-Prv)/p.Rtc  - (Prv-Ppa)/p.Rpv   # dVrv/dt
#     return nothing
# end

function residual(v, Vlv, Vrv, p, t)
    (Rmt, Rav, Rtc, Rpv, Rpul, Rsys, Ltc, Lpv, 
     Lmt, Lav, Vtot, Pth , period, P0pcd, V0pcd, 
     λpcd, Eeslvf, λlvf, P0lvf, Vdlvf, V0lvf, 
     Eesrvf, λrvf, P0rvf, Vdrvf, V0rvf, Eesspt, 
     Vdspt, P0spt, λspt, V0spt, Eespa, Vdpa, 
     Eespu, Vdpu, Eesao, Vdao, Eesvc, Vdvc, Vspt) = p
    t₁ = e(t)*Eesspt*(v-Vdspt)     + (1-e(t))*P0spt*(exp(λspt*(v-V0spt))-1.0)
    t₂ = e(t)*Eeslvf*(Vlv-v-Vdlvf) + (1-e(t))*P0lvf*(exp(λlvf*(Vlv-v-V0lvf))-1.0)
    t₃ = e(t)*Eesrvf*(Vrv-v-Vdrvf) + (1-e(t))*P0rvf*(exp(λrvf*(Vrv-v-V0rvf))-1.0)
    return -(t₁ - t₂ + t₃)
end

#  Minimal haemodynamic system model including ventricular interaction and valve dynamics https://doi.org/10.1016/j.medengphy.2003.10.
function System(du,u,p,t)
    # unpack
    (Vpa,Vpu,Vlv,Vao,Vvc,Vrv,Qmt,Qav,Qtc,Qpv) = u
    (Rmt, Rav, Rtc, Rpv, Rpul, Rsys, Ltc, Lpv, 
     Lmt, Lav, Vtot, Pth, period, P0pcd, V0pcd, 
     λpcd, Eeslvf, λlvf, P0lvf, Vdlvf, V0lvf, 
     Eesrvf, λrvf, P0rvf, Vdrvf, V0rvf, Eesspt, 
     Vdspt, P0spt, λspt, V0spt, Eespa, Vdpa, 
     Eespu, Vdpu, Eesao, Vdao, Eesvc, Vdvc, Vspt) = p

    # pericardium P-V relationship
    Vpcd = Vlv + Vrv
    Ppcd = P0pcd * (exp(λpcd*(Vpcd-V0pcd)) - 1.0)
    Pperi = Ppcd + Pth

    # find the septum volume via root finding
    @show Vlv, Vrv, Vspt, t
    Vspt = find_zero(v->residual(v, Vlv, Vrv, p, t), Vspt)
    p[end] = Vspt

    # ventricular  pressure-volume relationship
    Vlvf = Vlv + Vspt
    Vrvf = Vrv - Vspt
    Plvf = e(t)*Eeslvf*(Vlvf-Vdlvf) + (1-e(t))*P0lvf*(exp(λlvf*(Vlvf-V0lvf))-1.0)
    Prvf = e(t)*Eesrvf*(Vrvf-Vdrvf) + (1-e(t))*P0rvf*(exp(λrvf*(Vrvf-V0rvf))-1.0)
    Pspt = e(t)*Eesspt*(Vspt-Vdspt) + (1-e(t))*P0spt*(exp(λspt*(Vspt-V0spt))-1.0)
    Pspt = Plvf - Prvf
    Plv = Plvf + Pperi
    Prv = Prvf + Pperi

    # peripheral chamber pressure volume relationship
    Ppa = Eespa * (Vpa - Vdpa)
    Pao = Eesao * (Vao - Vdao)
    Ppu = Eespu * (Vpu - Vdpu)
    Pvc = Eesvc * (Vvc - Vdvc)

    # flow rates with inertial effects
    dQmtdt = (Ppu - Plv - Qmt*Rmt) / Lmt
    dQavdt = (Plv - Pao - Qav*Rav) / Lav
    dQtcdt = (Pvc - Prv - Qtc*Rtc) / Ltc
    dQpvdt = (Prv - Pao - Qpv*Rpv) / Lpv

    # populate the dX/dt vector, dV/dt = Qin - Qout
    du[1] = max(0,(Prv-Ppa)/Rpv)  - (Ppa-Ppu)/Rpul  # dVpa/dt
    du[2] = max(0,(Ppa-Ppu)/Rpul) - (Ppu-Plv)/Rmt   # dVpu/dt
    du[3] = max(0,(Ppu-Plv)/Rmt)  - (Plv-Pao)/Rav   # dVlv/dt
    du[4] = max(0,(Plv-Pao)/Rav)  - (Pao-Pvc)/Rsys  # dVao/dt
    du[5] = max(0,(Pao-Pvc)/Rsys) - (Pvc-Prv)/Rtc   # dVvc/dt
    du[6] = max(0,(Pvc-Prv)/Rtc)  - (Prv-Ppa)/Rpv   # dVrv/dt
    # flow rates
    # Qmt = (Prv-Ppa)/Rpv
    # Qav = (Ppa-Ppu)/Rpul
    # Qtc = (Ppu-Plv)/Rmt
    # Qpv = (Plv-Pao)/Rav
    # Qpa = (Pao-Pvc)/Rsys
    # Qvc = (Pvc-Prv)/Rtc
    du[7] = 0.0
    du[8] = 0.0
    du[9] = 0.0
    du[10] = 0.0
    return nothing
end


# initial state variables
Vlv = 94.6812       # "V_lv in component left_ventricle (mL)"
Vrv = 90.7302       # "V_rv in component right_ventricle (mL)"
Qmt = 245.5813      # "Q_mt in component flow (mL_per_second)"
Qav = 0.0           # "Q_av in component flow (mL_per_second)"
Qtc = 190.0661      # "Q_tc in component flow (mL_per_second)"
Qpv = 0.0           # "Q_pv in component flow (mL_per_second)"
Vpa = 43.0123       # "V_pa in component pulmonary_artery (mL)"
Vpu = 808.4579      # "V_pu in component pulmonary_vein (mL)"
Vao = 133.3381      # "V_ao in component aorta (mL)"
Vvc = 329.7803      # "V_vc in component vena_cava (mL)"

# initial conditions
u₀ = [Vlv, Vrv, Vpa, Vpu, Vlv, Vao, Vvc, Qmt, Qav, Qtc, Qpv]
tspan = (0.0, 100.0)

# params = PropDict("Rmt" => 0.0158,       # "R_mt in component heart_parameters (kPa_second_per_mL)"
#              "Rav" => 0.0180,       # "R_av in component heart_parameters (kPa_second_per_mL)"
#              "Rtc" => 0.0237,       # "R_tc in component heart_parameters (kPa_second_per_mL)"
#              "Rpv" => 0.0055,       # "R_pv in component heart_parameters (kPa_second_per_mL)"
#              "Rpul" => 0.1552,      # "R_pul in component heart_parameters (kPa_second_per_mL)"
#              "Rsys" => 1.0889,      # "R_sys in component heart_parameters (kPa_second_per_mL)"
#              "Ltc" => 8.0093e-5,    # "L_tc in component heart_parameters (kPa_second2_per_mL)"
#              "Lpv" => 1.4868e-4,    # "L_pv in component heart_parameters (kPa_second2_per_mL)"
#              "Lmt" => 7.6968e-5,    # "L_mt in component heart_parameters (kPa_second2_per_mL)"
#              "Lav" => 1.2189e-4,    # "L_av in component heart_parameters (kPa_second2_per_mL)"
#              "Vtot" => 5.5,         # "V_tot in component heart_parameters (mL)"
#              "Pth" => -4.0,         # "P_th in component heart_parameters (kPa)"
#              "period" => 0.75,      # "period in component driver_function (second)"
#              "P0pcd" => 0.5003,     # "P_0_pcd in component pericardium (kPa)"
#              "V0pcd" => 200,        # "V_0_pcd in component pericardium (mL)"
#              "λpcd" => 0.03,        # "lambda_pcd in component pericardium (per_mL)"
#              "Eeslvf" => 2.8798,    # "E_es_lvf in component lvf_calculator (kPa_per_mL)"
#              "λlvf" => 0.033,       # "lambda_lvf in component lvf_calculator (per_mL)"
#              "P0lvf" => 0.1203,     # "P_0_lvf in component lvf_calculator (kPa)"
#              "Vdlvf" => 0.0,        # "V_d_lvf in component lvf_calculator (mL)"
#              "V0lvf" => 0.0,        # "V_0_lvf in component lvf_calculator (mL)"
#              "Eesrvf" => 0.585,     # "E_es_rvf in component rvf_calculator (kPa_per_mL)"
#              "λrvf" => 0.023,       # "lambda_rvf in component rvf_calculator (per_mL)"
#              "P0rvf" => 0.2157,     # "P_0_rvf in component rvf_calculator (kPa)"
#              "Vdrvf" => 0.0,        # "V_d_rvf in component rvf_calculator (mL)"
#              "V0rvf" => 0.0,        # "V_0_rvf in component rvf_calculator (mL)"
#              "Eesspt" => 48.754,    # "E_es_spt in component septum (kPa_per_mL)"
#              "Vdspt" => 2.0,        # "V_d_spt in component septum (mL)"
#              "P0spt" => 1.1101,     # "P_0_spt in component septum (kPa)"
#              "λspt" => 0.435,       # "lambda_spt in component septum (per_mL)"
#              "V0spt" => 2.0,        # "V_0_spt in component septum (mL)"
#              "Eespa" => 0.369,      # "E_es_pa in component pulmonary_artery (kPa_per_mL)"
#              "Vdpa" => 0.0,         # "V_d_pa in component pulmonary_artery (mL)"
#              "Eespu" => 0.0073,     # "E_es_pu in component pulmonary_vein (kPa_per_mL)"
#              "Vdpu" => 0.0,         # "V_d_pu in component pulmonary_vein (mL)"
#              "Eesao" => 0.6913,     # "E_es_ao in component aorta (kPa_per_mL)"
#              "Vdao" => 0.0,         # "V_d_ao in component aorta (mL)"
#              "Eesvc" => 0.0059,     # "E_es_vc in component vena_cava (kPa_per_mL)"
#              "Vdvc" => 0.0,         # "V_d_vc in component vena_cava (mL)"
#              "Vspt" => 0.10)         # "V_spt in component septum (mL)"


# parameters
Rmt = 0.0158
Rav = 0.0180
Rtc = 0.0237
Rpv = 0.0055
Rpul = 0.1552
Rsys = 1.0889
Ltc = 8.0093e-5
Lpv = 1.4868e-4
Lmt = 7.6968e-5
Lav = 1.2189e-4
Vtot = 5.5
Pth = -4.0
period = 0.75
P0pcd = 0.5003
V0pcd = 200
λpcd = 0.03
Eeslvf = 2.8798
λlvf = 0.033
P0lvf = 0.1203
Vdlvf = 0.0
V0lvf = 0.0
Eesrvf = 0.585
λrvf = 0.023
P0rvf = 0.2157
Vdrvf = 0.0
V0rvf = 0.0
Eesspt = 48.754
Vdspt = 2.0
P0spt = 1.1101
λspt = 0.435
V0spt = 2.0
Eespa = 0.369
Vdpa = 0.0
Eespu = 0.0073
Vdpu = 0.0
Eesao = 0.6913
Vdao = 0.0
Eesvc = 0.0059
Vdvc = 0.0
Vspt = 0.1

params = [Rmt, Rav, Rtc, Rpv, Rpul, Rsys, Ltc, Lpv, 
          Lmt, Lav, Vtot, Pth , period, P0pcd, V0pcd, 
          λpcd, Eeslvf, λlvf, P0lvf, Vdlvf, V0lvf, 
          Eesrvf, λrvf, P0rvf, Vdrvf, V0rvf, Eesspt, 
          Vdspt, P0spt, λspt, V0spt, Eespa, Vdpa, 
          Eespu, Vdpu, Eesao, Vdao, Eesvc, Vdvc, Vspt]


# Pass to solver
prob = ODEProblem(System, u₀, tspan, params)
# https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts
res = solve(prob, Tsit5(), dtmax=0.02)

plot(res.t,computePLV.(res.t, getindex.(sol.u, 1)),label="P_\\ LV",lw=2)
plot!(res, idxs = [2] ,linewidth = 2, title = "Windkessel model",
      xaxis = "Time (t/T)", yaxis = "Pressure (mmHg)", label = "P_\\ AO")