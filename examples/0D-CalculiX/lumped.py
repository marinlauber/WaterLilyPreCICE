# Size of variable arrays:
sizeAlgebraic = 22
sizeStates = 10
sizeConstants = 43
from math import *
from numpy import *

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
 
    constants[0] = 0.0158                    	# "R_mt in component heart_parameters (kPa_second_per_mL)"
    constants[1] = 0.0180                    	# "R_av in component heart_parameters (kPa_second_per_mL)"
    constants[2] = 0.0237                    	# "R_tc in component heart_parameters (kPa_second_per_mL)"
    constants[3] = 0.0055                    	# "R_pv in component heart_parameters (kPa_second_per_mL)"
    constants[4] = 0.1552                    	# "R_pul in component heart_parameters (kPa_second_per_mL)"
    constants[5] = 1.0889                   	# "R_sys in component heart_parameters (kPa_second_per_mL)"
    constants[6] = 8.0093e-5                  # "L_tc in component heart_parameters (kPa_second2_per_mL)"
    constants[7] = 1.4868e-4                  # "L_pv in component heart_parameters (kPa_second2_per_mL)"
    constants[8] = 7.6968e-5                  # "L_mt in component heart_parameters (kPa_second2_per_mL)"
    constants[9] = 1.2189e-4                  # "L_av in component heart_parameters (kPa_second2_per_mL)"
    constants[10] = 5.5                    		# "V_tot in component heart_parameters (mL)"
    constants[11] = -4                    		# "P_th in component heart_parameters (kPa)"
    constants[12] = 1                    			# "A in component driver_function (dimensionless)"
    constants[13] = 80                    		# "B in component driver_function (per_second2)"
    constants[14] = 0.375                    	# "C in component driver_function (second)"
    constants[15] = 0.75                    	# "period in component driver_function (second)"
    constants[16] = 0.5003                    # "P_0_pcd in component pericardium (kPa)"
    constants[17] = 200                    		# "V_0_pcd in component pericardium (mL)"
    constants[18] = 0.03                    	# "lambda_pcd in component pericardium (per_mL)"
    constants[19] = 2.8798                    # "E_es_lvf in component lvf_calculator (kPa_per_mL)"
    constants[20] = 0.033                    	# "lambda_lvf in component lvf_calculator (per_mL)"
    constants[21] = 0.1203                    # "P_0_lvf in component lvf_calculator (kPa)"
    constants[22] = 0                    			# "V_d_lvf in component lvf_calculator (mL)"
    constants[23] = 0                    			# "V_0_lvf in component lvf_calculator (mL)"
    constants[24] = 0.585                    	# "E_es_rvf in component rvf_calculator (kPa_per_mL)"
    constants[25] = 0.023                    	# "lambda_rvf in component rvf_calculator (per_mL)"
    constants[26] = 0.2157                    # "P_0_rvf in component rvf_calculator (kPa)"
    constants[27] = 0                    			# "V_d_rvf in component rvf_calculator (mL)"
    constants[28] = 0                    			# "V_0_rvf in component rvf_calculator (mL)"
    constants[29] = 48.754                    # "E_es_spt in component septum (kPa_per_mL)"
    constants[30] = 2                    			# "V_d_spt in component septum (mL)"
    constants[31] = 1.1101                    # "P_0_spt in component septum (kPa)"
    constants[32] = 0.435                    	# "lambda_spt in component septum (per_mL)"
    constants[33] = 2                    			# "V_0_spt in component septum (mL)"
    constants[34] = 1                    			# "one in component septum (dimensionless)"
    constants[35] = 0.369                    	# "E_es_pa in component pulmonary_artery (kPa_per_mL)"
    constants[36] = 0                    			# "V_d_pa in component pulmonary_artery (mL)"
    constants[37] = 0.0073                    # "E_es_pu in component pulmonary_vein (kPa_per_mL)"
    constants[38] = 0                    			# "V_d_pu in component pulmonary_vein (mL)"
    constants[39] = 0.6913                    # "E_es_ao in component aorta (kPa_per_mL)"
    constants[40] = 0                    			# "V_d_ao in component aorta (mL)"
    constants[41] = 0.0059                    # "E_es_vc in component vena_cava (kPa_per_mL)"
    constants[42] = 0                    			# "V_d_vc in component vena_cava (mL)"

    states[0] = 94.6812                    		# "V_lv in component left_ventricle (mL)"
    states[1] = 90.7302                    		# "V_rv in component right_ventricle (mL)"
    states[2] = 245.5813                    	# "Q_mt in component flow (mL_per_second)"
    states[3] = 0                    					# "Q_av in component flow (mL_per_second)"
    states[4] = 190.0661                    	# "Q_tc in component flow (mL_per_second)"
    states[5] = 0                    					# "Q_pv in component flow (mL_per_second)"
    states[6] = 43.0123                    		# "V_pa in component pulmonary_artery (mL)"
    states[7] = 808.4579                    	# "V_pu in component pulmonary_vein (mL)"
    states[8] = 133.3381                    	# "V_ao in component aorta (mL)"
    states[9] = 329.7803                    	# "V_vc in component vena_cava (mL)"

    return (states, constants)

def computeRates(toi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    rates[0] = custom_piecewise([less(states[2] , 0.00000) & less(states[3] , 0.00000), 0.00000 , less(states[2] , 0.00000), -states[3] , less(states[3] , 0.00000), states[2] , True, states[2]-states[3]])
    rates[1] = custom_piecewise([less(states[4] , 0.00000) & less(states[5] , 0.00000), 0.00000 , less(states[4] , 0.00000), -states[5] , less(states[5] , 0.00000), states[4] , True, states[4]-states[5]])
    
    #C1-3
    # Vpcd = V_lv + V_rv == algebraic[2]
    algebraic[2] = states[0]+states[1]
    # Ppcd == algebraic[3]
    algebraic[3] = constants[16]*(exp(constants[18]*(algebraic[2]-constants[17]))-1.00000)
    #  Pperi == algebraic[4]
    algebraic[4] = algebraic[3]+constants[11]
    
    #C4
    # t*HR/60 == algebraic[0]
    algebraic[0] = custom_piecewise([less_equal(toi , constants[15]), toi , less_equal(toi , constants[15]*2.00000), toi-constants[15] , less_equal(toi , constants[15]*3.00000), toi-constants[15]*2.00000 , less_equal(toi , constants[15]*4.00000), toi-constants[15]*3.00000 , less_equal(toi , constants[15]*5.00000), toi-constants[15]*4.00000 , less_equal(toi , constants[15]*6.00000), toi-constants[15]*5.00000 , less_equal(toi , constants[15]*7.00000), toi-constants[15]*6.00000 , less_equal(toi , constants[15]*8.00000), toi-constants[15]*7.00000 , less_equal(toi , constants[15]*9.00000), toi-constants[15]*8.00000 , less_equal(toi , constants[15]*10.0000), toi-constants[15]*9.00000 , less_equal(toi , constants[15]*11.0000), toi-constants[15]*10.0000 , less_equal(toi , constants[15]*12.0000), toi-constants[15]*11.0000 , less_equal(toi , constants[15]*13.0000), toi-constants[15]*12.0000 , True, float('nan')])
    # e(t) == algebraic[1]
    algebraic[1] = constants[12]*exp(-constants[13]*(power(algebraic[0]-constants[14], 2.00000)))
    # what is this doing
    algebraic = rootfind_0(toi, constants, states, algebraic)
    
    # C5
    # Vlvf - Vdlvf
    algebraic[6] = states[0]-algebraic[5]
    # Eeslvf * (Vlvf - Vdlvf)
    algebraic[7] = constants[19]*(algebraic[6]-constants[22])
    # P0lvf * (exp(lambda_lvf * (Vlvf - V0lvf)) - 1)
    algebraic[8] = constants[21]*(exp(constants[20]*(algebraic[6]-constants[23]))-1.00000)
    # Plvf = e(t) * Eeslvf * (Vlvf - Vdlvf) + (1-e(t)) *P0lvf * (exp(lambda_lvf * (Vlvf - V0lvf)) - 1)
    algebraic[9] = algebraic[1]*algebraic[7]+(1.00000-algebraic[1])*algebraic[8]

    # C11
    algebraic[10] = algebraic[9]+algebraic[4]

    # C15
    algebraic[17] = constants[39]*(states[8]-constants[40])

    rates[3] = custom_piecewise([less(algebraic[10]-algebraic[17] , 0.00000) & less(states[3] , 0.00000), 0.00000 , True, ((algebraic[10]-algebraic[17])-states[3]*constants[1])/constants[9]])
    
    algebraic[11] = states[1]+algebraic[5]
    algebraic[12] = constants[24]*(algebraic[11]-constants[27])
    algebraic[13] = constants[26]*(exp(constants[25]*(algebraic[11]-constants[28]))-1.00000)
    algebraic[14] = algebraic[1]*algebraic[12]+(1.00000-algebraic[1])*algebraic[13]
    algebraic[15] = algebraic[14]+algebraic[4]
    
    # C13
    algebraic[16] = constants[35]*(states[6]-constants[36])+constants[11]
    rates[5] = custom_piecewise([less(algebraic[15]-algebraic[16] , 0.00000) & less(states[5] , 0.00000), 0.00000 , True, ((algebraic[15]-algebraic[16])-states[5]*constants[3])/constants[7]])
    
    # C14
    algebraic[18] = constants[37]*(states[7]-constants[38])+constants[11]
    rates[2] = custom_piecewise([less(algebraic[18]-algebraic[10] , 0.00000) & less(states[2] , 0.00000), 0.00000 , True, ((algebraic[18]-algebraic[10])-states[2]*constants[0])/constants[8]])
    
    # 16
    algebraic[19] = constants[41]*(states[9]-constants[42])
    rates[4] = custom_piecewise([less(algebraic[19]-algebraic[15] , 0.00000) & less(states[4] , 0.00000), 0.00000 , True, ((algebraic[19]-algebraic[15])-states[4]*constants[2])/constants[6]])
    
    # (Eespa*(Vpa-Vdpa) - Eespu*(Vpu-Vdpu))/Rpul
    algebraic[20] = (algebraic[16]-algebraic[18])/constants[4]
    # dVpa/dt = Qin - Qout = Qpv - Qul
    rates[6] = custom_piecewise([less(states[5] , 0.00000), -algebraic[20] , True, states[5]-algebraic[20]])
    rates[6] = max(states[5], 0) - algebraic[20]
    
    
    rates[7] = custom_piecewise([less(states[2] , 0.00000), algebraic[20] , True, algebraic[20]-states[2]])
    rates[7] = algebraic[20] - max(states[2], 0)

    algebraic[21] = (algebraic[17]-algebraic[19])/constants[5]
    
    rates[8] = custom_piecewise([less(states[3] , 0.00000), -algebraic[21] , True, states[3]-algebraic[21]])
    rates[8] = max(states[3], 0) - algebraic[21]

    rates[9] = custom_piecewise([less(states[4] , 0.00000), algebraic[21] , True, algebraic[21]-states[4]])
    rates[9] = algebraic[21] - max(states[4], 0)
    return(rates)

initialGuess0 = None
def rootfind_0(toi, constants, states, algebraic):
    """Calculate value of algebraic variable for DAE"""
    from scipy.optimize import fsolve
    global initialGuess0
    if (initialGuess0 == None): initialGuess0 = 0.1
    if not iterable(toi):
        algebraic[5] = fsolve(residualSN_0, initialGuess0, args=(algebraic, toi, constants, states), xtol=1E-6)
        initialGuess0 = algebraic[5]
    else:
        for (i,t) in enumerate(toi):
            algebraic[5][i] = fsolve(residualSN_0, initialGuess0, args=(algebraic[:,i], toi[i], constants, states[:,i]), xtol=1E-6)
            initialGuess0 = algebraic[5][i]
    return algebraic

def residualSN_0(algebraicCandidate, algebraic, toi, constants, states):
    algebraic[5] = algebraicCandidate
    return (0.00000) - ((((algebraic[1]*constants[29]*(algebraic[5]-constants[30])+(constants[34]-algebraic[1])*constants[31]*(exp(constants[32]*(algebraic[5]-constants[33]))-constants[34]))-algebraic[1]*constants[19]*(states[0]-algebraic[5]))-(1.00000-algebraic[1])*constants[21]*(exp(constants[20]*(states[0]-algebraic[5]))-1.00000))+algebraic[1]*constants[24]*(states[1]+algebraic[5])+(1.00000-algebraic[1])*constants[26]*(exp(constants[25]*(states[1]+algebraic[5]))-1.00000))

def custom_piecewise(cases):
    """Compute result of a piecewise function"""
    return select(cases[0::2],cases[1::2])

def solve_model(LPN_params=False):
    """Solve model with ODE solver"""
    from scipy.integrate import ode
    # Initialise constants and state variables
    (init_states, constants) = initConsts()
    # Allow user to update the LPN parameters
    if LPN_params:                       
        constants=LPN_params

    # Set timespan to solve over
    toi = linspace(0, 7.5, 750)

    # Construct ODE object to solve
    r = ode(computeRates)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
    r.set_initial_value(init_states, toi[0])
    r.set_f_params(constants)

    # Solve model
    states = array([[0.0] * len(toi)] * sizeStates)
    states[:,0] = init_states
    for (i,t) in enumerate(toi[1:]):
        if r.successful():
            r.integrate(t)
            states[:,i+1] = r.y
        else:
            break

    # Compute algebraic variables
    algebraic = computeAlgebraic(constants, states, toi)
    return (toi, states, algebraic)