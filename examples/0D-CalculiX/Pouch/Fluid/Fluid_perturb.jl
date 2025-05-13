using WaterLilyPreCICE,WriteVTK,StaticArrays,Interpolations
using Plots,JLD2,OrdinaryDiffEq

function Elastance(t;a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)
end
function Gaussian(t;μ=0,σ=0.1)
    exp(-((t%1 - μ) / σ)^2 / 2)
end

function Windkessel_P!(du,u,p,t)
    #unpack
    (Q,Plv,Pa,Pv) = u
    (Ra,Ca,Rv,Cv,Rp) = p

    # dVdt is a prescribed function
    dVdt = Q
    
    # what way is the flow going in the aortic valve?
    Qa = dVdt < -eps() ? -dVdt : (Pa - Plv)/1e5 # diastole, very small flow
    Qv = dVdt >  eps() ?  dVdt : (Plv - Pv)/1e5 # systole, very small flow

    # pressure in the ventricle
    Plv = dVdt > eps() ? Pv - dVdt*Rv : (dVdt < -eps() ? -dVdt*Ra + Pa : Plv)
    abs(dVdt) > 1e-6 && (u[1] = Plv) # store to access after
    
    # rates
    du[1] = 0 # empty
    du[2] = 0                     # dPlv/dt  
    du[3] = Qa/Ca - (Pa-Pv)/(Rp*Ca)  # dPa/dt
    du[4] = (Pa-Pv)/(Rp*Cv) - Qv/Cv  # dVv/dt  
end

tspan = (0,90)
u = [0, 60, 70, 8] # initial conditions for Plv, Pa, Pv
Ra = 8e6 * 1.333e-8
Rp = 3e8 * 1.333e-8
Rv = 1e6 * 1.333e-8
Ca = 8e-9 * 1.333e8
Cv = 5e-8 * 1.333e8
p = (Ra,Ca,Rv,Cv,Rp)

prob = ODEProblem(Windkessel_P!, u, tspan, p)
integrator = init(prob, Tsit5(), dtmax=1e-3, reltol=1e-6, abstol=1e-9,
                  save_everystep=false)

# parameters
Emax = 2      # mmHg/ml; slope of the ESPVR
Emin = 0.05   # mmHg/ml
HR = 60       # heart rate in beats/min

# this are the loading curves inside CalculiX
A1 = linear_interpolation([0.,10.,100], [0.,1.,1.])
A2 = linear_interpolation([0.,10.,100.], [0.,0.,10.])
A3 = linear_interpolation([0.,10.,100.], [0.,1.,-9.])
A4(ts) = ts<=10 ? 0 : 2Gaussian((ts-10)/20;μ=0.35,σ=0.08) #Elastance((ts-10)/20)
A5(ts) = ts<=10 ? ts/10 : A1(ts) - 2Gaussian((ts-10)/20;μ=0.35,σ=0.08) #Elastance((ts-10)/20)

function static_inflation(i,t,interface)
    i==1 && return  0.35*A1(t)
    i==6 && return -0.35*A1(t)
    i in [2,3] && return 0.35*A5(t) # 0.38*A3(t)
    i in [4,5] && return 0.35*A4(t) # 0.38*A2(t)
    i in [7,8] && return  -0.35*A5(t) # -0.38*A3(t)
    i in [9,10] && return -0.35*A4(t) # -0.38*A2(t) 
end
kPa2mmHg = 7.50062
function dynamic_inflation(i,t,interface;mmHg2Kpa=0.133322,scale=12.5,Pa=0.35/(mmHg2Kpa*scale))
    i==1 && return  interface.P[end]*mmHg2Kpa*scale
    i==6 && return -interface.P[end]*mmHg2Kpa*scale
    i in [2,3] && return (interface.P[end] - Pa*A4(t))*mmHg2Kpa*scale
    i in [4,5] && return Pa*A4(t)*mmHg2Kpa*scale
    i in [7,8] && return  -(interface.P[end] - Pa*A4(t))*mmHg2Kpa*scale
    i in [9,10] && return -Pa*A4(t)*mmHg2Kpa*scale
end

no_inflation(args...) = 0.0

# vtk attributes
vtk_srf(a::LumpedInterface) = Float32[el[1] for el in a.srf_id]
vtk_center(a::LumpedInterface) = WaterLilyPreCICE.center.(a.mesh)
vtk_normal(a::LumpedInterface) = WaterLilyPreCICE.normal.(a.mesh)
vtk_dS(a::LumpedInterface) = WaterLilyPreCICE.dS.(a.mesh)
vtk_u(a::LumpedInterface) = a.deformation
vtk_f(a::LumpedInterface) = a.forces
vtk_vis(a::LumpedInterface) = Float32[ifelse(el[1] ∈ [1,6], 0, ifelse(el[1] ∈ [2,3,7,8], 1, 2)) for el in a.srf_id]

# vtk attributes
custom = Dict("SRF" =>vtk_srf, "center"=>vtk_center, "normal"=>vtk_normal,
              "dS" => vtk_dS, "u" => vtk_u, "f"=>vtk_f, "vis"=>vtk_vis)

# coupling interface
interface = LumpedInterface(surface_mesh="../Solid/geom.inp",
                            func=dynamic_inflation,
                            integrator=integrator)

WaterLilyPreCICE.get_forces!(interface) # since we restart in the FEA, we need to compute the forces
interface.forces .= 0.0

# make the writer
wr = vtkWriter("pouch"; attrib=custom)
write!(wr,interface)
v = []; p = []; t = []; Qs = [] # storage for the volume

global iteration = 1

#@TODO, I know this from the inital conditions, but I should get it from the interface
global target_vol = 0.10136861655795992
global old_P = 0.0

while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    readData!(interface) # sets sum(Δt) = t+Δt
    @show sum(interface.deformation,dims=1)
        
    # @TODO pressure perturbation
    # check what iteration we are in, first one we perturb the pressure and the use the normal on
    # if iteration==1
        # println("iteration 1, perturbing the pressure")
        # push!(interface.P, 0.25*A1(sum(interface.Δt))) # 0.35kPa ventricular pressure
    # else
        # println("iteration > 1, using the correct pressure")
        # push!(interface.P, 0.21*A1(sum(interface.Δt))) # 0.35kPa ventricular pressure
    # end

    # @TODO constant flow inside the ventricle
    Cp = 1/100.
    vol = WaterLilyPreCICE.volume(interface)
    @show target_vol, vol
    
    # how much we need to change the pressure
    new_P = old_P + Cp*(target_vol - vol)

    # first ever iteration the volume is unknown yet, we simply prescribe a pressure
    if iteration==1
        new_P = 0.21*A1(sum(interface.Δt))
    end
    push!(interface.P, new_P)
    
    #@TODO classic constant pressure inflation
    # push!(interface.P, 0.21*A1(sum(interface.Δt))) # 0.35kPa ventricular pressure
    @show sum(interface.Δt), interface.P[end]
    @show new_P - old_P
    # @show interface.P[end] - old_P
 
    
    # update the ODEs, compute the forces on the mesh
    WaterLilyPreCICE.update!(interface; integrate=false)
    @show sum(interface.forces,dims=1)
    push!(interface.sol, [1,1,1]) # prevents the solver from crashing

    # write data to the other participant and advance the coupling
    writeData!(interface)
    global iteration += 1
    # global target_vol = vol
    global old_P = new_P

    # if we have converged, save if required
    if PreCICE.isTimeWindowComplete()

        global target_vol = target_vol + interface.Δt[end]*0.05/10 # flow rate desired to reach 0.05/10

        # global iteration = 1
        # (length(interface.Δt)+1)%1==0 && 
        write!(wr,interface)
         # some post-processing stuff
        push!(v, WaterLilyPreCICE.volume(interface))
        push!(p, A4(sum(interface.Δt)))
        push!(t, sum(interface.Δt))
        push!(Qs, WaterLilyPreCICE.get_Q(interface))
        println("")
        # save everytime so we avoid not having data
        jldsave("pouch_volume.jld2";volume=v,pressure=p,time=t,q=Qs,sol=interface.sol,
                                    intP=interface.P,intdt=interface.Δt,intV=interface.V)
    end
end
close(wr)
PreCICE.finalize()