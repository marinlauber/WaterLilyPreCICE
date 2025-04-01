using WaterLilyPreCICE,WriteVTK,StaticArrays,Interpolations
using Plots,JLD2,OrdinaryDiffEq

# set model parameters
# Cardiac parameters
V0 = 20                     #ml; intercept with volume axis of the ESPVR
Pfilling = 5                #mmHg; venous filling pressure
Pout = 5                    #mmHg: 
EDV = 120                   #ml; end-diastolic volume. We will use EDV with Pvenous to calculate Emin

# Valve resistances
Rmv_fwd = 0.002             #mmHg/ml/s; resistance in forward flow direction
Rmv_bwd = 1e10              #mmHg/ml/s; leak resistance
Rao_fwd = 0.002             #mmHg/ml/s; resistance in forward flow direction
Rao_bwd = 1e10              #mmHg/ml/s; leak resistance

# Arterial model parameters
R_WK2 = 1                   #mmHg/ml/s
C_WK2 = 2                   #ml/mmHg

function Elastance(t;Emin,Emax,a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
    (Emax-Emin) * α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)  + Emin   
end

function Windkessel_P!(du,u,p,t)
    # unpack
    (Q, Pao) = u
    (Pfill,Rmv_fwd,Rmv_bwd,Rao_fwd,Rao_bwd,R,C,Pout)  = p

    # dVdt is a prescribed function
    dVdt = Q

    # pressure in the ventricle
    PLV = dVdt ≥ 0 ? Pfill - dVdt*Rmv_fwd  : -dVdt*Rao_fwd + Pao
    u[1] = PLV # we store the pressure in the ventricle so that the interface can access it

    # what way is the flow going in the aortic valve?
    Qao = dVdt ≤ 0 ? -dVdt : (Pao-PLV)/Rao_bwd # diastole, very small flow

    # rates
    du[1] = 0                        # no rate of change
    du[2] = Qao/C - (Pao-Pout)/(R*C) # dPao/dt
end

tspan = (0,90)
u = [8, 50] # initial conditions
p = [Pfilling, Rmv_fwd, Rmv_bwd, Rao_fwd, Rao_bwd, R_WK2, C_WK2, Pout]

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
A4(ts) = ts<=10 ? 0 : Elastance((ts-10)/20;Emin,Emax)
A5(ts) = ts<=10 ? ts/10 : A1(ts) - Elastance((ts-10)/20;Emin,Emax)

function static_inflation(i,t,interface)
    i==1 && return  0.38*A1(t)
    i==6 && return -0.38*A1(t)
    i in [2,3] && return 0.38*A5(t) # 0.38*A3(t)
    i in [4,5] && return 0.38*A4(t) # 0.38*A2(t)
    i in [7,8] && return  -0.38*A5(t) # -0.38*A3(t)
    i in [9,10] && return -0.38*A4(t) # -0.38*A2(t) 
end

function dynamic_inflation(i,t,interface)
    i==1 && return  0.38*interface.P[end]
    i==6 && return -0.38*interface.P[end]
    i in [2,3] && return 0.38*(interface.P[end] -A4(t))
    i in [4,5] && return 0.38*A4(t)
    i in [7,8] && return  -0.38*(interface.P[end] - A4(t))
    i in [9,10] && return -0.38*A4(t) 
end

# vtk attributes
vtk_srf(a::LumpedInterface) = Float32[el[1] for el in a.srf_id]
vtk_center(a::LumpedInterface) = WaterLilyPreCICE.center.(a.mesh)
vtk_normal(a::LumpedInterface) = WaterLilyPreCICE.normal.(a.mesh)
vtk_dS(a::LumpedInterface) = WaterLilyPreCICE.dS.(a.mesh)
vtk_u(a::LumpedInterface) = a.deformation
vtk_f(a::LumpedInterface) = a.forces

# vtk attributes
custom = Dict("SRF" =>vtk_srf, "center"=>vtk_center, "normal"=>vtk_normal,
              "dS" => vtk_dS, "u" => vtk_u, "f"=>vtk_f)

# coupling interface
interface = LumpedInterface(;surface_mesh="../Solid/geom.inp",
                            func=static_inflation, integrator=integrator)

# make the writer
wr = vtkWriter("pouch"; attrib=custom)
v = []; p = []; t = []; Qs = [] # storage for the volume

while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    readData!(interface)

    # compute the pressure forces
    WaterLilyPreCICE.update!(interface)

    # some post-processing stuff
    push!(v, WaterLilyPreCICE.volume(interface))
    push!(p, A4(sum(interface.dt)))
    push!(t, sum(interface.dt))
    push!(Qs, WaterLilyPreCICE.get_Q(interface))

    # write data to the other participant
    writeData!(interface)

    # do the coupling
    (length(interface.dt)+1)%5==0 && write!(wr,interface)
end
close(wr)
jldsave("pouch_volume.jld2";volume=v,pressure=p,time=t,q=Qs,sol=interface.sol,
                            intP=interface.P,intdt=interface.dt,intV=interface.V)
PreCICE.finalize()


# jldopen("pouch_volume.jld2") do file
#     v = file["volume"]
#     p = file["pressure"]
#     t = file["time"]
#     plot(t,v,label="Volume",lw=2)
#     plot!(t,p,label="Pressure",lw=2,xlims=(0,90))
# end