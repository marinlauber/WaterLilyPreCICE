"""
TODO:
  - [ ] monitor every step, the pressures and volume, save with time step and interation counter to be able to plot
  - [ ] hold the displacement until convergence, remove force critertion sice it's not used
  - [ ] match the ventricular volume to the actuaiton volume when we know it must balance (roughly)
"""

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
A1 = linear_interpolation([0.,5.,100], [0.,1.,1.])
A2 = linear_interpolation([0.,5.,100.], [0.,0.,10.])
A3 = linear_interpolation([0.,5.,100.], [0.,1.,-9.])
A4(ts) = ts<=5 ? 0 : 2Gaussian((ts-5)/20;μ=0.35,σ=0.08) #Elastance((ts-10)/20)
A5(ts) = ts<=5 ? ts/5 : A1(ts) - 2Gaussian((ts-5)/20;μ=0.35,σ=0.08) #Elastance((ts-10)/20)

kPa2mmHg = 7.50062
function dynamic_inflation(i,t,interface;mmHg2Kpa=0.133322,scale=12.5,Pa=0.35/(mmHg2Kpa*scale))
    i==1 && return  interface.P[end]*mmHg2Kpa*scale
    i==6 && return -interface.P[end]*mmHg2Kpa*scale
    i in [2,3] && return (interface.P[end] - Pa*A4(t))*mmHg2Kpa*scale
    i in [4,5] && return Pa*A4(t)*mmHg2Kpa*scale
    i in [7,8] && return  -(interface.P[end] - Pa*A4(t))*mmHg2Kpa*scale
    i in [9,10] && return -Pa*A4(t)*mmHg2Kpa*scale
end

# vtk attributes
vtk_srf(a::LumpedInterface) = Float32[el[1] for el in a.srf_id]
vtk_center(a::LumpedInterface) = WaterLilyPreCICE.center.(a.mesh)
vtk_normal(a::LumpedInterface) = WaterLilyPreCICE.normal.(a.mesh)
vtk_dS(a::LumpedInterface) = WaterLilyPreCICE.dS.(a.mesh)
vtk_u(a::LumpedInterface) = a.deformation
vtk_f(a::LumpedInterface) = a.forces
vtk_vis(a::LumpedInterface) = Float32[ifelse(el[1] ∈ [1,6], 0, ifelse(el[1] ∈ [2,3,7,8], 1, 2)) for el in a.srf_id]
# vtk_du(a::LumpedInterface) =

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
run(`rm -rf vtk_data/ pouch.pvd`) #clean
wr = vtkWriter("pouch"; attrib=custom)

save!(wr,interface)
v = []; p = []; t = []; Qs = [] # storage for the volume

storage_step = [] # store the iteration, the pressure and the volume at every coupling step

global iteration = 1
global step = 1

# TODO, I know this from the inital conditions, but I should get it from the interface
global target_vol = 0.10136861655795992
@show interface.V
pop!(interface.V) # remove the first volume
push!(interface.V, target_vol) # add the target volume
global old_P = 0.0

# global Ps = []
# global Vs = []

while PreCICE.isCouplingOngoing()

    # read the data from the other participant
    readData!(interface) # sets sum(Δt) = t+Δt, and updates mesh
    vol = WaterLilyPreCICE.volume(interface)
    push!(interface.V, vol)

    # TODO start by a pressure ramp until we fill the ventricle to EDP
    if sum(interface.Δt) < 50
        new_P = 0.21*A1(sum(interface.Δt))
    # we are in the isovolumetric contraction
    else
        println("isovolumetric iteration")
        # no flow inside the ventricle requires a large compliance
        Cp = 1/100.
        dVdt = sum(interface.V[end-2:end].*SA[1.,-4.,3.])/2interface.Δt[end] # this eventually goes to zero
        dPdt = sum(interface.P[end-2:end].*SA[1.,-4.,3.])/2interface.Δt[end]
        dPdV = dPdt*(inv(dVdt)+eps())
        # Cp_test = diff(interface.P[end-1:end]) / diff(interface.V[end-1:end])
        @show dPdV, Cp
        # vol = WaterLilyPreCICE.volume(interface)
        # how much we need to change the pressure
        new_P = old_P + Cp*(target_vol - interface.V[end])
    end
    @show target_vol, vol
    @show target_vol - vol
    @show old_P, new_P
    @show old_P - new_P
    push!(interface.P, new_P)

    # store the data for that step
    push!(storage_step, [step, iteration, new_P, vol])
    global iteration += 1 # increment

    # update the ODEs, compute the forces on the mesh
    # WaterLilyPreCICE.update!(interface; integrate=false)
    WaterLilyPreCICE.get_forces!(interface)
    push!(interface.sol, [1,1,1]) # prevents the solver from crashing

    # write data to the other participant and advance the coupling
    writeData!(interface)
    #TODO should this be in a check for reading checkpoint?
    global old_P = new_P

    # if we have converged, save if required
    if PreCICE.isTimeWindowComplete()
        # reset iteration
        global iteration = 1
        global step += 1 # increment the step
        # save every n step
        (length(interface.Δt)+1)%10==1 && save!(wr,interface)

        global target_vol = volume(interface) # same as `vol``

         # some post-processing stuff
        push!(v, WaterLilyPreCICE.volume(interface))
        push!(p, A4(sum(interface.Δt)))
        push!(t, sum(interface.Δt))
        push!(Qs, WaterLilyPreCICE.get_Q(interface))
        println("")
        # save everytime so we avoid not having data
        jldsave("pouch_volume.jld2";volume=v,pressure=p,time=t,q=Qs,sol=interface.sol,
                                    intP=interface.P,intdt=interface.Δt,intV=interface.V,
                                    step=storage_step)
    end
end
close(wr)
PreCICE.finalize()