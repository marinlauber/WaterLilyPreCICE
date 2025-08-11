using CirculatorySystemModels
using ModelingToolkit
using OrdinaryDiffEq
using Plots

## Cycle time in seconds
τ = 0.85

# Resistances and Compliances
R_s = 1.11
C_sv = 11.0
MCFP = 7.0

# Set up the model elements
@variables t

# 
@named Rs = Resistor(R=R_s)

# the sphere is just an elastance
@named Sphere = Elastance(V₀=1.0, E=1.0, inP=false)

# constant pressure source
@named pressure = ConstantPressure(P=1.0)
# @named flow = ConstantFlow(Q=1.0)

circ_eqs = [
    connect(pressure.out, Rs.in),
    connect(Rs.out, Sphere.in),
    connect(Sphere.out, pressure.in),
]
# Add the component equations
@named _circ_model = ODESystem(circ_eqs, t)

@named circ_model = compose(_circ_model, [Rs, Csv])

# Simplify the ODE system
circ_sys = structural_simplify(circ_model)

# initial conditions
u0 = [Csv.p => MCFP]

# Then we can define the problem:
tspan = (0, 20)
prob = ODEProblem(circ_sys, u0, tspan)
integrator = init(prob, Vern7(), reltol=1e-6, abstol=1e-9, save_everystep=true)

dt = 1
for i in 1:20
    # step exactly to `t+dt`
    step!(integrator, dt, true)
end

set_u!(integrator, u0)
# reinit!(integrator, u0)


# # # # Simulate
# @which step!(integrator, dt, false)
# @time sol = solve(prob, Vern7(), reltol=1e-6, abstol=1e-9, saveat=range(20.9τ,21.9τ,length=1000))

# tspan = (20.9τ,21.9τ)
# p1 = plot(sol, idxs=[LV.p, Csa.in.p], tspan=tspan, xlabel = "Time [s]", ylabel = "Pressure [mmHg]",  hidexaxis = nothing) # Make a line plot
# p2 = plot(sol, idxs=[LV.V], tspan=tspan,xlabel = "Time [s]", ylabel = "Volume [ml]",  linkaxes = :all)
# p3 = plot(sol, idxs=[Csa.in.q,Csv.in.q], tspan=tspan,xlabel = "Time [s]", ylabel = "Flow rate [ml/s]", linkaxes = :all)
# p4 = plot(sol, idxs=(LV.V, LV.p),xlabel = "Volume [ml]", ylabel = "Pressure [mmHg]", linkaxes = :all)

# img = plot(p1, p2, p3, p4; layout=@layout([a b; c d]), legend = true)
# # # savefig(img,"single_chamber_model.png")
