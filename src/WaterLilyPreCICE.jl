module WaterLilyPreCICE

using WaterLily
using StaticArrays
using PreCICE

using Reexport
@reexport using WaterLily
@reexport using PreCICE

abstract type AbstractInterface end
export AbstractInterface

# include helpers
include("utils.jl")

"""
"""
mutable struct CoupledSimulation <: AbstractSimulation
    sim :: Simulation
    int :: AbstractInterface
    store :: Store
end
# overlead properties
Base.getproperty(f::CoupledSimulation, s::Symbol) = s in propertynames(f) ? getfield(f, s) : getfield(f.sim, s)

function CoupledSimulation(args...; T=Float64, mem=Array, interface=:CalculiXInterface, # WL specific
                           fname="geom.inp", map=(x,t)->x, scale=1.f0, # CalculiX specific
                           dir=nothing, curves=nothing, center=SA[0.,0.], # Gismo specific
                           func=(i,t)->0, prob=nothing,                   # Lumped specific
                           kwargs...) # args and kwargs are passed to Simulation

    # check that the interface exists
    @assert interface in [:GismoInterface,:CalculiXInterface,
                          :LumpedInterface] "The interface specified does not exist"
    if interface==:CalculiXInterface
        @assert length(args[1])==3 "Pure 2D simulations are not support for CalculiX coupling, please make it 3D (Nx,Ny)->(Nx,Ny,1)"
    elseif interface==:GismoInterface
        @assert length(args[1])==2 "3D simulations are not support for Gismo coupling"
    end
   
     # keyword aguments might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # create coupling
    PreCICE.createParticipant("WaterLily", configFileName, 0, 1)

    # generate the coupling interface from the coupling partner
    int, body = eval(interface)(T; fname, map ,scale, dir, curves, center, func, prob)

    # the simulation
    sim = Simulation(args...; mem, T, body, kwargs...)

    # storage for iteration
    store =  Store(sim)
    
    # return coupled sim
    CoupledSimulation(sim, int, store)
end


# overload sim_step!
import WaterLily: sim_step!
function sim_step!(sim::CoupledSimulation; kwargs...)
    # update the this participant this is type-specialized
    WaterLilyPreCICE.update!(sim.int, sim; kwargs...)
    # measure and momentum step
    measure!(sim); mom_step!(sim.flow, sim.pois)
    # get forces, this is type-specialized
    get_forces!(sim.int, sim.flow, sim.body; kwargs...)
end

readData!(sim::CoupledSimulation) = readData!(sim.int, sim.sim, sim.store)
function readData!(interface::AbstractInterface,sim::Simulation,store::Store)
    # set time step
    dt_precice = PreCICE.getMaxTimeStepSize()
    interface.dt[end] = dt_precice # min(sim.flow.Δt[end]*sim.U/sim.L, dt_precice) # min physical time step
    sim.flow.Δt[end] = interface.dt[end]*sim.L/sim.U # numerical time step

    if PreCICE.requiresWritingCheckpoint()
        store!(store,sim)
    end

    # Read control point displacements, type specific
    readData!(interface)
end

writeData!(sim::CoupledSimulation) = writeData!(sim.int, sim.sim, sim.store)
function writeData!(interface::AbstractInterface,sim::Simulation,store::Store)
    # write the force at the integration points
    writeData!(interface)

    # do the coupling    
    interface.dt[end] = PreCICE.advance(interface.dt[end])
    
    # read checkpoint if required or move on
    if PreCICE.requiresReadingCheckpoint()
        revert!(store,sim)
    end
end
export CoupledSimulation,readData!,writeData!

# CalculiX specific interface functions
include("CalculiXInterface.jl")
export CalculiXInterface,MeshBody,Tree,Bbox
@reexport using GeometryBasics

# G+Smo specific interface functions
include("GismoInterface.jl")
export GismoInterface

include("LumpedInterface.jl")
export LumpedInterface

end # module