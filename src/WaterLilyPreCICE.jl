module WaterLilyPreCICE

using WaterLily
using StaticArrays

using Reexport
@reexport using WaterLily

# include helpers
include("utils.jl")


using PreCICE
@reexport using PreCICE
mutable struct CoupledSimulation <: AbstractSimulation
    sim :: Simulation
    int :: AbstractInterface
    store :: Store
    function CoupledSimulation(args...;  mem=Array, interface=:CalculiXInterface, kwargs...)
        
        # check that the interface exists
        @assert interface in [:GismoInterface,:CalculiXInterface,
                              :LumpedInterface] "The interface specified does not exist"

        # keyword aguments might be specified
        if size(ARGS, 1) < 1
            configFileName = "precice-config.xml"
        else
            configFileName = ARGS[1]
        end

        # create coupling
        PreCICE.createParticipant("Fluid", configFileName, 0, 1)

        # generate the coupling interface from the coupling partner
        interface, body = eval(interface)(U,L;kwargs...)

        # WaterLily Simulation
        sim = Simulation(args...; mem, body, kwargs...)
        # storage for iteration
        store =  Store(sim)
        # return sim type
        new(sim, interface, store)
    end
end
# overlead properties
Base.getproperty(f::CoupledSmiulation, s::Symbol) = s in propertynames(f) ? getfield(f, s) : getfield(f.sim, s)

abstract type AbstractInterface end
export AbstractInterface

# function initialize!(U,L;interface=:CalculiXInterface,kwargs...)
    
#     # check that the interface exists
#     @assert interface in [:GismoInterface,:CalculiXInterface,
#                           :LumpedInterface] "The interface specified does not exist"

#     # keyword aguments might be specified
#     if size(ARGS, 1) < 1
#         configFileName = "precice-config.xml"
#     else
#         configFileName = ARGS[1]
#     end

#     # create coupling
#     PreCICE.createParticipant("Fluid", configFileName, 0, 1)

#     # generate the coupling interface from the coupling partener
#     return eval(interface)(U,L;kwargs...)
# end

# function readData!(interface::AbstractInterface,sim::Simulation,store::Store)
function readData!(sim::CoupledSimulation)
    # set time step
    dt_precice = PreCICE.getMaxTimeStepSize()
    sim.interface.dt[end] = dt_precice # min(sim.flow.Δt[end]*sim.U/sim.L, dt_precice) # min physical time step
    sim.flow.Δt[end] = sim.interface.dt[end]*sim.L/sim.U # numerical time step

    if PreCICE.requiresWritingCheckpoint()
        store!(sim.store,sim)
    end

    # Read control point displacements, type specific
    readData!(sim.interface)
end

# macro notimplemented(message="not implemented")
    # quote
        # error($esc(message))
    # end
# end

# function update!(interface::AbstractInterface,args...;kwargs...)
    # @notimplemented "`update!` not implemented for `AbstractInterface`, it must to be (type) specialized"
# end

import WaterLily: sim_step!
function sim_step!(sim::CoupledSimulation)
    mom_step!(flow,pois)
    getInterfaceForces!(interface,flow,body)
end

# function writeData!(interface::AbstractInterface,sim::Simulation,store::Store)
function writeData!(sim::CoupledSimulations)
    # write the force at the integration points
    writeData!(interface)

    # do the coupling    
    interface.dt[end] = PreCICE.advance(interface.dt[end])
    
    # read checkpoint if required or move on
    if PreCICE.requiresReadingCheckpoint()
        revert!(store,sim)
    end
end
export initialize!,Store,readData!,update!,step!,writeData!

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