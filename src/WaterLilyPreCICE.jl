module WaterLilyPreCICE

using WaterLily
using StaticArrays

# structure to store fluid state
mutable struct Store
    uˢ:: AbstractArray
    pˢ:: AbstractArray
    b :: AbstractBody
    function Store(sim::Simulation)
        new(copy(sim.flow.u),copy(sim.flow.p),copy(sim.body))
    end
end
function store!(s::Store,sim::Simulation)
    s.uˢ .= sim.flow.u; s.pˢ .= sim.flow.p
    s.b = copy(sim.body)
end
function revert!(s::Store,sim::Simulation)
    sim.flow.u .= s.uˢ; sim.flow.p .= s.pˢ; pop!(sim.flow.Δt)
    pop!(sim.pois.n); pop!(sim.pois.n) # pop predictor and corrector
    sim.body = s.b # nice and simple
end

abstract type AbstractInterface end
export AbstractInterface

using PreCICE
function initialize!(U,L;interface=:CalculiXInterface,kwargs...)
    
    # check that the interface exists
    @assert interface in [:GismoInterface,:CalculiXInterface] "The interface specified does not exist"

    # keyword aguments might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # create coupling
    PreCICE.createParticipant("Fluid", configFileName, 0, 1)

    # generate the coupling interface from the coupling partener
    return eval(interface)(U,L;kwargs...)
end

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

function update!(interface::AbstractInterface,sim::Simulation;kwargs...)
    @warn "not implemented"
end

function step!(flow::Flow,pois::AbstractPoisson,body,interface::AbstractInterface)
    mom_step!(flow,pois)
    getInterfaceForces!(interface,flow,body)
end

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
export initialize!,Store,readData!,update!,step!,writeData!

# CalculiX specific interface functions
include("CalculiXInterface.jl")
export CalculiXInterface

# G+Smo specific interface functions
include("GismoInterface.jl")
export GismoInterface

end # module