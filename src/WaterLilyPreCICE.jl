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

# CalculiX specific interface functions
include("CalculiXInterface.jl")
export CalculiXInterface,MeshBody,Tree,Bbox
@reexport using GeometryBasics

# G+Smo specific interface functions
include("GismoInterface.jl")
export GismoInterface

include("LumpedInterface.jl")
export LumpedInterface

"""
"""
mutable struct CoupledSimulation <: AbstractSimulation
    sim :: Simulation
    int :: AbstractInterface
    store :: Store
end
# overlead properties
Base.getproperty(f::CoupledSimulation, s::Symbol) = s in propertynames(f) ? getfield(f, s) : getfield(f.sim, s)

function CoupledSimulation(args...; fname="geom.inp", map=(x,t)->x, scale=1.f0, T=Float64,
                            mem=Array, interface=:CalculiXInterface, kwargs...)
    # check that the interface exists
    @assert interface in [:GismoInterface,:CalculiXInterface,
    :LumpedInterface] "The interface specified does not exist"

    # generate the coupling interface from the coupling partner
    # interface, body = CalculiXInterface(T; fname, map ,scale, kwargs...)
    body = MeshBody(fname ; map, scale, T)
    sim = Simulation(args...; mem, T, body, kwargs...)
    # storage for iteration
    store =  Store(sim)
   
     # keyword aguments might be specified
    if size(ARGS, 1) < 1
        configFileName = "precice-config.xml"
    else
        configFileName = ARGS[1]
    end

    # create coupling
    PreCICE.createParticipant("WaterLily", configFileName, 0, 1)

    numberOfVertices, dimensions = length(body.mesh.position), 3
    vertices = Array{Float64,2}(undef, numberOfVertices, dimensions)
    for i in 1:numberOfVertices
    vertices[i,:] .= body.mesh.position[i].data
    end
    ControlPointsID = PreCICE.setMeshVertices("Fluid-Mesh", vertices)
    # ControlPointsID = vertices

    # storage arrays
    forces = zeros(typeof(scale), size(vertices))
    deformation = zeros(typeof(scale), size(vertices))
    map_id = zeros(10)

    # initilise PreCICE
    PreCICE.initialize()
    dt = PreCICE.getMaxTimeStepSize()
    dt = 0.5
    interface = CalculiXInterface(ControlPointsID, vertices, forces, deformation, map_id, [dt])

    # new(sim, store)
    CoupledSimulation(sim, interface, store)
end


# # overload sim_step!
# import WaterLily: sim_step!
# function sim_step!(sim::CoupledSimulation;kwargs...)
#     # update the this participant this is type-specialized
#     WaterLilyPreCICE.update!(sim.interface, sim; kwargs...)
#     # measure and momentum step
#     measure!(sim); mom_step!(sim.flow,sim.pois)
#     # get forces, this is type-specialized
#     get_forces!(sim.interface,sim.flow,sim.body; kwargs...)
# end


function initialize!(U,L;interface=:CalculiXInterface,kwargs...)
    
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
    PreCICE.createParticipant("WaterLily", configFileName, 0, 1)

    # generate the coupling interface from the coupling partener
    return eval(interface)(;kwargs...)
end

macro notimplemented(message="not implemented")
    quote
        error($esc(message))
    end
end

function update!(interface::AbstractInterface, args...; kwargs...)
    @notimplemented "`update!` not implemented for `AbstractInterface`, it must to be (type) specialized"
end

function step!(flow::Flow,pois::AbstractPoisson,body,interface::AbstractInterface)
    mom_step!(flow,pois)
    getInterfaceForces!(interface,flow,body)
end

readData!(sim::CoupledSimulation) = readData!(sim.interface, sim.sim, sim.store)
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

writeData!(sim::CoupledSimulation) = writeData!(sim.interface, sim.sim, sim.store)
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
export CoupledSimulation,initialize!,Store,readData!,update!,step!,writeData!

end # module