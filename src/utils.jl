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