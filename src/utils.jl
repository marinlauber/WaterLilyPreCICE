# structure to store fluid state
mutable struct Store
    uˢ:: AbstractArray
    pˢ:: AbstractArray
    b :: AbstractBody
    function Store(sim::Simulation)
        new(copy(sim.flow.u),copy(sim.flow.p),deepcopy(sim.body))
    end
end
# Base.copy(::WaterLily.NoBody) = WaterLily.NoBody()
function store!(s::Store,sim::Simulation)
    s.uˢ .= sim.flow.u; s.pˢ .= sim.flow.p
    s.b = deepcopy(sim.body)
end
function revert!(s::Store,sim::Simulation)
    sim.flow.u .= s.uˢ; sim.flow.p .= s.pˢ; pop!(sim.flow.Δt)
    pop!(sim.pois.n); pop!(sim.pois.n) # pop predictor and corrector
    sim.body = s.b # nice and simple
end

# unpack subarray of increasing values [0,0.5,1,0,0.33,0.66,1] -> [[0,0.5,1],[0,0.33,0.66,1]]
function unpack(a)
    tmp=[a[1]]; ks=Vector{Number}[]
    for i ∈ 2:length(a)
        if a[i]>=a[i-1]
            push!(tmp,a[i])
        else
            push!(ks,tmp); tmp=[a[i]]
        end
    end
    push!(ks,tmp)
    return ks
end
function knotVectorUnpack(knots)
    knots = reshape(knots,reverse(size(knots)))[1,:]
    unpack(knots)
end

function getControlPoints(points, knots)
    points = reshape(points,reverse(size(points)))
    ncp = [length(knot)-count(i->i==0,knot) for knot in knots]
    @assert sum(ncp) == size(points,2) "Number of control points does not match the number of points"
    return [points[:,1+sum(ncp[1:i])-ncp[1]:sum(ncp[1:i])] for i in 1:length(ncp)]
end
function quadPointUnpack(quadPoints)
    quadPoints = reshape(quadPoints,reverse(size(quadPoints)))
    quadPoints = [filter(!isone,filter(!iszero,quadPoints[:,i]))[1] for i in 1:size(quadPoints,2)]
    unpack(quadPoints)
end
function getDeformation(points,knots)
    ncp = [length(knot)-count(i->i==0,knot) for knot in knots]
    return [points[:,1+sum(ncp[1:i])-ncp[1]:sum(ncp[1:i])] for i in 1:length(ncp)]
end