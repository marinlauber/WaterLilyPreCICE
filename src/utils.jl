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

import WaterLily: @loop
function measure_CPU!(μ₀,μ₁,V,σ::AbstractArray{T,N},body::AbstractBody;t=zero(T),ϵ=1) where {T,N}
    V .= zero(T); μ₀ .= one(T); μ₁ .= zero(T); d²=(2+ϵ)^2
    @fastmath @inline function fill!(μ₀,μ₁,V,d,I)
        d[I] = sdf(body,loc(0,I,T),t,fastd²=d²)
        if d[I]^2<d²
            for i ∈ 1:N
                dᵢ,nᵢ,Vᵢ = measure(body,loc(i,I,T),t,fastd²=d²)
                V[I,i] = Vᵢ[i]
                μ₀[I,i] = WaterLily.μ₀(dᵢ,ϵ)
                for j ∈ 1:N
                    μ₁[I,i,j] = WaterLily.μ₁(dᵢ,ϵ)*nᵢ[j]
                end
            end
        elseif d[I]<zero(T)
            for i ∈ 1:N
                μ₀[I,i] = zero(T)
            end
        end
    end
    @loop fill!(μ₀,μ₁,V,σ,I) over I ∈ inside(σ)
    BC!(μ₀,zeros(SVector{N,T}),false)
    BC!(V ,zeros(SVector{N,T}),false)
end