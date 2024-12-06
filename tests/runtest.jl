# inside bbox or not
outside(x::SVector,bbox::Rect) = !(all(bbox.origin .≤ x) && all(x .≤ bbox.origin+bbox.widths)) # 1.679 ns (0 allocations: 0 bytes)
rect = Rect(0,0,0,1,1,1) # origin and widths
# @assert !inside(SA[0.5,1,2.5],rect) && inside(SA[0.5,0.5,0.5],rect)
# 
# distance to box center
dist(x::SVector,bbox::Rect) = √sum(abs2,x.-bbox.origin-0.5bbox.widths) # 1.707 ns (0 allocations: 0 bytes)
rect = Rect(0,0,0,1,1,1) # origin and widths
# @assert dist(SA[1.0,1.0,1.0],rect) == √0.75
# @assert dist(SA[1.5,1.0,1.0],rect) == √1.5
# @assert dist(SA[1.5,1.5,1.0],rect) == √2.25
# @assert dist(SA[1.5,1.5,1.5],rect) == √3.0
