using Test
using StaticArrays
using WaterLily

@testset "KDTree.jl" begin
    # test the box
    b = Bbox(SA[0.0,0.0,0.0],SA[0.5,0.5-eps(),0.5-eps()])
    @assert b.C == SA[0.0,0.0,0.0] && b.R == SA[0.5,0.5-eps(),0.5-eps()]
    @assert !b.leaf && b.indices == nothing
    # test the split of the box
    left,right = split(b)
    @assert all(left.C.≈SA[-0.25,0.,0.]) && all(left.R.≈SA[0.25,0.5,0.5])
    @assert all(right.C.≈SA[0.25,0.,0.]) && all(right.R.≈SA[0.25,0.5,0.5])
    # test inside
    @assert !inside([0.5+eps(),0.5,0.5],b) && !inside([-0.5-eps(),0.5,0.5],b)
    # test the filer
    points = rand(3,10)
    box = bounding_box(points)
    left, right = split(box)
    @assert length(p_left) + length(p_right) == size(points,2)
    @assert !(p_left in p_right) # check that there is no point in both boxes

    # test the tree
    points = rand(3,100)
    t = Tree(points)
    @assert length(t) == 1
    @assert t[1].C == SA[0.0,0.0,0.0] && t[1].R == SA[0.5,0.5-eps(),0.5-eps()]
    @assert !t[1].leaf && t[1].indices == nothing
end
@testset "MeshBodies.jl" begin
    # inside bbox or not
    rect = Rect(0,0,0,1,1,1) # origin and widths
    @assert !inside(SA[0.5,1,2.5],rect) && inside(SA[0.5,0.5,0.5],rect)
 
    # distance to box center
    rect = Rect(0,0,0,1,1,1) # origin and widths
    @assert dist(SA[1.0,1.0,1.0],rect) == √0.75
    @assert dist(SA[1.5,1.0,1.0],rect) == √1.5
    @assert dist(SA[1.5,1.5,1.0],rect) == √2.25
    @assert dist(SA[1.5,1.5,1.5],rect) == √3.0
end