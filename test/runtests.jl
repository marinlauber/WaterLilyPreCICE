using Test
using StaticArrays
using WaterLilyPreCICE

@testset "MeshBodies.jl" begin
    # inside bbox or not
    rect = Rect(0,0,0,1,1,1) # origin and widths
    @test WaterLilyPreCICE.outside(SA[0.5,1,2.5],rect) && 
           !WaterLilyPreCICE.outside(SA[0.5,0.5,0.5],rect)
    # distance to box center
    rect = Rect(0,0,0,1,1,1) # origin and widths
    @test WaterLilyPreCICE.outside(SA[-1,0,0],rect)
    @test !WaterLilyPreCICE.outside(SA[ 1,0,0],rect)
    @test WaterLilyPreCICE.dist(SA[1.0,1.0,1.0],rect) == √0.75
    @test WaterLilyPreCICE.dist(SA[1.5,1.0,1.0],rect) == √1.5
    @test WaterLilyPreCICE.dist(SA[1.5,1.5,1.0],rect) == √2.25
    @test WaterLilyPreCICE.dist(SA[1.5,1.5,1.5],rect) == √3.0
    # make a single tri mesh
    mesh = GeometryBasics.Mesh(Point{3, Float64}.([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.]]), 
                               TriangleFace{Int}(1, 2, 3))
    # test tri operation
    @test all(WaterLilyPreCICE.locate(mesh[1], SA[0.,0.,0.]) .≈ [0.,0.,0.]) 
    @test all(WaterLilyPreCICE.locate(mesh[1], SA[1.,1.,1.]) .≈ [1.0,1.0,0.0])
    @test all(WaterLilyPreCICE.locate(mesh[1], SA[0.,1.,0.]) .≈ [.5,.5,0.])
    z = 5-10rand() # anywhere above should project onto it
    @test all(WaterLilyPreCICE.locate(mesh[1], SA[0.33,0.33,z]) .≈ SA[0.33,0.33,0.])
    # test normal and oriented surface
    @test all(WaterLilyPreCICE.normal(mesh[1]) .≈ [0.,0.,1.]) # should not throw
    @test all(WaterLilyPreCICE.dS(mesh[1]) .≈ [0.,0.,.5])
    @test WaterLilyPreCICE.area(mesh[1]) ≈ .5
    @test all(WaterLilyPreCICE.center(mesh[1]) .≈ [2/3,1/3,0.])
    # tests distance
    @test WaterLilyPreCICE.d²_fast(mesh[1], SA[0.,0.,0.]) ≈ sum(abs2,[2/3,1/3,0.]) # fast uses center
    @test WaterLilyPreCICE.d²(mesh[1], SA[0.,0.,0.]) ≈ 0.0 # this one is accurate
    @test WaterLilyPreCICE.d²(mesh[1], SA[0.,0.,z]) ≈ z^2 # this one is accurate
    @test WaterLilyPreCICE.d²(mesh[1], SA[0.,0.,2z]) ≈ 4z^2 # this one is accurate
    # the MeshBody
    b1 = MeshBody(mesh,(1,1),(x,t)->x,Rect(0.,0.,0.,1.,1.,1.),1.0,1/2,true)
    b2 = copy(b1)
    @test b1.mesh == b2.mesh
    @test all(b1.map(SA[0.,0.,0.],0.) .≈ b2.map(SA[0.,0.,0.],0.))
    @test sdf(b1,SA[2/3,1/3,0.],0.) ≈ 0. # in the bbox, accurate
    @test all(measure(b1,SA[2/3,1/3,0.],0.) .≈ (0., [0.,0.,1.],[0.,0.,0.])) # in the bbox, accurate
    # more complex outside of box
    # @test sdf(b1,SA[2/3,1/3,2.],0.0) ≈ 1.5 # outside the bbox, accurate
    @test all(measure(b1.mesh,SA[2/3,1/3,2.],0.) .≈ (2.0, [0.,0.,1.])) # measuring the mesh is always correct
    # @test all(measure(b1,SA[2/3,1/3,2.],0.) .≈ (1.5, [0.,0.,0.], [0.,0.,0.])) # outside the bbox, not accurate
    # the single trinagle has zero flux contribution as it is flat on the xy-plane
    @test all(volume(b1.mesh) .≈ volume(b1) .≈ [0.,0.,0.])
    # test the important bits
    L = 64 # sphere and domain size
    mesh = triangle_mesh(Sphere(Point3f(L/2,L/2,L/2), L/4)) # triangle mesh from a sphere
    sphere = MeshBody(mesh,(1,1),(x,t)->x,Rect(mesh.position),1.0,1/2,true)
    flow = Flow((L,L,L),(0.,0.,0.))
    apply!((x,t)->x[1],flow.p) # hydrostatic pressure
    # @test all(WaterLilyPreCICE.volume(sphere) .≈ 4/3*π*(L/4)^3)
    # @test all(sum(WaterLilyPreCICE.forces(sphere,flow)) .≈ [0.,0.,1.75])
end
@testset "utils.jl" begin
    # make a sim with a NoBody()
    sim = Simulation((16,16,16),(0.,0.,0.),16)
    s = WaterLilyPreCICE.Store(sim)
    # save the data
    WaterLilyPreCICE.store!(s,sim)
    sim_step!(sim)
    sim.flow.u .-= 123.0; sim.flow.p .-= 123.0
    WaterLilyPreCICE.revert!(s,sim)
    @test all(sim.flow.u .≈ 0.0) && all(sim.flow.p .≈ 0.0)
    # test the unpack function
    a = [0,1,0,1]
    @test all(WaterLilyPreCICE.unpack(a) .== [[0,1],[0,1]])
    a = [0,0,0,0.5,1,1,1,0,0,0,0.25,0.5,0.75,1,1,1]
    @test all(WaterLilyPreCICE.unpack(a) .== [[0,0,0,0.5,1,1,1],[0,0,0,0.25,0.5,0.75,1,1,1]])
    # test the knotVectorUnpack
    a = zeros(4,2); a[3,:] .= 1
    @test WaterLilyPreCICE.knotVectorUnpack(a) == [[0,1],[0,1]]
end
@testset "KDTree.jl" begin
    # test the box
    b = WaterLilyPreCICE.Bbox(SA[0.0,0.0,0.0],SA[0.5,0.5-eps(),0.5-eps()])
    @test b.C == SA[0.0,0.0,0.0] && b.R == SA[0.5,0.5-eps(),0.5-eps()]
    @test !b.leaf && b.indices == nothing
    # test the split of the box
    left,right = WaterLilyPreCICE.split(b)
    @test all(left.C.≈SA[-0.25,0.,0.]) && all(left.R.≈SA[0.25,0.5,0.5])
    @test all(right.C.≈SA[0.25,0.,0.]) && all(right.R.≈SA[0.25,0.5,0.5])
    # test inside
    @test !inside([0.5+eps(),0.5,0.5],b) && !inside([-0.5-eps(),0.5,0.5],b)
    # test the filer
    points = rand(3,10)
    box = WaterLilyPreCICE.bounding_box(points)
    left, right = WaterLilyPreCICE.split(box)
    # @test length(left) + length(right) >= size(points,2)
    # @test !(left in right) # check that there is no point in both boxes

    # test the tree
    # points = rand(3,100)
    # t = WaterLilyPreCICE.Tree(points)
    # @test length(t) > 1
    # @test t[1].C == SA[0.0,0.0,0.0] && t[1].R == SA[0.5,0.5-eps(),0.5-eps()]
    # @test !t[1].leaf && t[1].indices == nothing
end