using Test,WaterLilyPreCICE,StaticArrays

@testset "Utils" begin
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
import LinearAlgebra: cross
@testset "MeshBodies" begin
    for T in [Float32,Float64]
        # inside bbox or not
        rect = Rect(0,0,0,1,1,1) # origin and widths
        @test WaterLilyPreCICE.outside(SA{T}[0.5,1,2.5],rect) &&
             !WaterLilyPreCICE.outside(SA{T}[0.5,0.5,0.5],rect)
        # distance to box center
        rect = Rect(0,0,0,1,1,1) # origin and widths
        @test WaterLilyPreCICE.outside(SA{T}[-1,0,0],rect)
        @test !WaterLilyPreCICE.outside(SA{T}[ 1,0,0],rect)
        # make a single tri mesh
        mesh = GeometryBasics.Mesh(Point{3, T}.([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.]]),
                                   [TriangleFace{Int}(1, 2, 3)])
        # test tri operation
        tri = mesh[1];
        @test all(WaterLilyPreCICE.locate(tri, SA{T}[0.,0.,0.]) .≈ [0.,0.,0.])
        @test all(WaterLilyPreCICE.locate(tri, SA{T}[1.,1.,1.]) .≈ [1.0,1.0,0.0])
        @test all(WaterLilyPreCICE.locate(tri, SA{T}[0.,1.,0.]) .≈ [.5,.5,0.])
        z = 5-10rand() # anywhere above should project onto it
        @test all(WaterLilyPreCICE.locate(tri, SA{T}[0.33,0.33,z]) .≈ [0.33,0.33,0.])
        # test normal and oriented surface
        @test all(WaterLilyPreCICE.normal(tri) .≈ [0.,0.,1.]) # should not throw
        @test all(WaterLilyPreCICE.dS(tri) .≈ [0.,0.,.5])
        @test WaterLilyPreCICE.area(tri) ≈ .5
        @test WaterLilyPreCICE.area(tri) ≈ √sum(abs2,WaterLilyPreCICE.dS(tri))
        @test WaterLilyPreCICE.area(tri) ≈ 0.5f0*√sum(abs2,cross(tri.points[2]-tri.points[1],tri.points[3]-tri.points[1]))
        @test all(WaterLilyPreCICE.center(tri) .≈ [2/3,1/3,0.])
        # tests distance
        @test WaterLilyPreCICE.d²_fast(tri, SA{T}[0.,0.,0.]) ≈ sum(abs2,[2/3,1/3,0.]) # fast uses center
        @test WaterLilyPreCICE.d²(tri, SA{T}[0.,0.,0.]) ≈ 0.0 # this one is accurate
        @test WaterLilyPreCICE.d²(tri, SA{T}[0.,0.,z]) ≈ z^2 # this one is accurate
        @test WaterLilyPreCICE.d²(tri, SA{T}[0.,0.,2z]) ≈ 4z^2 # this one is accurate
        # test quad operation
        quad_points = Point{3, T}.([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.]])
        quad_face = [QuadFace{Int}(1, 2, 3, 4)]
        quad = GeometryBasics.Mesh(quad_points,quad_face)[1]
        @test all(WaterLilyPreCICE.locate(quad, SA{T}[0.,0.,0.]) .≈ [0.,0.,0.])
        @test all(WaterLilyPreCICE.locate(quad, SA{T}[1.,1.,0.]) .≈ [1.,1.,0.])
        @test all(WaterLilyPreCICE.locate(quad, SA{T}[.5,.5,1.]) .≈ [.5,.5,0.])
        # test normal and oriented surface
        @test all(WaterLilyPreCICE.normal(quad) .≈ [0.,0.,1.]) # should not throw
        @test all(WaterLilyPreCICE.dS(quad) .≈ [0.,0.,1.])
        @test all(WaterLilyPreCICE.subdS(quad) .≈ ([0.,0.,.5],[0.,0.,.5]))
        @test WaterLilyPreCICE.area(quad) ≈ 1
        @test WaterLilyPreCICE.area(quad) ≈ √sum(abs2,WaterLilyPreCICE.dS(quad))
        @test all(WaterLilyPreCICE.center(quad) .≈ [.5,.5,0.])
        # tests distance
        @test WaterLilyPreCICE.d²_fast(quad, SA{T}[0.,0.,0.]) ≈ sum(abs2,[0.5,0.5,0.]) # fast uses center
        @test WaterLilyPreCICE.d²(quad, SA{T}[0.,0.,0.]) ≈ 0.0 # this one is accurate
        @test WaterLilyPreCICE.d²(quad, SA{T}[0.,0.,z]) ≈ z^2 # this one is accurate
        @test WaterLilyPreCICE.d²(quad, SA{T}[0.,0.,2z]) ≈ 4z^2 # this one is accurate
        # the MeshBody
        b1 = MeshBody(mesh;map=(x,t)->x,scale=1.0,T)
        b2 = copy(b1)
        @test b1.mesh == b2.mesh     # are the same
        @test !(b1.mesh === b2.mesh) # not the same object
        @test all(b1.map(SA{T}[0.,0.,0.],0.) .≈ b2.map(SA{T}[0.,0.,0.],0.))
        @test sdf(b1,SA{T}[2/3,1/3,0.],0.) ≈ zero(T) # in the bbox, accurate
        @test all(measure(b1,SA{T}[2/3,1/3,0.],0.) .≈ (0., [0.,0.,1.],[0.,0.,0.])) # in the bbox, accurate
        # more complex outside of box
        @test sdf(b1,SA{T}[2/3,1/3,4.1],0.0) ≈ 8 # outside the bbox, not accurate
        @test all(measure(b1.mesh,b1.velocity,SA{T}[2/3,1/3,4.1]) .≈ (4.1, [0.,0.,1.], [0.,0.,0.])) # measuring the mesh is always correct
        @test all(measure(b1,SA{T}[2/3,1/3,4.1],0.) .≈ (8.f0, [0.,0.,0.], [0.,0.,0.])) # outside the bbox, not accurate
        # the single triangle has zero flux contribution as it is flat on the xy-plane
        @test all(volume(b1.mesh) .≈ volume(b1) .≈ [0.,0.,0.])
        # test velocity update
        b3 = MeshBody(mesh;map=(x,t)->x.+t,scale=1.0,T)
        update!(b3,one(T),one(T)) # provides time and dt such that the map is x.+1
        @test all(measure(b3,SA{T}[5/3,4/3,1],0.) .≈ (0., [0.,0.,1.], [1.,1.,1.]))
        @test all(measure(b3,SA{T}[5/3,4/3,0],0.) .≈ (-1., [0.,0.,1.], [1.,1.,1.]))
        @test all(measure(b3,SA{T}[5/3,4/3,5.1],0.) .≈ (8.f0, [0.,0.,0.], [0.,0.,0.])) # outside the bbox, not accurate
    end
end
@testset "GismoInterface" begin
    Qs = collect(0:0.1:1.0)
    @test Qs==Qs
end
@testset "Interface" begin
    # @test_nowarn
    # # clean-up
    # @test_nowarn rm("TEST_DIR",recursive=true)
    # @test_nowarn rm("test_vtk_reader_$D.pvd")
    for T in [Float32,Float64]
        # test the important bits on tris
        sphere = MeshBody("../meshes/sphere.stl";T)
        numberOfVertices = length(sphere.mesh.position)
        interface = Interface([one(T)],[zero(T),zero(T)],zeros(T,numberOfVertices,3),[[1,2,3]],[T(0.25)],
                              "rw_mesh","read_data","write_data")
        @test_nowarn update!(sphere,interface,one(T))
        @test all(interface.deformation .≈ zero(T))
        @test_broken update!(interface)
    end
end
@testset "WaterLilyPreCICE" begin
    for T in [Float32,Float64]
        # test the important bits on tris
        L = 64 # sphere and domain size, the radius is L/4
        sphere = MeshBody("../meshes/sphere.stl";scale=L/2,map=(x,t)->x.+L/2.f0,T)
        flow = Flow((L,L,L),(0.,0.,0.);T)
        apply!((x,t)->x[1],flow.p) # hydrostatic pressure
        @test all(isapprox.(WaterLilyPreCICE.volume(sphere)./(L/4)^3,4/3*π,atol=0.06))
        @test all(isapprox.(abs.(sum(WaterLilyPreCICE.forces(sphere,flow)))/(4π/3*(L/4)^3),[0.,0.,1.],atol=5e-2))
        # test the important bits on quad
        cube = MeshBody("../meshes/cube_S4.inp";scale=L/2,map=(x,t)->x.+L/2.f0,T)
        @test all(isapprox.(WaterLilyPreCICE.volume(cube)./(L/2)^3,1,atol=1e-4))
        @test all(isapprox.(abs.(sum(WaterLilyPreCICE.forces(cube,flow,δ=0)))/(L/2)^3,[0.,0.,1.],atol=1e-4))
    end
end