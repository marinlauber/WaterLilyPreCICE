using Plots
include("body_bvh_clean.jl")

# dummy measure all the things
function direct_measure(body::MeshBody,x::SVector{D},t;kwargs...) where D
    mesh = body.mesh
    u=1; a=b=d²_fast(mesh[1],x) # fast method
    for I in 2:length(mesh)
        b = d²_fast(mesh[I],x)
        b<a && (a=b; u=I) # Replace current best
    end
    n,p = normal(mesh[u]),SVector(locate(mesh[u],x))
    s = x-p # signed Euclidian distance
    d = sign(sum(s.*n))*√sum(abs2,s)
    return d,n,zero(x)
end

# main code, map is inverted compared to classical WL
map(x,t) = x.+SA[4,4,4]
body = MeshBody("bvh/obj/sphere_centered.obj";map,T=Float32,lvl=4,mem=Array);
# what's the bounding bo of the mesh
bb = Rect(body.mesh.position)

bvh = body.bvh;

N = 24
sdf = zeros(N,N,N)
measure_sdf!(sdf, body, 0.0)
@inside sdf[I] = direct_measure(body,loc(0,I),0.0)[1]
flood(sdf[:,:,N÷2])
contour!(sdf[:,:,N÷2]',levels=[0.], lw=5)