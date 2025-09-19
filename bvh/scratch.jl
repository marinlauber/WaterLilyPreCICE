using GeometryBasics
using MeshIO, FileIO

mesh=load("obj/dragon_15k.obj")


Rec=Rect(mesh)


function split_longest(rect::Rect)
    axis = argmax(rect.widths)  # 1=x, 2=y, 3=z
    value = rect.origin[axis] + rect.widths[axis] / 2
    return GeometryBasics.split(rect, axis, value)
end


function change(x)
    x[4]=8
end