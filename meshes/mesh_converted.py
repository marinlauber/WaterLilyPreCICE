import numpy as np
import argparse

def to_float(ls):
    return [float(ls[i]) for i in range(len(ls))]
def to_int(ls):
    return [int(ls[i]) for i in range(len(ls))]

def stl2inp(fname):
    vertex=[]; nodes_id=[]
    stl = open(fname, "r")
    print(stl)
    lines = stl.readlines()[1:-1]
    tris=np.zeros((int(len(lines)//7),3,3))
    connectivity=np.zeros((int(len(lines)//7),3))
    for i in range(len(lines)//7):
        for k in range(3):
            vertex_coords = to_float(lines[i*7+k+2].strip().split()[1:])
            vertex.append(vertex_coords)
            nodes_id.append(3*i+k)
            tris[i,k,:] = vertex_coords
            connectivity[i,k]=3*i+k
    unique= np.unique(np.array(vertex), axis=0)
    idx = np.arange(2,unique.shape[0]+2)
    stl.close()
    ID=[]
    for i in range(tris.shape[0]):
        tri=[]
        for j in range(3):
            d = np.argmin(np.sum((tris[i,j,:]-unique)**2, axis=1))
            tri.append(idx[d])
        ID.append(tri)
    ID=np.array(ID)
    mx = max(idx)+1
    inp = open(fname[:-3]+"inp", "w")
    inp.write("*NODE, NSET=Ninterface")
    for i in range(len(idx)):
        inp.write("\n%d, %.8f, %.8f, %.8f" % (idx[i],*unique[i,:]))
    inp.write("\n*ELEMENT, type=S3, ELSET=PLATE\n")
    for i in range(ID.shape[0]):
        inp.write("%d, %d, %d, %d\n" % (i+mx,*ID[i,:]))
    inp.close()
    return None

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename")
    stl2inp(parser.parse_args().filename)