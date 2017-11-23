""" Subdivide the unit cube into a mesh of cubes """

from netgen.meshing import Mesh, Element2D, Element3D, MeshPoint, Pnt
from netgen.meshing import FaceDescriptor
import ngsolve as ngs


def GenerateCubeMesh(n):  # create a n x n x n mesh of cubes

    mesh = Mesh()
    pids = []
    for i in range(n+1):
        for j in range(n+1):
            for k in range(n+1):
                pids.append(mesh.Add(MeshPoint(Pnt(i/n, j/n, k/n))))

    for i in range(n):
        for j in range(n):
            for k in range(n):
                base = i * (n+1)*(n+1) + j*(n+1) + k
                baseup = base+(n+1)*(n+1)
                pnum = [base, base+1, base+(n+1)+1, base+(n+1),
                        baseup, baseup+1, baseup+(n+1)+1, baseup+(n+1)]
                elpids = [pids[p] for p in pnum]
                mesh.Add(Element3D(1, elpids))

    def AddSurfEls(p1, dx, dy, facenr):
        for i in range(n):
            for j in range(n):
                base = p1 + i*dx+j*dy
                pnum = [base, base+dx, base+dx+dy, base+dy]
                elpids = [pids[p] for p in pnum]
                mesh.Add(Element2D(facenr, elpids))

    # Face X = 0
    mesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    AddSurfEls(0, 1, n+1, facenr=1)

    # Face Y = 0
    mesh.Add(FaceDescriptor(surfnr=2, domin=1, bc=2))
    AddSurfEls(0, (n+1)*(n+1), 1, facenr=2)

    # Face Z = 0
    mesh.Add(FaceDescriptor(surfnr=3, domin=1, bc=3))
    AddSurfEls(0, n+1, (n+1)*(n+1), facenr=3)

    # Face X = 1
    mesh.Add(FaceDescriptor(surfnr=4, domin=1, bc=4))
    AddSurfEls((n+1)**3-1, -(n+1), -1, facenr=4)

    # Face Z = 1
    mesh.Add(FaceDescriptor(surfnr=5, domin=1, bc=5))
    AddSurfEls((n+1)**3-1, -(n+1)*(n+1), -(n+1), facenr=5)

    # Face Y = 1
    mesh.Add(FaceDescriptor(surfnr=6, domin=1, bc=6))
    AddSurfEls((n+1)**3-1, -1, -(n+1)*(n+1), facenr=6)

    # Set names for the above boundary parts (0-based indexing!)
    mesh.SetBCName(0, 'X0')
    mesh.SetBCName(1, 'Y0')
    mesh.SetBCName(2, 'Z0')
    mesh.SetBCName(3, 'X1')
    mesh.SetBCName(4, 'Z1')
    mesh.SetBCName(5, 'Y1')

    return mesh


mesh = GenerateCubeMesh(9)
ngs.Draw(ngs.Mesh(mesh))
