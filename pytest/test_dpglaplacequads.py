""" DPG for Laplacian using rectangular elements """

from ngsolve import *
from netgen.geom2d import SplineGeometry, unit_square
from math import pi
from numpy import log

def solvedpg(p, h,
             U = 16*x*(1-x)*y*(1-y),
             minusDeltaU = 32*(y*(1-y)+x*(1-x))) :

    mesh = Mesh(unit_square.GenerateMesh(quad_dominated=True, maxh=h))
    
    X = H1(mesh, order=p+1, dirichlet=[1,2,3,4])
    Q = HDiv(mesh, order=p, flags={"orderinner": 0})
    Y = L2(mesh, order=p+2)
    XY = FESpace( [X,Q,Y])

    u,q,e = XY.TrialFunction()
    w,r,d = XY.TestFunction()
    n = specialcf.normal(mesh.dim)

    a = BilinearForm(XY, symmetric=True, flags={"eliminate_internal" : True})
    a+= SymbolicBFI(grad(u) * grad(d))
    a+= SymbolicBFI(grad(e) * grad(w))
    a+= SymbolicBFI(q*n*d, element_boundary=True)
    a+= SymbolicBFI(e*r*n, element_boundary=True)
    a.components[2] += Laplace(1.0) 
    a.components[2] += Mass(1.0)    

    f = LinearForm(XY) 
    f.components[2] += Source(minusDeltaU)

    c = Preconditioner(a, type="direct")
    uqe = GridFunction(XY)

    a.Assemble()
    f.Assemble()
    bvp = BVP(bf=a, lf=f, gf=uqe, pre=c).Do()

    Uh = uqe.components[0]
    l2err = sqrt(Integrate((Uh - U)*(Uh-U), mesh))
    return l2err

def in_range(err,expected_errors,k,l):
    eerr = expected_errors[k][l]
    return err <= eerr * 1.01 and err >= eerr * 0.99

def test_dpglaplacequads():
    ngsglobals.msg_level = 1
    expected_errors = [ [1.19e-2,   1.61e-3,    2.05e-4,    2.57e-5],
                        [1.91e-3,   1.24e-4,    7.81e-6,    4.89e-7],
                        [7.42e-5,   2.22e-6,    6.84e-8,    2.13e-9],
                        [8.70e-6,   1.40e-7,    2.2e-9,     3.44e-11]]

    maxlevel = 4
    maxorder = 4
    for k in range(1,maxlevel+1):
        for l in range(0,maxorder):
            print(k,l)
            hl = 2**(-l)/2
            err = solvedpg(p = k, h = hl,
                                U = sin(pi*x)*sin(pi*y),
                                minusDeltaU = 2*pi*pi*sin(pi*x)*sin(pi*y))
            assert in_range(err,expected_errors,k-1,l) #, "Error outside expected range for p=%d, k=%d" % (k,l)

if __name__ == "__main__":
    test_dpglaplacequads()
