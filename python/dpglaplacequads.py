""" DPG for Laplacian using rectangular elements """

# For command line run:
# Terminal>>  python3 dpglaplacequads.py 

from ngsolve import *
from netgen.geom2d import SplineGeometry, unit_square
from math import pi
from numpy import log

ngsglobals.msg_level = 1

def solvedpg(p, h,
             U = 16*x*(1-x)*y*(1-y),
             minusDeltaU = 32*(y*(1-y)+x*(1-x))) :

    mesh = Mesh(unit_square.GenerateMesh(quad_dominated=True, maxh=h))
    
    X = H1(mesh, order=p+1, dirichlet=[1,2,3,4])
    Q = HDiv(mesh, order=p, orderinner=0)
    Y = L2(mesh, order=p+2)
    XY = FESpace( [X,Q,Y])

    u,q,e = XY.TrialFunction()
    w,r,d = XY.TestFunction()
    n = specialcf.normal(mesh.dim)

    a = BilinearForm(XY, symmetric=True, eliminate_internal=True)
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
    print ("L2-error:", l2err)
    return l2err

def test_solvedpg():    
    success = abs(solvedpg(p=2,h=0.125)) < 1.e-14
    assert success, "Biquadratic not exactly solved with p=2"

    
def hconvergencetable(errors,maxr):
    print("========================")
    print(" Mesh    Errors    Rate")
    print("-----------------------")
    rates = []
    for i in range(maxr):
        rates.append('  *  ')
    for i in range(1,maxr):
        if abs(errors[i])>1.e-15:
            rates[i]=format(log(errors[i-1]/errors[i])/log(2),'+5.2f')
    for i in range(maxr):
        print(" h/%-4d %8.2e   %s "%(pow(2,i+1), errors[i],rates[i]))
    print("========================")
    
def collecterrors(k,maxr):
    l2e = []
    for l in range(0,maxr):
        hl = 2**(-l)/2
        l2e.append(solvedpg(p = k, h = hl,
                            U = sin(pi*x)*sin(pi*y),
                            minusDeltaU = 2*pi*pi*sin(pi*x)*sin(pi*y)))
    return l2e


maxlevels = 5
er = collecterrors(4,maxlevels)
hconvergencetable(er,maxlevels)
