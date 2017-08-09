""" Automatic adaptivity for Poisson equation using DPG method """

from ngsolve import *
from netgen.geom2d import SplineGeometry

def setup():
    ngsglobals.msg_level=1

    geom = SplineGeometry("../pde/square.in2d")
    mesh = Mesh( geom.GenerateMesh(maxh=0.5))


    lam = 1+1j
    f = CoefficientFunction( lam*exp(-100.0*(x*x+y*y)) )

    p = 3

    Xo = H1(mesh,order=p+1,dirichlet=[1],complex=True)
    Xf = HDiv(mesh,order=p,complex=True,flags={"orderinner":1})
    Y  = L2(mesh,order=p+2,complex=True)
    XY = FESpace([Xo,Xf,Y], flags={"complex":True})

    u,q,e = XY.TrialFunction()
    w,r,d = XY.TestFunction()

    n = specialcf.normal(mesh.dim)

    a = BilinearForm(XY,symmetric=True,flags={"eliminate_internal" : True})
    a+= SymbolicBFI(grad(u) * grad(d) + grad(e) * grad(w))
    a+= SymbolicBFI(q*n*d, element_boundary=True)
    a+= SymbolicBFI(e*r*n, element_boundary=True)
    a.components[2] += Mass(1.0)
    a.components[2] += Laplace(1.0)

    b = LinearForm(XY)
    b+= SymbolicLFI(f*d)

    uqe = GridFunction(XY)
    c = Preconditioner(a, type="direct")
    bvp = BVP(bf=a, lf=b, gf=uqe, pre=c)
    return mesh, XY, uqe, a, b, bvp

def SolveBVP(XY, uqe, a, b, bvp):
    XY.Update()
    uqe.Update()
    a.Assemble()
    b.Assemble()
    bvp.Do()


def CalcErrorMark(mesh, uqe):
    e = uqe.components[2]
    elerr = Integrate(e*Conj(e) + grad(e)*Conj(grad(e)),
                      mesh, VOL, element_wise=True)
    maxerr = max(abs(elerr.NumPy()))
    
    for el in mesh.Elements():
        mesh.SetRefinementFlag(el, bool(abs(elerr[el.nr]) > 0.5*maxerr))

    return elerr

def adaptivestep(XY, uqe, a, b, bvp):
    SolveBVP(XY, uqe, a, b, bvp)
    elerr = CalcErrorMark()
    globalerr = sqrt(abs(sum(elerr)))
    mesh.Refine()
    return globalerr
    

def test_laplaceadaptive():
    globalerr = 1
    itcount = 0
    mesh, XY, uqe, a, b, bvp = setup()
    while XY.ndof<30000 and globalerr > 1.e-6:

        itcount += 1
        mesh.Refine()
        SolveBVP(XY, uqe, a, b, bvp)
        elerr = CalcErrorMark(mesh,uqe)
        globalerr = sqrt(abs(sum(elerr)))
        
    assert itcount <= 9, "More iterations required than expected"
    assert globalerr <= 1.e-6, "Error higher than expected"

if __name__ == "__main__":
    test_laplaceadaptive()
