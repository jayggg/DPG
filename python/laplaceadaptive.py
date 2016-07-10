""" Automatic adaptivity for Poisson equation using DPG method """

# Run using:
# Terminal>>   netgen laplaceadaptive.py


from ngsolve import *
from netgen.geom2d import SplineGeometry

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


def SolveBVP():
    XY.Update()
    uqe.Update()
    a.Assemble()
    b.Assemble()
    bvp.Do()
    Draw(uqe.components[0])
    Redraw(blocking=True)


def CalcErrorMark():
    e = uqe.components[2]
    elerr = Integrate(e*Conj(e) + grad(e)*Conj(grad(e)),
                      mesh, VOL, element_wise=True)
    maxerr = max(abs(elerr.NumPy()))
    
    for el in mesh.Elements():
        mesh.SetRefinementFlag(el, bool(abs(elerr[el.nr]) > 0.5*maxerr))

    return elerr

def adaptivestep():
    
    print("Adaptive step")
    SolveBVP()
    elerr = CalcErrorMark()
    globalerr = sqrt(abs(sum(elerr)))
    print("Total estimated error = ", globalerr)
    mesh.Refine()
    return globalerr
    

globalerr = 1
itcount = 0

while XY.ndof<30000 and globalerr > 1.e-6:

    itcount += 1
    mesh.Refine()
    print("Adaptive step ", itcount)

    SolveBVP()
    elerr = CalcErrorMark()
    globalerr = sqrt(abs(sum(elerr)))
    print("Total estimated error=%g, ndofs=%d"%(globalerr,XY.ndof))
    
