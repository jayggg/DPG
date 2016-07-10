""" Primal DPG for 3D Maxwell equations without using libDPG """

import sympy as sm
from ngsolve import *
import netgen

# Set a "manufactured" exact solution and RHS using sympy:

X,Y,Z=sm.symbols('X Y Z')
Es = ( (1-Y) * Y*Z * (1-Z), (1-X) * X * (1-Z) * Z, (1-X) * X * (1-Y) * Y )

def smcurl(M):        # symbolic curl
    return ( sm.diff( M[2], Y) - sm.diff( M[1], Z),
             sm.diff( M[0], Z) - sm.diff( M[2], X),
             sm.diff( M[1], X) - sm.diff( M[0], Y) )

def str2coef(symbol): # symbolic string to ngsolve coefficient function 
    print(symbol)
    return CoefficientFunction(eval(str(symbol)
                                    .replace("X","x")
                                    .replace("Y","y")
                                    .replace("Z","z")))
Eexact = str2coef(Es)
curlcurlEexact = str2coef(smcurl(smcurl(Es)))

k = 1                 # wave number & RHS: 
F = curlcurlEexact - k*k*Eexact


# Compute numerical solution using DPG forms

mesh = Mesh (netgen.csg.unit_cube.GenerateMesh(maxh=0.5))

p = 3;
Xo= HCurl(mesh, order=p+1,dirichlet=[1,2,3,4,5,6],complex=True)
Xh= HCurl(mesh, order=p,  complex=True,flags={"orderinner": 1})
Y = HCurl(mesh, order=p+4,complex=True,flags={"discontinuous": True})
XY = FESpace([Xo,Xh,Y], flags={"complex":True})

E,M,e = XY.TrialFunction()
G,W,d = XY.TestFunction()

def cross(G,N):
    return CoefficientFunction( ( G[1]*N[2] - G[2]*N[1],
                                  G[2]*N[0] - G[0]*N[2],
                                  G[0]*N[1] - G[1]*N[0] ) )

n = specialcf.normal(mesh.dim)

# Set bilinear and linear forms using NGSpy's symbolic forms

a = BilinearForm(XY, symmetric=False, flags={"eliminate_internal" : True})
a+= SymbolicBFI(curl(E) * curl(d) - k*k*E*d)
a+= SymbolicBFI(curl(e) * curl(G) - k*k*e*G)
a+= SymbolicBFI(M * cross(d,n), element_boundary=True)
a+= SymbolicBFI(cross(e,n) * W, element_boundary=True)
a+= SymbolicBFI(curl(e) * curl(d) + e*d)


f = LinearForm(XY)
f+= SymbolicLFI(F*d)

EMe = GridFunction(XY)

# Solve the linear system

cdirect = Preconditioner(a, type="direct")
a.Assemble(heapsize=int(1e8))
f.Assemble()
bvp = BVP(bf=a, lf=f, gf=EMe, pre=cdirect).Do()

Draw (EMe.components[0])


# Compare with exact solution

print ("L2-error:",
       sqrt(Integrate((EMe.components[0]-Eexact)*(EMe.components[0]-Eexact),
                      mesh)))

