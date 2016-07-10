import sympy as sm
from ngsolve import *
import netgen

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

k = 1
F = curlcurlEexact - k*k*Eexact


mesh = Mesh (netgen.csg.unit_cube.GenerateMesh(maxh=0.5))

# # exact electric field
# exu = CoefficientFunction( ( (1-y)*y*z*(1-z),
#                              (1-x)*x*(1-z)*z,
#                              (1-x)*x*(1-y)*y) )

# # curl curl 
# F = CoefficientFunction( (2*((1-z)*z + (1-y)*y),
#                           2*((1-z)*z + (1-x)*x),
#                           2*((1-x)*x + (1-y)*y) ) )
# k = 1 
# F = F - k*k* exu

p = 3;
V = HCurl(mesh, order=p+1,dirichlet=[1,2,3,4,5,6])
Q = HCurl(mesh, order=p,  flags={"orderinner": 1})
Y = HCurl(mesh, order=p+4,flags={"discontinuous": True} )

XY = FESpace( [V,Q,Y])

E,M,e = XY.TrialFunction()
G,W,d = XY.TestFunction()


def cross(G,N):
    return CoefficientFunction( ( G[1]*N[2] - G[2]*N[1],
                                  G[2]*N[0] - G[0]*N[2],
                                  G[0]*N[1] - G[1]*N[0] ) )

n = specialcf.normal(mesh.dim)

a = BilinearForm(XY, symmetric=False, flags={"eliminate_internal" : True})
a+= SymbolicBFI(curl(E) * curl(d) - k*k*E*d)
a+= SymbolicBFI(curl(e) * curl(G) - k*k*e*G)
a+= SymbolicBFI( M * cross(d,n), element_boundary=True)
a+= SymbolicBFI( cross(e,n) * W, element_boundary=True)
a+= SymbolicBFI( curl(e) * curl(d) + e*d )

f = LinearForm(XY)
f+= SymbolicLFI(F*d)

EMe = GridFunction(XY)

## Direct solve:
##
# a.Assemble(heapsize=int(1e8))
# f.Assemble()
# uqe.vec.data = a.mat.Inverse(XY.FreeDofs()) * f.vec
# Draw (uMe.components[0])

## Iterative solution:
##
# clocal = Preconditioner(a, type="local")
cdirect = Preconditioner(a, type="direct")
a.Assemble(heapsize=int(1e8))
f.Assemble()
bvp = BVP(bf=a, lf=f, gf=EMe, pre=cdirect).Do()

Draw (EMe.components[0])

print ("L2-error:",
       sqrt(Integrate((EMe.components[0]-Eexact)*(EMe.components[0]-Eexact),
                      mesh)))

