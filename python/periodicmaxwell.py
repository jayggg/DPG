""" Periodic Maxwell system, solved using libDPG:

TEST PROBLEM DETAILS:

We want to solve the boundary value problem 

 curl ( curl E ) - k * k * E = f,   on Omega = [0, 1]^3,

with boundary condition 

 n x ( curl E ) + i * k * n x (E x n) = g, on {z=0} and {z=1},

and periodic b.c. elsewhere.  (Here we are using the -i*om*t 
time-harmonic convention.) Introducing the flux interface variable
M = -i * curl(E), the DPG variational formulation reads

   Y(e;v)          +  b(E,M; v)    = L(v),      for all v
   conj(b(F,W; e)) +  c(E,M ; F,W) = G(F,W),    for all F, W

where,

  Y(e;v)       = (curl e, curl v) + (e,v)
  b(E,M; v)    = (curl E, curl v) - (k*k*E,v) + i*<<n x M, v>> 
  c(E,M; F,W)  = - <n x M + k * n x (E x n), n x W + k * n x (F x n)> 
  L(v)         = (f,v)
  G(F,W)       = - <g, i * (n x W + k * n x (F x n)) > 
               = <i * g, n x W + k * n x (F x n) >. 

To measure errors, we consider the case with the contrived exact solution:

 E = [ z         ] 
     [ z * z     ]
     [ z * z * z ]
 

                                 {z=1} boundary, bc=3  
                                 _____________________
         Omega = [0, 1]^3       /                    /|     
     z                         /____________________/ |     
     ^                         |                    | |      
     |                         |  material 2        | |  all   
     |                         |                    | |  periodic 
     |                         |  wavenumber k2     |/|  boundaries, 
     |------>y        {z=0.5}, |____________________/ |  bc=1   
    /                 bc=4     |                    | |     
   /                           |  material 1        | |      
  x                            |                    | |     
                               |  wavenumber k1     |/      
                               |____________________/       
                                          
                                {z=0} boundary, bc=2   

"""

from ngsolve import *
from ctypes import CDLL

libDPG = CDLL("../libDPG.so")
mesh = Mesh("../pde/periodiclayers.vol.gz")


symbolic = False   # If True, use NGSolve's symbolic forms.
                   # Else, use libDPG's precompiled forms. 
                   # (The solution should be the same in either case.)

k1 = 1.0           # bot material 
k2 = 10.0          # top material

k  = CoefficientFunction( [k1, k2] )
bdry  = CoefficientFunction( [0.0, 1.0, 1.0, 0.0] )
kbdry = CoefficientFunction( [0.0, k1,  k2,  0.0] )

f  = CoefficientFunction( (-k*k*z, -2-k*k*z*z, -k*k*z*z*z) )
Eex= CoefficientFunction( ( z,      z*z,        z*z*z    ) )
igxn = CoefficientFunction([( 0j           , 0j      , 0j ),
                            ( 0            , 1j      , 0 ),
                            ( -2j - k2     , 1j + k2 , 0 ),
                            ( 0            , 0       , 0 )])
ikbarg=CoefficientFunction([(0            ,  0                 , 0),
                            ( 1j * k1     ,  0                 , 0),
                            (-k2 * (1j+k2), -k2*(2*1j+k2)      , 0),
                            (0            ,  0                 , 0)])


p = 3 
S0 = FESpace("hcurlho", mesh, order=p+3, complex=True,
             flags={"discontinuous":True})
S1 = FESpace("hcurlho_periodic", mesh, order=p, complex=True,
             flags={'xends':[0,1], 'yends':[0,1] })
S2 = FESpace("hcurlho_periodic", mesh, order=p+1, complex=True,
             flags={"orderinner": 0, 'xends':[0,1], 'yends':[0,1]})
S = FESpace( [S0,S1,S2], flags={"complex":True})

e,E,M = S.TrialFunction()
v,F,W = S.TestFunction()


b = LinearForm(S)
a = BilinearForm(S, symmetric=False, flags={"eliminate_internal" : True})


if symbolic:
    
    b+= SymbolicLFI( f * v )                   # (f,v)
    b+= SymbolicLFI(ikbarg * F.Trace(), BND)   # <i*kbar* n x (g x n), F>
    b+= SymbolicLFI(igxn * W.Trace(),BND)      # <i * g x n, W>

    def cross(G,N):    # G x N 
        return CoefficientFunction( ( G[1]*N[2] - G[2]*N[1],
                                      G[2]*N[0] - G[0]*N[2],
                                      G[0]*N[1] - G[1]*N[0] ) )    
    def bvol( EE,  vv): 
        return ( curl(EE) * curl(vv) - k*k* EE * vv )
    
    n = specialcf.normal(mesh.dim)
    
    a+= SymbolicBFI( bvol(E,v) +                 # (curl E, curl v) - k*k(E,v)
                     Conj(bvol(e,F)) )           # + c.c.
    a+= SymbolicBFI( 1j * M * cross(v,n) - 1j * cross(e,n) * W,                 
                 element_boundary=True)          # i<<M, v x n>>           
    a+= SymbolicBFI(-kbdry * Conj(kbdry) * E.Trace() *              
                    F.Trace(), BND)              # -<k*kbar E x n, F x n>
    a+= SymbolicBFI(-bdry*M.Trace()*W.Trace(),   # -<M x n, W x n>
                    BND)
    a+= SymbolicBFI(kbdry * E.Trace() * cross(W.Trace(), n) +
                    Conj(kbdry) * cross( M.Trace(),n) * F.Trace(),
                    BND)                         # <k E, W x n> + c.c.
    a+= SymbolicBFI(e*v + curl(e)*curl(v))

else:

    b.components[0] += LFI("sourceedge", coef=f)
    b.components[2] += LFI("neumannedge", coef=igxn )
    b.components[1] += LFI("neumannedge", coef=ikbarg)

    a+= BFI("curlcurlpg", coef=[2,1,1])          # (curl E, curl v)
    a+= BFI("eyeeyeedge", coef=[2,1,-k*k])       # -(k*k E, v)
    a+= BFI("trctrcxn",   coef=[3,1,1j] )        # i<<M, v x n>>
    a+= BFI("xnbdry", coef=[2,3,kbdry])          # <k E, W x n>
    a.components[1]+=BFI("robinedge", coef=-kbdry * Conj(kbdry))
    a.components[2]+=BFI("robinedge", coef=-bdry)
    a.components[0]+= BFI("massedge", coef=1.0)
    a.components[0]+= BFI("curlcurledge", coef=1.0)

b.Assemble()

eEM = GridFunction(S)
c = Preconditioner(a, type="direct")
a.Assemble(heapsize=int(5e8))

bvp = BVP(bf=a, lf=b, gf=eEM, pre=c).Do()


Eh = eEM.components[1]
Draw(Eh)
print("Error = ", sqrt(Integrate( (Eh - Eex) * Conj(Eh - Eex), mesh)))
