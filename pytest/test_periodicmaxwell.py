from ngsolve import *
from ctypes import CDLL
import time
import numpy as np

libDPG = CDLL("../libDPG.so")

def setup():
    ngsglobals.msg_level = 0

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
    return p, k, bdry, kbdry, f, Eex, igxn, ikbarg

def solve_periodicmaxwell(p, k, bdry, kbdry, f, Eex,
                          igxn, ikbarg, symbolic):

    mesh = Mesh("../pde/periodiclayers.vol.gz")
    S0 = FESpace("hcurlho", mesh, order=p+3, complex=True,
                 discontinuous=True)
    S1 = FESpace("hcurlho_periodic", mesh, order=p, complex=True,
                 xends=[0,1], yends=[0,1] )
    S2 = FESpace("hcurlho_periodic", mesh, order=p+1, complex=True,
                 orderinner=0, xends=[0,1], yends=[0,1])
    S = FESpace( [S0,S1,S2], complex=True)

    e,E,M = S.TrialFunction()
    v,F,W = S.TestFunction()

    b = LinearForm(S)
    a = BilinearForm(S, symmetric=False, eliminate_internal=True)


    # If True, use NGSolve's symbolic forms.
    # Else, use libDPG's precompiled forms. 
    # (The solution should be the same in either case.)
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

    start = time.time()
    b.Assemble()
    end = time.time()
    print("time to assemble b, symbolic: ",symbolic,": ", end-start)

    eEM = GridFunction(S)
    c = Preconditioner(a, type="direct")
    SetHeapSize(int(5e8))
    
    start = time.time()
    a.Assemble()
    end = time.time()
    print("time to assemble a, symbolic: ",symbolic,": ", end-start)

    bvp = BVP(bf=a, lf=b, gf=eEM, pre=c).Do()

    Eh = eEM.components[1]
    err = sqrt(Integrate( (Eh - Eex) * Conj(Eh - Eex), mesh))
    return err

def test_periodicmaxwell():
    p, k, bdry, kbdry, f, Eex, igxn, ikbarg = setup()
    errSym = solve_periodicmaxwell(p, k, bdry, kbdry, f, Eex, 
                                   igxn, ikbarg, symbolic=True)

    errNonsym = solve_periodicmaxwell(p, k, bdry, kbdry, f, Eex,
                                      igxn, ikbarg, symbolic=False)
    assert np.abs(errSym) < 2.e-10
    assert np.abs(errNonsym) < 2.e-10


if __name__ == "__main__":
    test_periodicmaxwell()

