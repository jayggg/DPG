from ngsolve.la import InnerProduct
from math import sqrt

def pcg(A, B, b, x=None, tol=1.e-16, maxits=100, saveitfn=None):
    
    """Preconditioned Conjugate Gradient iterations for solving the
    B-preconditioned A-linear system

        B A x = B b 

    in the inv(B)-inner product.
    """

    if x == None: 
        x = b.CreateVector()     # if not given initial guess, 
        x[:] = 0.0               # set initial solution vector = 0
    r = b.CreateVector()         # r = residual
    p = b.CreateVector()         # p = search direction
    Br= b.CreateVector()         # for storing B*r and A*p 
    Ap= Br
    r.data  = b - A * x
    Br.data = B * r
    p.data  = Br
    rBr = [0,InnerProduct(r,Br)] # 2 consecutive residual B-norms
    
    for it in range(maxits):
        
        rBr[0]  = rBr[1]         # update residual's B-norm 
        Ap.data = A * p 
        pAp     = InnerProduct(p,Ap)
        alpha   = rBr[0]/pAp     # alpha = <r, B*r> / <p, A*p>
        x.data += alpha * p      # x = x + alpha * p 
        r.data += (-alpha) * Ap  # r1 = r0 - alpha * A*p
        
        Br.data = B * r          
        rBr[1]  = InnerProduct(r,Br)
        beta    = rBr[1]/rBr[0]  # beta = <r1, B*r1> / <r0, B*r0>
        p.data  = beta * p
        p.data += Br             # p = B*r1 + beta * p

        if saveitfn != None:
            saveitfn(x,it)

        if abs(rBr[1].imag) > 1.e-8 or abs(pAp.imag) > 1.e-8:
            print('\n*** rBr=', rBr, '\n*** pAp=', pAp)
            print('*** System not Hermitian!')
        else:
            if rBr[0].real * rBr[1].real < 0 :
                print('\n*** Preconditioner indefinite!')
                
        print('PCG%6d'%it, ': pAp=%12.9g %12.9g'
              %(pAp.real,rBr[1].real) )
        if abs(pAp) < tol or abs(rBr[1]) < tol:
            break;
            
    return x


