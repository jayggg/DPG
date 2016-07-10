############################################################
# Solve Helmholtz equation with impedance boundary condition
#
#      - Delta u - k*k u = f    on Omega
#       n.grad u - i*k u = g    on bdry
#
# using the DPG variational formulation 
#
#   Y(e;v)    +  b(u,q; v)    = F(v)
#   b(w,r; e) +  c(w,r ; u,q) = conj(G(w,r))
#
# Here conj denotes complex conjugate,
#
#  Y(e,v)       = (grad e, grad v) + k*k (e,v)
#  b(u,q; v)    = (grad u, grad v) - k*k*(u,v) - <<q.n, v>> 
#  c(u,q; w,r)  = - <q.n - ik u, r.n - ik w> 
#  F(v)         = (f,v)
#  G(w,r,w)     = - <g, r.n - ik w> 
#
# where
#
# (.,.)   is sum over all (complex) element L2 inner products,
# <<.,.>> is sum over all (complex) element boundary L2 inner products,
# <.,.>   is global bdry (complex) L2 inner product.
#
# The new terms are contained in the form c, which is a 
# Hermitian positive addition to the whole system, which 
# imposes the impedance condition q.n - ik u = g.
############################################################

shared = "../libDPG"
constant heapsize = 10000000

geometry = square4bdry.in2d
mesh = square4bdry4.vol.gz

constant one   = 1
constant minus = -1.0

# approx number of waves in the (unit-sized) domain 
constant nwav = 2

# propagation angle 
constant theta = (pi/16.0)

# wavenumber
constant k  = (2.0*pi*nwav)


constant kc = (k*cos(theta))
constant ks = (k*sin(theta))
constant ksqr = (k*k) 
constant minusksqr = (-k*k) 

# exact solution ex = exp(I*k.x)  for error computation
coefficient ex     (exp(I*(kc*x + ks*y)))     


# Note that bc markers for the given square geometry are:
#              3
#        +-----------+
#        |           |
#        |           |
#      4 |           | 2
#        |           |
#        +-----------+
#              1
#  bc=1:  g = grad ex . [0,-1] - ik ex 
#  bc=2:  g = grad ex . [ 1,0] - ik ex  
#  bc=3:  g = grad ex . [ 0,1] - ik ex  
#  bc=4:  g = grad ex . [-1,0] - ik ex  

coefficient minusg   # -g
(I*(ks+k)*ex), (-I*(kc-k)*ex), (-I*(ks-k)*ex), (I*(kc+k)*ex)

coefficient minusgik # -g * ik
(I*k*I*(ks+k)*ex), (-I*k*I*(kc-k)*ex), (-I*k*I*(ks-k)*ex), (I*k*I*(kc+k)*ex)


coefficient ik          (I*k)
coefficient minusik     (-I*k)


# Finite element spaces                      (p = 0,1,2,...)
fespace fs1 -type=l2ho   -order=4 -complex   # e, v, deg p+2
fespace fs2 -type=h1ho   -order=3 -complex   # u, w, deg p+1
fespace fs3 -type=hdivho -order=2            # q, r, deg p 
                         -orderinner=1 
			 -complex
fespace fs  -type=compound -spaces=[fs1,fs2,fs3] -complex


# Forms
#   RHS: We use standard NGSolve source integrators.
linearform lf  -fespace=fs
neumannhdiv   minusg     -comp=3    # -<g,r.n>
neumann       minusgik   -comp=2    # +<g,ik*wh> = -<g*ik,wh>
#   LHS: We use standard and DPG integrators to make the
#   composite sesquilinear form
#    a( e,u,q,uh ; v,w,r,wh )  
#          = Y(e;v) + b(u,q; v)
#                   + conj( b(w,r; e) +  c(w,r ; u,q) ).
bilinearform dpg -fespace=fs -linearform=lf 
       -nonsym      
       -eliminate_internal   
gradgrad    (2) (1) one              # (grad u, grad v)
eyeeye      (2) (1) minusksqr        # - k*k*(u,v) 
flxtrc      (3) (1) minus            # - <<q.n, v>> 
flxflxbdry  (3) (3) minus            # - <q.n, r.n>
flxtrcbdry  (3) (2) minusik          # + <q.n, ik w>
trctrcbdry  (2) (2) minusksqr        # - <ik u, ik w> 
laplace             one   -comp=1
mass                ksqr  -comp=1


# Solve iteratively:
gridfunction euq -fespace=fs

#preconditioner c  -type=direct -bilinearform=dpg 
#preconditioner c  -type=local -bilinearform=dpg 
#preconditioner c  -type=vertexschwarz -addcoarse  -bilinearform=dpg 
preconditioner c  -type=vertexschwarz  -bilinearform=dpg 

numproc bvp n2 -bilinearform=dpg -linearform=lf 
        -gridfunction=euq 
	-preconditioner=c
	-solver=cg 
	-innerproduct=hermitean 
	-prec=1.e-10 -maxsteps=1000

# Compute error
gridfunction uu -fespace=fs2 -addcoef 
numproc getcomp ngc3 -comp=2 -compoundgf=euq -componentgf=uu
coefficient err ( abs(uu-ex) )
numproc integrate absL2error -coefficient=err

# Visualize
numproc visualization see_real_part 
        -scalarfunction=euq.2:1 -subdivision=4  -nolineartexture 

