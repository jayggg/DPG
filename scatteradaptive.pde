# Use the DPG method to adaptively compute a wave scattered 
# from a triangular scatterer. Adaptivity puts relatively more
# elements where the beam-like scattered wave is present.
#
# One can use any of the different possible implementation 
# techniques for the Helmholtz equation with impedance bc.
# (The results were observed not to depend on which.)
#
################################################################
# Compute scattered wave from a triangular scatterer. 
# The problem for utotal = uincident + uscattered is 
#
#  -Delta utotal - k*k utotal = 0,    outside scatterer
#                      utotal = 0     on scatterer boundary 
#                                     (sound-soft b.c.).
#
# Given uincident, we compute the scattered wave by solving:
#
#  -Delta uscattered - k*k uscattered = 0,  outside scatterer,
#                 uscattered = -uincident,  on scatterer boundary,
#   n.grad uscattered - ik uscattered = 0,  on rest of boundary.
#
# Press the Solve button repeatedly to proceed with successive 
# adaptive iterations.
################################################################

shared = libDPG
constant heapsize = 1000000000

geometry = triangularscatterer.in2d
mesh = triangularscatterer.vol.gz

# Propagation angle 
constant theta = (pi/3.0)

# Wavenumber
constant k  = (5*pi)

# Sound-soft (Dirichlet) conditions by penalty
constant penalty = 1.e5

constant kc = (k*cos(theta))
constant ks = (k*sin(theta))
constant ksqr  = (k*k) 

# The incident wave
coefficient uincident   (exp(I*(kc*x + ks*y)))

constant one   = 1
constant minus = -1.0
constant dd    = 1.0
constant cc    = 1.0
constant minusksqr = (-k*k) 

coefficient ik          (I*k)
coefficient minusik     (-I*k)
coefficient cik         (cc*I*k)
coefficient minuscik    (-cc*I*k)
coefficient diksqr      ((dd*I*k)*(dd*I*k))
coefficient minusdiksqr (-(dd*I*k)*(dd*I*k))


# The boundary has 5 parts, the fifth being the scatterer boundary. 
# We impose outgoing impedance bc everywhere, except the scatterer 
# boundary, where sound-soft boundary condition is imposed. 

coefficient notsoft
(-1.0), (-1.0), (-1.0), (-1.0), 0

coefficient notsoftik
(-I*k), (-I*k), (-I*k), (-I*k), 0


coefficient notsoftksqr
(-k*k),(-k*k),(-k*k),(-k*k),0

# Dirichlet bc imposed on scatterer boundary (bc=5) by penalty
coefficient soft
0,0,0,0,penalty

coefficient inc  
0,0,0,0,-penalty * exp(I*(kc*x + ks*y))

# Finite element spaces                              (p = 0,1,2,...)
fespace fs1 -type=l2ho   -order=4 -complex   # e, v, deg p+2
fespace fs2 -type=l2ho   -order=3 -complex   # u, w, deg p+1
fespace fs3 -type=hdivho -order=2            # q, r, deg p 
                         -orderinner=1 -complex
fespace fs4 -type=h1ho   -order=3            # uh, wh, deg p+1
                         -orderinner=1 -complex
fespace fs  -type=compound -spaces=[fs1,fs2,fs3,fs4] -complex


# Forms 
#   RHS:
linearform lf  -fespace=fs
neumann     inc       -comp=4        # - penalty <uincident, v> 
#   LHS:
bilinearform dpg -fespace=fs -linearform=lf 
       -nonsym      
       -eliminate_internal   
gradgrad    (2) (1) one              # (grad u, grad v)
eyeeye      (2) (1) minusksqr        # - k*k*(u,v) 
flxtrc      (3) (1) minus            # - <<q.n, v>> 
trctrc      (4) (1) cik              # + <<c * ik uh, v>>
trctrc      (2) (1) minuscik         # - <<c * ik u,  v>>
trctrc      (4) (4) diksqr           # - <<d * ik uh, d * ik wh>> 
trctrc      (4) (2) minusdiksqr      # + <<d * ik uh, d * ik w >> 
trctrc      (2) (2) diksqr           # - <<d * ik u,  d * ik w >> 
flxflxbdry  (3) (3) notsoft          # - <q.n, r.n>
flxtrcbdry  (3) (4) notsoftik        # + <q.n, ik wh>
trctrcbdry  (4) (4) notsoftksqr      # - <ik uh, ik wh> 
robin               soft  -comp=4    # + penalty <u, v> 
laplace             one   -comp=1    
mass                ksqr  -comp=1

# Solve:
gridfunction euqf -fespace=fs

preconditioner c  -type=direct -bilinearform=dpg 
#preconditioner c  -type=vertexschwarz -bilinearform=dpg 
#preconditioner c  -type=local -bilinearform=dpg 

numproc bvp n2 -bilinearform=dpg -linearform=lf 
        -gridfunction=euqf
	-preconditioner=c
	-solver=cg 
	-innerproduct=hermitean 
	-prec=1.e-10 -maxsteps=1000


# Estimate error & Mark elements for local refinement
fespace dg0 -type=l2ho  -order=0
gridfunction eestim -fespace=dg0

numproc enormsc estimate_using_e_norm -bilinearform=dpg -fespace=fs 
        -solution=euqf 
	-estimator=eestim      # output element error estimate values 
	-yintegrators=[13,14]  # integrators forming Y-innerproduct 
	-yspaces=[1]           # the Y-space

numproc markelements mark_large_error_elements
        -error=eestim -minlevel=1 -factor=0.25 


# Visualize
numproc visualization see_real_part 
        -scalarfunction=euqf.2:1 -subdivision=4  -nolineartexture 


