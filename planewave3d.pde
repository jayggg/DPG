# Solve Helmholtz equation with impedance boundary condition
#
#      - Delta u - k*k u = f    on Omega
#       n.grad u - i*k u = g    on bdry
#
# on a 3D domain. 


shared = libDPG
constant heapsize = 1000000000

geometry = cube6bc.geo
mesh     = cube6bc4.vol.gz


constant one   = 1
constant minus = -1.0
constant dd    = 1.0
constant cc    = 1.0

# number of waves in a unit-sized domain 
constant nwav = 2

# propagation angle 
constant theta = (pi/11.0)
constant phi   = (pi/3.0)

# wavenumber
constant k  = (2.0*pi*nwav)
constant k1 = (k*cos(theta)*sin(phi))
constant k2 = (k*sin(theta)*sin(phi))
constant k3 = (k*cos(phi))

constant ksqr  = (k*k) 
constant minusksqr = (-k*k) 

coefficient ik          (I*k)
coefficient minusik     (-I*k)
coefficient cik         (cc*I*k)
coefficient minuscik    (-cc*I*k)
coefficient diksqr      ((dd*I*k)*(dd*I*k))
coefficient minusdiksqr (-(dd*I*k)*(dd*I*k))


# exact solution ex = exp(I*k.x)  for error computation
coefficient ex     (exp(I*(k1*x + k2*y + k3*z)))     

# grad ex = I * [k1,k2,k3] * exp( I * [k1,k2,k3] . [x,y,z] )

# We need to set b.c. using g = n.grad ex - ik ex  on each face:
#   
#  bc=1:  x=0,  g = grad ex . [-1, 0, 0] - ik ex 
#                 = -i (k1 + k)  ex
#  bc=2:  y=0,  g = grad ex . [ 0,-1, 0] - ik ex  
#                 = -i (k2 + k) ex
#  bc=3:  z=0,  g = grad ex . [ 0, 0,-1] - ik ex  
#                 = -i (k3 + k) ex 
#  bc=4:  x=1,  g = grad ex . [ 1, 0, 0] - ik ex  
#                 =  i (k1 - k) ex
#  bc=5:  y=1,  g = grad ex . [ 0, 1, 0] - ik ex  
#                 =  i (k2 - k) ex
#  bc=6:  z=1,  g = grad ex . [ 0, 0, 1] - ik ex  
#                 =  i (k3 - k) ex
#


coefficient minusg   # -g
(I*(k1+k)*ex), (I*(k2+k)*ex), (I*(k3+k)*ex), (I*(k-k1)*ex), (I*(k-k2)*ex), (I*(k-k3)*ex)

coefficient minusgik # -g * ik
(I*k*I*(k1+k)*ex), (I*k*I*(k2+k)*ex), (I*k*I*(k3+k)*ex), (I*k*I*(k-k1)*ex), (I*k*I*(k-k2)*ex), (I*k*I*(k-k3)*ex)

coefficient f 
(0.0)


# finite element spaces                              (p = 0,1,2,...)
fespace fs1 -type=l2ho   -order=6 -complex   # e, v, deg p+2
fespace fs2 -type=h1ho   -order=5 -complex   # u, w, deg p+1
fespace fs3 -type=hdivho -order=4            # q, r, deg p 
                         -orderinner=1 
			 -complex

fespace fs  -type=compound -spaces=[fs1,fs2,fs3] -complex


# forms 
linearform lf     -fespace=fs
source        f          -comp=1 
neumannhdiv   minusg     -comp=3        # -<g,r.n>
neumann       minusgik   -comp=2        # +<g,ik*wh> = -<g*ik,wh>

# We need to make the sesquilinearform
#  a( e,u,q,uh ; v,w,r,wh ) = 
#   Y(e;v)  +  b(u,q,uh; v) + conj( b(w,r,wh; e) +  c(w,r,wh ; u,q,uh) )

bilinearform dpg -fespace=fs -linearform=lf 
       -nonsym      # This option is necessary for Hermitian too.
       -eliminate_internal   
gradgrad    (2) (1) one              # (grad u, grad v)
eyeeye      (2) (1) minusksqr        # - k*k*(u,v) 
flxtrc      (3) (1) minus            # - <<q.n, v>> 
flxflxbdry  (3) (3) minus            # - <q.n, r.n>
flxtrcbdry  (3) (2) minusik          # + <q.n, ik w>
trctrcbdry  (2) (2) minusksqr        # - <ik u, ik w> 
laplace             one   -comp=1
mass                ksqr  -comp=1


# solve:
gridfunction euqf -fespace=fs

#preconditioner c  -type=direct -bilinearform=dpg 
#preconditioner c  -type=local -bilinearform=dpg 
#preconditioner c  -type=vertexschwarz -addcoarse  -bilinearform=dpg 
preconditioner c  -type=vertexschwarz -bilinearform=dpg 

numproc bvp n2 -bilinearform=dpg -linearform=lf 
        -gridfunction=euqf -preconditioner=c
	-solver=cg -innerproduct=hermitean 
	-prec=1.e-10 -maxsteps=1000


gridfunction uu -fespace=fs2 -addcoef 
numproc getcomp ngc3 -comp=2 -compoundgf=euqf -componentgf=uu
coefficient err ( abs(uu-ex) )
numproc integrate absL2error -coefficient=err

# Visualize
numproc visualization see_real_part 
        -scalarfunction=euqf.2:1 -subdivision=4  -nolineartexture 

