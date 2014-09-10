# Primal DPG method for real valued Laplace eq. in 2D
# 
##########################################################
# The BVP: 
#     -Delta u + u = f    on Omega 
#                u = 0    on bdry.
#
# The DPG weak formulation: 
#  
#   Y(e;v)       +  b(u,q; v)  = (f,v)
#   b(w,r; e)                  = 0
#
# where conj denotes complex conjugate,
#
#  Y(e,v)            = (grad e, grad v) + (u,v)
#  b(u,q,uh; v)      = (grad u, grad v) + (u,v) - <<q.n, v>> 
#
# and where
#
# (.,.)   is sum over all element L2 inner products,
# <<.,.>> is sum over all element boundary L2 inner products.
#
# The spaces: 
#     u, w     in H1
#     e, v     in broken H1
#     q, r     element boundary traces of H(div)
##########################################################


shared = libDPG
constant heapsize = 10000000

geometry = square.in2d
mesh = square2.vol.gz

constant one   = 1
constant minus = -1.0

# Exact solution for error computation
coefficient uex     ( x*(1-x)*y*(1-y) )
coefficient graduex ( ( (1-2*x)*y*(1-y), x*(1-x)*(1-2*y) ) )
coefficient f       ( 2*x*(1-x)+2*y*(1-y) + uex )

# DPG's compound finite element space of index p, where p=0,1,2,...
# has these component spaces:
fespace fs1 -type=h1ho   -order=4 -dirichlet=[1]   # p+1
fespace fs2 -type=hdivho -order=3 -orderinner=1    # p
fespace fs3 -type=l2ho   -order=5                  # p+2
fespace fs  -type=compound -spaces=[fs1,fs2,fs3]

# Forms: Specify a dpg integrator operating on component 
# spaces I and J as:   <integrator name> <I> <J> <coeff>.
bilinearform dpg -fespace=fs 
                 -symmetric           # real symmetric matrices
		 -eliminate_internal  # statically condense out 
gradgrad (1) (3) one      # (grad u, grad v) + transpose 
flxtrc   (2) (3) minus    # - << q.n, v >>   + transpose
eyeeye   (1) (3) one      # (u,v)            + transpose 
laplace  one -comp=3      # (grad e, grad v)
mass     one -comp=3      # (e,v) 


linearform lf    -fespace=fs
source      f    -comp=3

# Solve:  After static condensation, the DPG system is guaranteed 
# to be symmetric and positive definite, so we use CG to solve. 
gridfunction uqe -fespace=fs


preconditioner c  -type=local -bilinearform=dpg 
# # If you want to use a direct solver instead, use this:
#define preconditioner c  -type=direct -bilinearform=dpg 

numproc bvp n2 -bilinearform=dpg -linearform=lf 
        -gridfunction=uqe -preconditioner=c
	-solver=cg -innerproduct=hermitean 
	-prec=1.e-10 -maxsteps=1000

# Separate u and q components:
gridfunction u -fespace=fs1      # Flag "addcoef" allows it to be  
               -addcoef          # treated like a coefficient.
gridfunction q -fespace=fs2   -addcoef  	
numproc getcomp get_u  -comp=1 -compoundgf=uqe -componentgf=u 
numproc getcomp get_q  -comp=2 -compoundgf=uqe -componentgf=q


# Compute L2 and H1 errors in u:            Calculate || u - U ||.
coefficient errl2 (sqrt((u-uex)*(u-uex))) # Here we need u to be 
                                          # treatable as a coefficient.
numproc integrate  calcL2err -coefficient=errl2  
coefficient errh1 ( sqrt ((grad_u-graduex)*(grad_u-graduex)) )  
numproc integrate  calcH1err -coefficient=errh1  

# Compute approx H^(-1/2) norm of error in q:
fespace RT -type=hdivho -order=6   # Space to extend numerical traces,
gridfunction qRT -fespace=RT       # of degree >= deg(fs2).

bilinearform hdivipe -fespace=RT   # H(div) inner product.
divdivhdiv one
masshdiv   one

gridfunction qex -fespace=RT
fespace      dg0 -type=l2ho   -order=0
gridfunction qerrsqr -fespace=dg0

numproc setvalues interpolate_flux # Interp/proj exact Q=graduex.
	-gridfunction=qex 
	-coefficient=graduex
numproc setvalues interp_q_to_qRT  
	-gridfunction=qRT          # Interp/proj q to RT space (here 
	-coefficient=q             # we need q treated as coefficient).

numproc fluxerr  calc_fluxerror_fracnorm  # Calculate ||q - Q||.
	-exactq=qex -discreteq=qRT -extensionspace=RT 
	-fespace=fs -hdivproduct=hdivipe -errorsquareq=qerrsqr  

# Write error values into file (one line per refinement & solve):
numproc writefile tablulate_errors
	-variables=[mesh.levels,fes.fs.ndof,fes.fs1.ndof,integrate.calcL2err.value,integrate.calcH1err.value,fluxerr.calc_fluxerror_fracnorm.value] 
	-filename=tablerror.out

# Visualize
numproc visualization see_real_part 
        -scalarfunction=uqe.1 -subdivision=4  -nolineartexture 
