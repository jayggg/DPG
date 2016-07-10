# Using adaptivity in DPG. A simple laplace equation example.

shared = "../libDPG"
constant heapsize = 10000000

geometry = square.in2d
mesh = square2.vol.gz

constant one = 1
constant minus = -1.0
coefficient lam (1+I)

# Source (exact solution unknown)
coefficient f       ( (1+I)*exp( -100.0*(x*x+y*y) ) )

# Compound finite element space:
fespace fs1 -type=h1ho   -order=3 -dirichlet=[1] -complex
fespace fs2 -type=hdivho -order=2 -orderinner=1  -complex
fespace fs3 -type=l2ho   -order=4                -complex
fespace fs  -type=compound -spaces=[fs1,fs2,fs3] -complex

# Forms:
bilinearform dpg -fespace=fs -eliminate_internal		 
gradgrad (1) (3) lam      # (grad u, grad v) + Hermitian transpose 
flxtrc   (2) (3) minus    # - << q.n, v >>   + Hermitian transpose
laplace  one -comp=3      # (grad e, grad v)
mass     one -comp=3      # (e,v) 

linearform lf    -fespace=fs
source      f    -comp=3

# Solve:
gridfunction uqe -fespace=fs
numproc bvp n2 -bilinearform=dpg -linearform=lf 
        -gridfunction=uqe -solver=direct

# Estimate & Mark
fespace dg0 -type=l2ho   -order=0

gridfunction eestim -fespace=dg0

numproc enormsc estimate_using_e_norm
        -solution=uqe -estimator=eestim -fespace=fs -bilinearform=dpg
	-yintegrators=[2,3] -yspaces=[3]

numproc markelements mark_large_error_elements
        -error=eestim -minlevel=1 -factor=0.5 


# Visualize imaginary part 
numproc visualization see_solution
        -scalarfunction=uqe.1:1 -subdivision=4  -nolineartexture
