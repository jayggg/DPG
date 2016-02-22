# curl ( curl E ) - k * k * E = f,   on Omega = [0, 1]^3
#
# n x ( curl E ) + i * k * n x (E x n) = g, on {z=0} and {z=1},
#
# and periodic b.c. elsewhere.
#
# (Using -i*om*t time-harmonic convention.)
#
# Flux variable M = -i * curl(E).
#
# DPG variational formulation: 
#
#   Y(e;v)          +  b(E,M; v)    = L(v),      for all v
#   conj(b(F,W; e)) +  c(E,M ; F,W) = G(F,W),    for all F, W
#
# where,
#
#  Y(e;v)       = (curl e, curl v) + (e,v)
#  b(E,M; v)    = (curl E, curl v) - (k*k*E,v) + i*<<n x M, v>> 
#  c(E,M; F,W)  = - <n x M + k * n x (E x n), n x W + k * n x (F x n)> 
#  L(v)         = (f,v)
#  G(F,W)       = - <g, i * (n x W + k * n x (F x n)) > 
#               = <i * g, n x W + k * n x (F x n) > 
#
# and
#
# Exact solution 
# E = [ z         ] 
#     [ z * z     ]
#     [ z * z * z ]
# 
#
#                                 {z=1} boundary, bc=3  
#                                 _____________________
#         Omega = [0, 1]^3       /                    /|     
#     z                         /____________________/ |     
#     ^                         |                    | |      
#     |                         |  material 2        | |  all   
#     |                         |                    | |  periodic 
#     |                         |  wavenumber k2     |/|  boundaries, 
#     |------>y        {z=0.5}, |____________________/ |  bc=1   
#    /                 bc=4     |                    | |     
#   /                           |  material 1        | |      
#  x                            |                    | |     
#                               |  wavenumber k1     |/      
#                               |____________________/       
#                                          
#                                {z=0} boundary, bc=2   
#

shared = libDPG
constant heapsize = 500000000

geometry = periodiclayers.geo       
mesh     = periodiclayers.vol.gz


constant k1 = (1.0)
constant k1bar= (1.0)

constant k2 = (10.0)
constant k2bar = (10.0)


# Source
coefficient f 
( -k1 * k1 * z         , 
  -2 - k1 * k1 * z * z ,
  -k1 * k1 * z * z * z )

( -k2 * k2 * z         , 
  -2 - k2 * k2 * z * z ,
  -k2 * k2 * z * z * z )

# Exact solution
coefficient E_ex 
( z         , 
  z * z     ,
  z * z * z )

( z         , 
  z * z     ,
  z * z * z )


# Coefficients of integrators

#coefficient g
#(0, 0, 0),                              # periodic boundaries
#(0, 0, 0),                              # {z=0} boundary
#( I*k2*x*(1-x)*(y*y/2-y*y*y/3)*(1/6), I*k2*y*(1-y)*(x*x/2-x*x*x/3)*(1/6), 0 ), # {z=1} boundary
#(0, 0, 0)                               # surface between materials

coefficient igxn 
( 0          , 0     , 0 ),
( 0          , I     , 0 ),
( -2 * I - k2, I + k2, 0 ),
( 0          , 0     , 0)

coefficient ikbargt                
( 0                      , 0                          , 0 ),
( I * k1bar              , 0                          , 0 ),
( -I * k2bar - k2bar * k2, -2 * I * k2bar - k2bar * k2, 0 ),
(0                       , 0                          , 0 )          

constant one   = 1.0
constant minus = -1.0
coefficient minusksq (-k1*k1), (-k2*k2)
coefficient ii (I), (I)

# Coefficients of boundary integrators
coefficient minusksqbdry (0.0), (-k1*k1bar), (-k2*k2bar), (0.0)
coefficient minusbdry (0.0), (-1.0), (-1.0), (0.0)
coefficient kbdry (0.0), (k1), (k2), (0.0)


# Finite element spaces
fespace fs1 -type=hcurlho           -order=6 -complex        # e, v
                                             -discontinuous
fespace fs2 -type=hcurlho_periodic  -order=3 -complex        # E, F   
fespace fs3 -type=hcurlho_periodic  -order=4 -complex        # M, W
                                             -orderinner=1

fespace fs -type=compound -spaces=[fs1,fs2,fs3] -complex

# Forms

linearform l -fespace=fs
sourceedge f        -comp=1  # (f, v)
neumannedge ikbargt -comp=2  # <i * conj(k) * n x (g x n), F>
neumannedge igxn    -comp=3  # <i * g x n, W>

bilinearform b -fespace=fs
       -nonsym 
       -eliminate_internal
curlcurlpg (2) (1) one                   # (curl E, curl v)
eyeeyeedge (2) (1) minusksq              # - (k*k E, v)
trctrcxn   (3) (1) ii                    # i <<M, v x n>> 
robinedge  minusksqbdry  -comp=2         # - <k*kbar E x n, F x n>
robinedge  minusbdry     -comp=3         # - <M x n, W x n>
xnbdry     (2) (3) kbdry                 # <k E, W x n>
massedge     one -comp=1
curlcurledge one -comp=1

# Solve
gridfunction eEM -fespace=fs

numproc bvp np_solve -bilinearform=b -linearform=l 
        -gridfunction=eEM 
        -solver=direct

#preconditioner c  -type=local -bilinearform=b 

#numproc bvp np_solve -bilinearform=b -linearform=l 
#        -gridfunction=eEM  -preconditioner=c -maxsteps=2000
#        -innerproduct=hermitian    

# Compute error
gridfunction E -fespace=fs2 -addcoef -novisual 
numproc getcomp ng_get -comp=2 -compoundgf=eEM -componentgf=E
coefficient err ( abs(E-E_ex) )
numproc integrate absL2error -coefficient=err



# Visualize
fespace v -h1ho -dim=3 -order=4 -complex
gridfunction E_exact -fespace=v
numproc setvalues np_set -gridfunction=E_exact -coefficient=E_ex

numproc visualization npvis -scalarfunction=eEM.2:1 -subdivision=4 -nolineartexture 



##