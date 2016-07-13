""" A demo of the nanogap EOT project on very coarse meshes  """

from nanogapring import *
from ngsolve import *

thickness={'alox'  : 0.010,   # alox layer thickness
           'gold'  : 0.150,   # gold layer thickness
           'glass' : 0.500 }, # glass layer thickness

nlayers={'alox'  : 0,    
         'gold'  : 1,    
         'glass' : 0 }

ncylind=1
w = 0.1


meshfile = 'tmp.vol.gz'  # generate and save the mesh  
m = genmesh(w, nlayers, ncylind, savemeshfile=meshfile)
m.Curve(3)
Draw(m)

print('Solving...')     

# Etot = solve(meshfile=meshfile,localprec=True)
Etot = solve(meshfile=meshfile,localprec=False)

solfile = 'tmp.sol'      # save the computed solution 
Etot.Save(solfile)       # use load_solution.py to view solution later
