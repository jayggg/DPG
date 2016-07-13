from nanogapring import *

p = 1
solution_filename = 'tmp.sol'
mesh_filename = 'tmp.vol.gz'
etot = loadsol(solution_filename,mesh_filename,p)
