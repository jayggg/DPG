DPG in NGSolve
==============

This repository makes a shared library that can be added onto 
NGSolve [http://sourceforge.net/projects/ngsolve/].
Its aim is to provide facilities to experiment with Discontinuous
Petrov Galerkin (DPG) methods.

CONTRIBUTORS: Jay Gopalakrishnan, Lukas Kogler, Joachim Schoberl.

SEND comments or bug reports to  gjay@pdx.edu.

---------------------------------------------------------------

SUGGESTION: If you are new to DPG methods, or to this library, 
you may want to start visiting the input pde files (each of which 
is heavily commented) in this order:

1) primallaplace.pde  -  Simplest DPG method for Laplace eq.

2) laplaceadaptive.pde - See how to use adaptivity in the DPG
context. Load the file and press Solve button repeatedly to proceed to
next adaptive iteration.

3) helmholtz1.pde - simplest DPG method for Helmholtz eq w/impedance bc

   helmholtz2.pde - another way to implement impedance bc w/vol elts
   
   helmholtz3.pde - changing formulation to implement impedance bc
   
   helmholtz4.pde - yet another, completely hybridized formulation

4) scatteradaptive.pde - adaptively compute a scattered wave

5) planewave3d - Compute a 3D wave using DPG method


