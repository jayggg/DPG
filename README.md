DPG in NGSolve
=============

This repository makes a shared library that can be added onto 
NGSolve [http://sourceforge.net/projects/ngsolve/].
Its aim is to provide facilities to experiment with Discontinuous
Petrov Galerkin (DPG) methods.

CONTRIBUTORS: Jay Gopalakrishnan, Lukas Kogler, Nicole Olivares, Joachim Schoberl.

SEND comments or bug reports to  gjay@pdx.edu.

---------------------------------------------------------------

SUGGESTION: If you are new to DPG methods, or to this library, 
you may want to start visiting the input pde files (each of which 
is heavily commented) in this order:

1) primallaplace.pde  -  Simplest DPG method for Laplace eq.

2) laplaceadaptive.pde - See how to use adaptivity in the DPG
context. Load the file and press Solve button repeatedly to proceed to
next adaptive iteration.

3) helmholtz.pde - Solve the Helmholtz equation with impedance bc.

4) scatteradaptive.pde - Adaptively compute a scattered wave

5) planewave3d.pde - Compute a 3D wave using DPG method

6) periodicmaxwell.pde - A 3D Maxwell problem with periodic bc in the x- and y-directions

To use any of these pde files as input, first type "make" on a Mac or
Linux terminal, then "netgen", and then load the pde file. 
