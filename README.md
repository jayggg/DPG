# DPG methods in NGSolve

This repository provide facilities to experiment with Discontinuous Petrov Galerkin (DPG) methods, by making  a  shared library that can be added onto 
[NGSolve](http://sourceforge.net/projects/ngsolve/).

CONTRIBUTORS: Jay Gopalakrishnan, Lukas Kogler, Nicole Olivares, Joachim Schoberl.

SEND comments or bug reports to  Jay Gopalakrishnan <gjay@pdx.edu>.


## Install


- Do make sure you have the dependencies installed before proceeding: A working installation of NGSolve and Netgen is required. Please  follow the [instructions online](https://gitlab.asc.tuwien.ac.at/jschoeberl/ngsolve-docu/wikis/home) for installing the development version of these packages. (Please ensure that the compile script `ngscxx`  is in your path after a successful install of NGSolve.) 
- Clone this repository: `git clone https://github.com/jayggg/DPG`
- Navigate to the cloned folder `DPG` and type `make`. This should compile the provided sources on Linux and on Mac (with a working `make`) and create the shared library `libDPG.so`.


## Examples using PDE file interface

Please see comments in individual PDE files in folder `pde`.  Suggestions of pde-file examples to peruse:


1) primallaplace.pde  -  Simplest DPG method for Laplace eq.

2) laplaceadaptive.pde - See how to use adaptivity in the DPG
context. Load the file and press Solve button repeatedly to proceed to
next adaptive iteration.

3) helmholtz.pde - Solve the Helmholtz equation with impedance bc.

4) scatteradaptive.pde - Adaptively compute a scattered wave

5) planewave3d.pde - Compute a 3D wave using DPG method

6) periodicmaxwell.pde - A 3D Maxwell problem with periodic bc in the x- and y-directions


## Other software options for DPG  

If you don't like this code, try out the DPG implementations by our friends at: 

- [Camellia](https://github.com/CamelliaDPG/Camellia)
- [MFEM](https://github.com/mfem/mfem/blob/master/examples/ex8p.cpp)