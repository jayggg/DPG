# DPG methods in NGSolve

This repository provide facilities to experiment with Discontinuous Petrov Galerkin (DPG) methods, by making  a  shared library that can be added onto 
[NGSolve](http://sourceforge.net/projects/ngsolve/).

CONTRIBUTORS: Jay Gopalakrishnan, Lukas Kogler, Nicole Olivares, Joachim Schoberl.

SEND comments or bug reports to  Jay Gopalakrishnan <gjay@pdx.edu>.


## Install


- Do make sure you have the dependencies installed before proceeding: A working installation of NGSolve and Netgen is required. Please  follow the [instructions online](https://gitlab.asc.tuwien.ac.at/jschoeberl/ngsolve-docu/wikis/home) for installing the development version of these packages. (Please ensure that the compile script `ngscxx`  is in your path after a successful install of NGSolve.) 
- Clone this repository: `git clone https://github.com/jayggg/DPG`
- Navigate to the cloned folder `DPG` and type `make`. This should compile the provided sources on Linux or  Mac systems (where GNU or other `make` is already installed) and should create the shared library called `libDPG`.

## Some examples using python interface

Starting version 6.1, NGsolve provides an interface, called `NGSpy`, to many of its facilities using Python 3, including symbolic forms.  DPG methods can be implemented directly  using these new symbolic facilities, or by loading the the precompiled DPG library from python using CDLL. (The latter is at times faster.)  If you want to explore implementing DPG methods using NGSPy, start with these examples:

- [laplaceadaptive.py](./python/laplaceadaptive.py): In a terminal where PYTHONPATH is set to find the NGsolve libs, navigate to `python` folder and type `netgen  laplaceadaptive.py` to see a demo of automatic adaptivity using DPG methods for the Laplace equation. This example uses pure NGSpy, and there is no need to compile or load `libDPG`.
  
- [periodicmaxwell.py](./python/periodicmaxwell.py): Solve a 3D Maxwell problem, with x and y periodicity, using `libDPG` (which includes an implementation of periodic H(curl) spaces).

## Some examples using PDE file interface

NGSolve uses PDE files (with extension `.pde`) to specify inputs detailing boundary value problems. The DPG library can be utilized using this interface. 
Please see comments in individual `.pde` file examples  within folder `pde`. You may want to explore these pde-file examples:


- [primallaplace.pde](pde/primallaplace.pde)  -  Simplest DPG method for Laplace eq.
- [laplaceadaptive.pde](pde/laplaceadaptive.pde) - See how to use adaptivity in the DPG
context. Load the file and press Solve button repeatedly to proceed to
next adaptive iteration.
- [helmholtz.pde](pde/helmholtz.pde) - Solve the Helmholtz equation with impedance bc.
- [scatteradaptive.pde](pde/scatteradaptive.pde) - Adaptively compute a scattered wave.
- [planewave3d.pde](pde/planewave3d.pde) - Compute a 3D wave using DPG method
- [periodicmaxwell.pde](pde/periodicmaxwell.pde) - A 3D Maxwell problem with periodic boundary conditions in the x and y-directions.

## Index 

- Adaptivity: [Python example](./python/laplaceadaptive.py), [Pde file example](pde/laplaceadaptive.pde)
- Finite Elements:  [DG trace element](spaces/l2trace.cpp)
- Finite Element Spaces: [Periodic spaces](web/periodic.md), [DG trace space](spaces/l2trace.cpp)
- [Hexahedral mesh elements](web/prismhex.md) 
- [Periodic finite element spaces](web/periodic.md) 
- [Periodic meshes](web/periodic.md) 
- [Preconditioner: Schwarz on vertex patches](misc/vertexschwarz.cpp)
- [Prismatic mesh elements](web/prismhex.md) 
- [Quotient norm approximation by polynomial extension](misc/fluxerr.cpp)
- [Traces of DG spaces](spaces/l2trace.cpp)
- [Thin layers](web/prismhex.md) 


## Other software options for DPG methods 

Check out other DPG implementations by our friends:

- [Camellia](https://github.com/CamelliaDPG/Camellia)
- [MFEM](https://github.com/mfem/mfem/blob/master/examples/ex8p.cpp)